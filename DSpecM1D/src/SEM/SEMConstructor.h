#ifndef SEM_CONSTRUCTOR_H
#define SEM_CONSTRUCTOR_H

#include "SEM.h"
#include "../InputParametersNew.h"

namespace Full1D {

inline SEM::SEM(const InputParametersNew &paramsNew)
    : SEM(paramsNew.earthModel(), paramsNew.maxstep(), paramsNew.nq(),
          paramsNew.inputParameters().lmax()) {}

template <class model1d>
SEM::SEM(const model1d &inp_model, double maxstep, int NQ, int lmax)
    : m_mesh(inp_model, NQ, 1.0, maxstep, false),
      m_frequencyNorm{1.0 / inp_model.TimeNorm()}, m_lmax{lmax},
      _k2{lmax * (lmax + 1)},
      normint{1.0 / (inp_model.TimeNorm() * frequencynorm) *
              sqrt(inp_model.DensityNorm() / densitynorm)},
      m_lengthNorm{inp_model.LengthNorm()},
      bigg_db{6.67230 * std::pow(10.0, -11.0) /
              inp_model.GravitationalConstant()},
      m_momentNorm{
          inp_model.MassNorm() *
          std::pow(inp_model.LengthNorm() / inp_model.TimeNorm(), 2.0)},
      m_nq{NQ} {

  // std::cout << "Radius etc of planet: " << inp_model.OuterRadius() << " m"
  //           << std::endl;
  m_meshModel = MeshModel(m_mesh, inp_model);

  // -----------------------------------------------------------------------
  // boundary and offset information for the local-to-global map
  m_fsb = m_mesh.FS_Boundaries();
  {
    std::size_t totnum = 0;
    for (int idx = 0; idx < m_mesh.NE(); ++idx) {
      auto tmp = m_vecOffset[idx];
      if (totnum < m_fsb.size()) {
        if (idx == m_fsb[totnum]) {
          tmp += 1;
          totnum += 1;
        }
      }
      m_vecOffset.push_back(tmp);
    }
  }

  // -----------------------------------------------------------------------
  // derivative values of Lagrange polynomials at GLL nodes
  auto q = m_mesh.GLL();
  auto pleg =
      Interpolation::LagrangePolynomial(q.Points().begin(), q.Points().end());
  {
    for (int idxk = 0; idxk < q.N(); ++idxk) {
      std::vector<double> vec_tmp(q.N(), 0.0), vec_tmp1(q.N(), 0.0);
      for (int idxi = 0; idxi < q.N(); ++idxi) {
        vec_tmp[idxi] = pleg.Derivative(idxi, q.X(idxk));
        vec_tmp1[idxi] = pleg(idxi, q.X(idxk));
      }
      vec_lag_deriv.push_back(vec_tmp);
      vec_delta.push_back(vec_tmp1);
    };
  }

  // -----------------------------------------------------------------------
  // radial matrices
  {
    using T = Eigen::Triplet<double>;

    totlen = this->ltgS(2, m_mesh.NE() - 1, NQ - 1) + 1;

    // fluid / solid bookkeeping
    m_vecFluid = std::vector<int>(m_mesh.NE(), 0);
    m_vecDof = std::vector<bool>(m_mesh.NE(), false);
    for (int idxe = 0; idxe < m_mesh.NE(); ++idxe) {
      if (inp_model.IsFluid(m_mesh.LayerNumber(idxe))) {
        m_vecFluid[idxe] = 1;
        m_hasFluid = true;
      }
      if (idxe > 0 && (std::abs(m_vecFluid[idxe] - m_vecFluid[idxe - 1]) == 1))
        m_vecDof[idxe - 1] = true;
    }

    // find m_el and m_eu — toroidal DOF range
    // Toroidal modes live in the mantle: the solid region immediately above
    // the deepest fluid region (outer core). Walk from the centre outward to
    // find the deepest fluid region, ignoring shallower fluid layers (ocean).
    // For a purely solid model: m_el=0, m_eu=NE.
    m_el = 0;
    m_eu = m_mesh.NE();
    if (m_hasFluid) {
      // Walk from centre (idx=0) upward to find the first (deepest) fluid
      // region — this is the outer core for Earth models.
      int first_fluid_from_bottom = -1;
      for (int idx = 0; idx < m_mesh.NE(); ++idx) {
        if (m_vecFluid[idx] == 1) {
          first_fluid_from_bottom = idx;
          break;
        }
      }

      if (first_fluid_from_bottom >= 0) {
        // Find the top of this contiguous fluid region (outer core top)
        int top_of_deepest_fluid = first_fluid_from_bottom;
        for (int idx = first_fluid_from_bottom; idx < m_mesh.NE(); ++idx) {
          if (m_vecFluid[idx] == 1)
            top_of_deepest_fluid = idx;
          else
            break;
        }
        // m_el: first solid element above the outer core (mantle base)
        m_el = top_of_deepest_fluid + 1;
        // m_eu: next fluid element above m_el (e.g. ocean), or NE if none
        m_eu = m_mesh.NE();
        for (int idx = m_el; idx < m_mesh.NE(); ++idx) {
          if (m_vecFluid[idx] == 1) {
            m_eu = idx;
            break;
          }
        }
      }
    }
    // std::cout << "\nToroidal range: " << m_el << " to " << m_eu - 1
    //           << ", total NE: " << m_mesh.NE() << std::endl;
    // std::cout << "\nRadius range for toroidal modes: " << m_mesh.ELR(m_el)
    //           << " to " << m_mesh.EUR(m_eu - 1) << std::endl;
    // std::cout << "\nRadius range of uppermost fluid layer: "
    //           << m_mesh.ELR(m_mesh.NE() - 1) << " to "
    //           << m_mesh.EUR(m_mesh.NE() - 1) << std::endl;

    // radial inertia + stiffness matrices
    {
      auto num_radial_dof = this->ltgR(1, m_mesh.NE() - 1, NQ - 1) + 1;
      std::vector<T> tpl_in, tpl_ke, tpl_ke_atten;

      for (int idxe = 0; idxe < m_mesh.NE(); ++idxe) {
        double elem_width = m_mesh.EW(idxe);
        for (int i = 0; i < q.N(); ++i) {
          double xrad = m_mesh.NodeRadius(idxe, i);
          double tmp = elem_width / 2.0 * q.W(i) *
                       m_meshModel.Density(idxe, i) * xrad * xrad;
          auto idx_uu = this->ltgR(0, idxe, i);
          tpl_in.push_back(T(idx_uu, idx_uu, tmp));
        }
      }
      m_matInertia0.resize(num_radial_dof, num_radial_dof);
      m_matInertia0.setFromTriplets(tpl_in.begin(), tpl_in.end());
      m_matInertia0.makeCompressed();

      // diagonal stiffness terms
      for (int idxe = 0; idxe < m_mesh.NE(); ++idxe) {
        double elem_width = m_mesh.EW(idxe);
        for (int i = 0; i < q.N(); ++i) {
          double tmp0 = elem_width / 2.0 * q.W(i);
          auto crad = m_mesh.NodeRadius(idxe, i);
          auto crho = m_meshModel.Density(idxe, i);
          auto gi = m_meshModel.Gravity(idxe, i);
          auto Ai = m_meshModel.A(idxe, i);
          auto Ni = m_meshModel.N(idxe, i);
          auto Ai_a = m_meshModel.A_atten(idxe, i);
          auto Ni_a = m_meshModel.N_atten(idxe, i);
          double tmp_u =
              tmp0 * (4.0 * crho * (pi_db * bigg_db * crho * crad - gi) * crad +
                      4 * (Ai - Ni));
          double tmp_u_a = 4 * tmp0 * (Ai_a - Ni_a);
          auto idxtiu = this->ltgR(0, idxe, i);
          tpl_ke.push_back(T(idxtiu, idxtiu, tmp_u));
          tpl_ke_atten.push_back(T(idxtiu, idxtiu, tmp_u_a));
        }
      }
      // surface gravity boundary term
      int idxpb = this->ltgR(1, m_mesh.NE() - 1, NQ - 1);
      double rpb = m_mesh.NodeRadius(m_mesh.NE() - 1, NQ - 1);
      tpl_ke.push_back(T(idxpb, idxpb, rpb / (4.0 * pi_db * bigg_db)));

      // coupling stiffness terms
      for (int idxe = 0; idxe < m_mesh.NE(); ++idxe) {
        double elem_width = m_mesh.EW(idxe);
        auto e2 = 0.5 * elem_width;
        double d_val = 2.0 / elem_width;
        stdvvec mat_d(q.N(), stdvec(q.N(), 0.0));
        for (int i = 0; i < q.N(); ++i)
          for (int j = 0; j < q.N(); ++j)
            mat_d[i][j] = d_val * vec_lag_deriv[j][i];

        std::vector<double> FV(q.N()), FVA(q.N());
        for (int i = 0; i < q.N(); ++i) {
          FV[i] = m_meshModel.F(idxe, i);
          FVA[i] = m_meshModel.F_atten(idxe, i);
        }

        for (int i = 0; i < q.N(); ++i) {
          for (int j = 0; j < q.N(); ++j) {
            auto rj = m_mesh.NodeRadius(idxe, j);
            auto rhoj = m_meshModel.Density(idxe, j);
            auto tmp_uu = elem_width * q.W(j) * FV[j] * rj * mat_d[i][j];
            auto tmp_pdu = q.W(j) * rhoj * rj * rj * vec_lag_deriv[j][i];
            auto tmp_uu_a = elem_width * q.W(j) * FVA[j] * rj * mat_d[i][j];
            auto idx_u_i = this->ltgR(0, idxe, i);
            auto idx_u_j = this->ltgR(0, idxe, j);
            auto idx_p_i = this->ltgR(1, idxe, i);
            tpl_ke.push_back(T(idx_u_i, idx_u_j, tmp_uu));
            tpl_ke.push_back(T(idx_u_j, idx_u_i, tmp_uu));
            tpl_ke.push_back(T(idx_p_i, idx_u_j, tmp_pdu));
            tpl_ke.push_back(T(idx_u_j, idx_p_i, tmp_pdu));
            tpl_ke_atten.push_back(T(idx_u_i, idx_u_j, tmp_uu_a));
            tpl_ke_atten.push_back(T(idx_u_j, idx_u_i, tmp_uu_a));
          }
        }

        for (int i = 0; i < q.N(); ++i) {
          for (int j = 0; j < q.N(); ++j) {
            double tmp_uu = 0.0, tmp_pp = 0.0, tmp_uu_a = 0.0;
            for (int k = 0; k < q.N(); ++k) {
              auto rk = m_mesh.NodeRadius(idxe, k);
              auto ddrr = rk * rk * mat_d[i][k] * mat_d[j][k];
              tmp_uu += q.W(k) * m_meshModel.C(idxe, k) * ddrr;
              tmp_pp += q.W(k) * ddrr;
              tmp_uu_a += q.W(k) * m_meshModel.C_atten(idxe, k) * ddrr;
            }
            tmp_uu *= e2;
            tmp_pp *= e2 / (4.0 * pi_db * bigg_db);
            tmp_uu_a *= e2;
            auto idx_u_i = this->ltgR(0, idxe, i);
            auto idx_u_j = this->ltgR(0, idxe, j);
            auto idx_p_i = this->ltgR(1, idxe, i);
            auto idx_p_j = this->ltgR(1, idxe, j);
            tpl_ke.push_back(T(idx_u_i, idx_u_j, tmp_uu));
            tpl_ke.push_back(T(idx_p_i, idx_p_j, tmp_pp));
            tpl_ke_atten.push_back(T(idx_u_i, idx_u_j, tmp_uu_a));
          }
        }
      }

      m_matKe0.resize(num_radial_dof, num_radial_dof);
      m_matKe0.setFromTriplets(tpl_ke.begin(), tpl_ke.end());
      m_matKe0.makeCompressed();
      m_matKe0Atten.resize(num_radial_dof, num_radial_dof);
      m_matKe0Atten.setFromTriplets(tpl_ke_atten.begin(), tpl_ke_atten.end());
      m_matKe0Atten.makeCompressed();
    }
  }

  // -----------------------------------------------------------------------
  // toroidal base matrices
  {
    using T = Eigen::Triplet<double>;
    std::vector<T> tpl_in_0, tpl_ke_1, tpl_ke_2, tpl_ke_a1, tpl_ke_a2;
    auto ntdof = this->ltgT(m_eu - 1, NQ - 1) + 1;
    // std::cout << "Lower element: " << m_el << ", upper element: " << m_eu
    //           << std::endl;
    // std::cout << "Total toroidal DOFs: " << ntdof << std::endl;
    // std::cout << "Correct DOF range: " << this->ltgT(m_el, 0) << " to "
    // << this->ltgT(m_eu - 1, NQ - 1) << std::endl;
    for (int idxe = m_el; idxe < m_eu; ++idxe) {
      double elem_width = m_mesh.EW(idxe);
      for (int i = 0; i < q.N(); ++i) {
        double xrad = m_mesh.NodeRadius(idxe, i);
        double tmp = elem_width / 2.0 * q.W(i) * m_meshModel.Density(idxe, i) *
                     xrad * xrad;
        auto idx_ww = this->ltgT(idxe, i);
        tpl_in_0.push_back(T(idx_ww, idx_ww, tmp));
      }
    }
    m_matInTBase.resize(ntdof, ntdof);
    m_matInTBase.setFromTriplets(tpl_in_0.begin(), tpl_in_0.end());
    m_matInTBase.makeCompressed();

    for (int idxe = m_el; idxe < m_eu; ++idxe) {
      double elem_width = m_mesh.EW(idxe);
      double e2 = 0.5 * elem_width;
      double d_val = 2.0 / elem_width;
      stdvvec mat_d(q.N(), stdvec(q.N(), 0.0));
      for (int i = 0; i < q.N(); ++i)
        for (int j = 0; j < q.N(); ++j)
          mat_d[i][j] = d_val * vec_lag_deriv[j][i];

      std::vector<double> LV(q.N()), NV(q.N()), RR(q.N()), LVA(q.N()),
          NVA(q.N());
      std::vector<int> IDXW(q.N());
      for (int i = 0; i < q.N(); ++i) {
        LV[i] = m_meshModel.L(idxe, i);
        NV[i] = m_meshModel.N(idxe, i);
        RR[i] = m_mesh.NodeRadius(idxe, i);
        LVA[i] = m_meshModel.L_atten(idxe, i);
        NVA[i] = m_meshModel.N_atten(idxe, i);
        IDXW[i] = this->ltgT(idxe, i);
      }

      for (int i = 0; i < q.N(); ++i) {
        double tmp0 = elem_width / 2.0 * q.W(i);
        tpl_ke_1.push_back(T(IDXW[i], IDXW[i], tmp0 * (LV[i] - 2.0 * NV[i])));
        tpl_ke_2.push_back(T(IDXW[i], IDXW[i], tmp0 * NV[i]));
        tpl_ke_a1.push_back(
            T(IDXW[i], IDXW[i], tmp0 * (LVA[i] - 2.0 * NVA[i])));
        tpl_ke_a2.push_back(T(IDXW[i], IDXW[i], tmp0 * NVA[i]));
      }
      for (int i = 0; i < q.N(); ++i) {
        for (int j = 0; j < q.N(); ++j) {
          auto t = -e2 * q.W(j) * LV[j] * RR[j] * mat_d[i][j];
          auto ta = -e2 * q.W(j) * LVA[j] * RR[j] * mat_d[i][j];
          tpl_ke_1.push_back(T(IDXW[i], IDXW[j], t));
          tpl_ke_1.push_back(T(IDXW[j], IDXW[i], t));
          tpl_ke_a1.push_back(T(IDXW[i], IDXW[j], ta));
          tpl_ke_a1.push_back(T(IDXW[j], IDXW[i], ta));
        }
      }
      for (int i = 0; i < q.N(); ++i) {
        for (int j = 0; j < q.N(); ++j) {
          double tmp_ww = 0.0, tmp_ww_a = 0.0;
          for (int k = 0; k < q.N(); ++k) {
            auto ddrr = RR[k] * RR[k] * mat_d[i][k] * mat_d[j][k];
            tmp_ww += q.W(k) * LV[k] * ddrr;
            tmp_ww_a += q.W(k) * LVA[k] * ddrr;
          }
          tpl_ke_1.push_back(T(IDXW[i], IDXW[j], e2 * tmp_ww));
          tpl_ke_a1.push_back(T(IDXW[i], IDXW[j], e2 * tmp_ww_a));
        }
      }
    }

    auto make_smat = [&ntdof](std::vector<Eigen::Triplet<double>> &tpl) {
      Eigen::SparseMatrix<double> m(ntdof, ntdof);
      m.setFromTriplets(tpl.begin(), tpl.end());
      m.makeCompressed();
      return m;
    };
    m_vecKeTBase.push_back(make_smat(tpl_ke_1));
    m_vecKeTBase.push_back(make_smat(tpl_ke_2));
    m_vecKeTAtten.push_back(make_smat(tpl_ke_a1));
    m_vecKeTAtten.push_back(make_smat(tpl_ke_a2));
  }

  // -----------------------------------------------------------------------
  // spheroidal base matrices
  {
    m_vecKeSBase.reserve(3);
    m_vecKeSAtten.reserve(3);
    m_vecInSBase.reserve(2);
    auto totlen_S = this->ltgS(2, m_mesh.NE() - 1, NQ - 1) + 1;
    using T = Eigen::Triplet<double>;
    std::vector<T> tpl_in_0, tpl_in_1;
    std::vector<T> tpl_ke_0, tpl_ke_1, tpl_ke_2;
    std::vector<T> tpl_ke_a0, tpl_ke_a1, tpl_ke_a2;

    // inertia
    for (int idxe = 0; idxe < m_mesh.NE(); ++idxe) {
      for (int i = 0; i < q.N(); ++i) {
        double xrad = m_mesh.NodeRadius(idxe, i);
        double tmp = m_mesh.EW(idxe) / 2.0 * q.W(i) *
                     m_meshModel.Density(idxe, i) * xrad * xrad;
        auto idx_uu = this->ltgS(0, idxe, i);
        auto idx_vv = this->ltgS(1, idxe, i);
        tpl_in_0.push_back(T(idx_uu, idx_uu, tmp));
        tpl_in_1.push_back(T(idx_vv, idx_vv, tmp));
      }
    }

    // diagonal stiffness
    for (int idxe = 0; idxe < m_mesh.NE(); ++idxe) {
      double elem_width = m_mesh.EW(idxe);
      for (int i = 0; i < q.N(); ++i) {
        double tmp0 = elem_width / 2.0 * q.W(i);
        auto Li = m_meshModel.L(idxe, i);
        auto Ai = m_meshModel.A(idxe, i);
        auto Ni = m_meshModel.N(idxe, i);
        auto gi = m_meshModel.Gravity(idxe, i);
        auto crad = m_mesh.NodeRadius(idxe, i);
        auto crho = m_meshModel.Density(idxe, i);
        auto Li_a = m_meshModel.L_atten(idxe, i);
        auto Ai_a = m_meshModel.A_atten(idxe, i);
        auto Ni_a = m_meshModel.N_atten(idxe, i);
        auto idxtiu = this->ltgS(0, idxe, i);
        auto idxtiv = this->ltgS(1, idxe, i);
        auto idxtip = this->ltgS(2, idxe, i);

        tpl_ke_0.push_back(
            T(idxtiu, idxtiu,
              4.0 * tmp0 *
                  (crho * (pi_db * bigg_db * crho * crad - gi) * crad + Ai -
                   Ni)));
        tpl_ke_1.push_back(T(idxtiu, idxtiu, tmp0 * Li));
        tpl_ke_1.push_back(T(idxtiv, idxtiv, tmp0 * (Li - 2 * Ni)));
        tpl_ke_2.push_back(T(idxtiv, idxtiv, tmp0 * Ai));
        tpl_ke_1.push_back(T(idxtip, idxtip, tmp0 / (4.0 * pi_db * bigg_db)));
        tpl_ke_1.push_back(T(idxtip, idxtiv, tmp0 * crho * crad));
        tpl_ke_1.push_back(T(idxtiv, idxtip, tmp0 * crho * crad));
        tpl_ke_1.push_back(
            T(idxtiu, idxtiv, tmp0 * (crho * gi * crad - Li - 2 * (Ai - Ni))));
        tpl_ke_1.push_back(
            T(idxtiv, idxtiu, tmp0 * (crho * gi * crad - Li - 2 * (Ai - Ni))));

        tpl_ke_a0.push_back(T(idxtiu, idxtiu, 4.0 * tmp0 * (Ai_a - Ni_a)));
        tpl_ke_a1.push_back(T(idxtiu, idxtiu, tmp0 * Li_a));
        tpl_ke_a1.push_back(T(idxtiv, idxtiv, tmp0 * (Li_a - 2 * Ni_a)));
        tpl_ke_a2.push_back(T(idxtiv, idxtiv, tmp0 * Ai_a));
        tpl_ke_a1.push_back(
            T(idxtiu, idxtiv, -tmp0 * (Li_a + 2 * (Ai_a - Ni_a))));
        tpl_ke_a1.push_back(
            T(idxtiv, idxtiu, -tmp0 * (Li_a + 2 * (Ai_a - Ni_a))));
      }
    }

    // coupling stiffness
    for (int idxe = 0; idxe < m_mesh.NE(); ++idxe) {
      double elem_width = m_mesh.EW(idxe);
      double e2 = 0.5 * elem_width;
      double d_val = 2.0 / elem_width;
      stdvvec mat_d(q.N(), stdvec(q.N(), 0.0));
      for (int i = 0; i < q.N(); ++i)
        for (int j = 0; j < q.N(); ++j)
          mat_d[i][j] = d_val * vec_lag_deriv[j][i];

      std::vector<double> LV(q.N()), NV(q.N()), AV(q.N()), FV(q.N()), RR(q.N());
      std::vector<double> LVA(q.N()), NVA(q.N()), AVA(q.N()), FVA(q.N());
      for (int i = 0; i < q.N(); ++i) {
        LV[i] = m_meshModel.L(idxe, i);
        NV[i] = m_meshModel.N(idxe, i);
        AV[i] = m_meshModel.A(idxe, i);
        FV[i] = m_meshModel.F(idxe, i);
        RR[i] = m_mesh.NodeRadius(idxe, i);
        LVA[i] = m_meshModel.L_atten(idxe, i);
        NVA[i] = m_meshModel.N_atten(idxe, i);
        AVA[i] = m_meshModel.A_atten(idxe, i);
        FVA[i] = m_meshModel.F_atten(idxe, i);
      }

      for (int i = 0; i < q.N(); ++i) {
        auto ri = m_mesh.NodeRadius(idxe, i);
        auto rhoi = m_meshModel.Density(idxe, i);
        for (int j = 0; j < q.N(); ++j) {
          auto rj = m_mesh.NodeRadius(idxe, j);
          auto rhoj = m_meshModel.Density(idxe, j);
          auto tmp_uu_0 = elem_width * q.W(j) * FV[j] * rj * mat_d[i][j];
          auto tmp_vv_1 = -e2 * q.W(j) * LV[j] * rj * mat_d[i][j];
          auto tmp_uvd_1 = e2 * q.W(i) * LV[i] * ri * mat_d[j][i];
          auto tmp_udv_1 = -e2 * q.W(j) * FV[j] * rj * mat_d[i][j];
          auto tmp_pdu_0 = q.W(j) * rhoj * rj * rj * vec_lag_deriv[j][i];
          auto tmp_uu_0a = elem_width * q.W(j) * FVA[j] * rj * mat_d[i][j];
          auto tmp_vv_1a = -e2 * q.W(j) * LVA[j] * rj * mat_d[i][j];
          auto tmp_uvd_1a = e2 * q.W(i) * LVA[i] * ri * mat_d[j][i];
          auto tmp_udv_1a = -e2 * q.W(j) * FVA[j] * rj * mat_d[i][j];
          auto idx_u_i = this->ltgS(0, idxe, i),
               idx_u_j = this->ltgS(0, idxe, j);
          auto idx_v_i = this->ltgS(1, idxe, i),
               idx_v_j = this->ltgS(1, idxe, j);
          auto idx_p_i = this->ltgS(2, idxe, i);
          tpl_ke_0.push_back(T(idx_u_i, idx_u_j, tmp_uu_0));
          tpl_ke_0.push_back(T(idx_u_j, idx_u_i, tmp_uu_0));
          tpl_ke_1.push_back(T(idx_v_i, idx_v_j, tmp_vv_1));
          tpl_ke_1.push_back(T(idx_v_j, idx_v_i, tmp_vv_1));
          tpl_ke_1.push_back(T(idx_u_i, idx_v_j, tmp_uvd_1));
          tpl_ke_1.push_back(T(idx_v_j, idx_u_i, tmp_uvd_1));
          tpl_ke_1.push_back(T(idx_u_i, idx_v_j, tmp_udv_1));
          tpl_ke_1.push_back(T(idx_v_j, idx_u_i, tmp_udv_1));
          tpl_ke_0.push_back(T(idx_p_i, idx_u_j, tmp_pdu_0));
          tpl_ke_0.push_back(T(idx_u_j, idx_p_i, tmp_pdu_0));
          tpl_ke_a0.push_back(T(idx_u_i, idx_u_j, tmp_uu_0a));
          tpl_ke_a0.push_back(T(idx_u_j, idx_u_i, tmp_uu_0a));
          tpl_ke_a1.push_back(T(idx_v_i, idx_v_j, tmp_vv_1a));
          tpl_ke_a1.push_back(T(idx_v_j, idx_v_i, tmp_vv_1a));
          tpl_ke_a1.push_back(T(idx_u_i, idx_v_j, tmp_uvd_1a));
          tpl_ke_a1.push_back(T(idx_v_j, idx_u_i, tmp_uvd_1a));
          tpl_ke_a1.push_back(T(idx_u_i, idx_v_j, tmp_udv_1a));
          tpl_ke_a1.push_back(T(idx_v_j, idx_u_i, tmp_udv_1a));
        }
      }

      for (int i = 0; i < q.N(); ++i) {
        for (int j = 0; j < q.N(); ++j) {
          double tmp_uu_0 = 0.0, tmp_vv_1 = 0.0, tmp_pp_0 = 0.0;
          double tmp_uu_0a = 0.0, tmp_vv_1a = 0.0;
          for (int k = 0; k < q.N(); ++k) {
            auto rk = m_mesh.NodeRadius(idxe, k);
            auto ddrr = rk * rk * mat_d[i][k] * mat_d[j][k];
            tmp_uu_0 += q.W(k) * m_meshModel.C(idxe, k) * ddrr;
            tmp_vv_1 += q.W(k) * m_meshModel.L(idxe, k) * ddrr;
            tmp_pp_0 += q.W(k) * ddrr;
            tmp_uu_0a += q.W(k) * m_meshModel.C_atten(idxe, k) * ddrr;
            tmp_vv_1a += q.W(k) * m_meshModel.L_atten(idxe, k) * ddrr;
          }
          auto idx_u_i = this->ltgS(0, idxe, i),
               idx_u_j = this->ltgS(0, idxe, j);
          auto idx_v_i = this->ltgS(1, idxe, i),
               idx_v_j = this->ltgS(1, idxe, j);
          auto idx_p_i = this->ltgS(2, idxe, i),
               idx_p_j = this->ltgS(2, idxe, j);
          tpl_ke_0.push_back(T(idx_u_i, idx_u_j, e2 * tmp_uu_0));
          tpl_ke_1.push_back(T(idx_v_i, idx_v_j, e2 * tmp_vv_1));
          tpl_ke_0.push_back(
              T(idx_p_i, idx_p_j, e2 * tmp_pp_0 / (4.0 * pi_db * bigg_db)));
          tpl_ke_a0.push_back(T(idx_u_i, idx_u_j, e2 * tmp_uu_0a));
          tpl_ke_a1.push_back(T(idx_v_i, idx_v_j, e2 * tmp_vv_1a));
        }
      }
    }

    auto make_smat = [&totlen_S](std::vector<Eigen::Triplet<double>> &tpl) {
      Eigen::SparseMatrix<double> m(totlen_S, totlen_S);
      m.setFromTriplets(tpl.begin(), tpl.end());
      m.makeCompressed();
      return m;
    };
    auto make_in = [&totlen_S](std::vector<Eigen::Triplet<double>> &tpl) {
      Eigen::SparseMatrix<double> m(totlen_S, totlen_S);
      m.setFromTriplets(tpl.begin(), tpl.end());
      m.makeCompressed();
      return m;
    };
    m_vecInSBase.push_back(make_in(tpl_in_0));
    m_vecInSBase.push_back(make_in(tpl_in_1));
    m_vecKeSBase.push_back(make_smat(tpl_ke_0));
    m_vecKeSBase.push_back(make_smat(tpl_ke_1));
    m_vecKeSBase.push_back(make_smat(tpl_ke_2));
    m_vecKeSAtten.push_back(make_smat(tpl_ke_a0));
    m_vecKeSAtten.push_back(make_smat(tpl_ke_a1));
    m_vecKeSAtten.push_back(make_smat(tpl_ke_a2));
  }
};

}   // namespace Full1D

#endif   // SEM_CONSTRUCTOR_H
