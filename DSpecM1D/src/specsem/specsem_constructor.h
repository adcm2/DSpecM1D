#ifndef SPECSEM_CONSTRUCTOR_H
#define SPECSEM_CONSTRUCTOR_H

#include "specsem.h"

namespace Full1D {

template <class model1d>
specsem::specsem(const model1d &inp_model, double maxstep, int NQ, int lmax)
    : _mesh(inp_model, NQ, 1.0, maxstep, false),
      _freq_norm{1.0 / inp_model.TimeNorm()}, _lmax{lmax},
      _k2{lmax * (lmax + 1)},
      normint{1.0 / (inp_model.TimeNorm() * frequencynorm) *
              sqrt(inp_model.DensityNorm() / densitynorm)},
      _length_norm{inp_model.LengthNorm()},
      bigg_db{6.67230 * std::pow(10.0, -11.0) /
              inp_model.GravitationalConstant()},
      _moment_norm{
          inp_model.MassNorm() *
          std::pow(inp_model.LengthNorm() / inp_model.TimeNorm(), 2.0)},
      _NQ{NQ} {

  _mesh_model = MeshModel(_mesh, inp_model);

  // -----------------------------------------------------------------------
  // boundary and offset information for the local-to-global map
  fsb = _mesh.FS_Boundaries();
  {
    std::size_t totnum = 0;
    for (int idx = 0; idx < _mesh.NE(); ++idx) {
      auto tmp = vec_offset[idx];
      if (totnum < fsb.size()) {
        if (idx == fsb[totnum]) {
          tmp += 1;
          totnum += 1;
        }
      }
      vec_offset.push_back(tmp);
    }
  }

  // -----------------------------------------------------------------------
  // derivative values of Lagrange polynomials at GLL nodes
  auto q = _mesh.GLL();
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

    totlen = this->LtG_S(2, _mesh.NE() - 1, NQ - 1) + 1;

    // fluid / solid bookkeeping
    _vec_fluid = std::vector<int>(_mesh.NE(), 0);
    _vec_dof = std::vector<bool>(_mesh.NE(), false);
    for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
      if (inp_model.IsFluid(_mesh.LayerNumber(idxe))) {
        _vec_fluid[idxe] = 1;
        _has_fluid = true;
      }
      if (idxe > 0 && (std::abs(_vec_fluid[idxe] - _vec_fluid[idxe - 1]) == 1))
        _vec_dof[idxe - 1] = true;
    }

    // find _el and _eu — toroidal DOF range
    // Toroidal modes live in the mantle: the solid region immediately above
    // the deepest fluid region (outer core). Walk from the centre outward to
    // find the deepest fluid region, ignoring shallower fluid layers (ocean).
    // For a purely solid model: _el=0, _eu=NE.
    _el = 0;
    _eu = _mesh.NE();
    if (_has_fluid) {
      // Walk from centre (idx=0) upward to find the first (deepest) fluid
      // region — this is the outer core for Earth models.
      int first_fluid_from_bottom = -1;
      for (int idx = 0; idx < _mesh.NE(); ++idx) {
        if (_vec_fluid[idx] == 1) {
          first_fluid_from_bottom = idx;
          break;
        }
      }

      if (first_fluid_from_bottom >= 0) {
        // Find the top of this contiguous fluid region (outer core top)
        int top_of_deepest_fluid = first_fluid_from_bottom;
        for (int idx = first_fluid_from_bottom; idx < _mesh.NE(); ++idx) {
          if (_vec_fluid[idx] == 1)
            top_of_deepest_fluid = idx;
          else
            break;
        }
        // _el: first solid element above the outer core (mantle base)
        _el = top_of_deepest_fluid + 1;
        // _eu: next fluid element above _el (e.g. ocean), or NE if none
        _eu = _mesh.NE();
        for (int idx = _el; idx < _mesh.NE(); ++idx) {
          if (_vec_fluid[idx] == 1) {
            _eu = idx;
            break;
          }
        }
      }
    }

    // radial inertia + stiffness matrices
    {
      auto num_radial_dof = this->LtG_R(1, _mesh.NE() - 1, NQ - 1) + 1;
      std::vector<T> tpl_in, tpl_ke, tpl_ke_atten;

      for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
        double elem_width = _mesh.EW(idxe);
        for (int i = 0; i < q.N(); ++i) {
          double xrad = _mesh.NodeRadius(idxe, i);
          double tmp = elem_width / 2.0 * q.W(i) *
                       _mesh_model.Density(idxe, i) * xrad * xrad;
          auto idx_uu = this->LtG_R(0, idxe, i);
          tpl_in.push_back(T(idx_uu, idx_uu, tmp));
        }
      }
      mat_inertia_0.resize(num_radial_dof, num_radial_dof);
      mat_inertia_0.setFromTriplets(tpl_in.begin(), tpl_in.end());
      mat_inertia_0.makeCompressed();

      // diagonal stiffness terms
      for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
        double elem_width = _mesh.EW(idxe);
        for (int i = 0; i < q.N(); ++i) {
          double tmp0 = elem_width / 2.0 * q.W(i);
          auto crad = _mesh.NodeRadius(idxe, i);
          auto crho = _mesh_model.Density(idxe, i);
          auto gi = _mesh_model.Gravity(idxe, i);
          auto Ai = _mesh_model.A(idxe, i);
          auto Ni = _mesh_model.N(idxe, i);
          auto Ai_a = _mesh_model.A_atten(idxe, i);
          auto Ni_a = _mesh_model.N_atten(idxe, i);
          double tmp_u =
              tmp0 * (4.0 * crho * (pi_db * bigg_db * crho * crad - gi) * crad +
                      4 * (Ai - Ni));
          double tmp_u_a = 4 * tmp0 * (Ai_a - Ni_a);
          auto idxtiu = this->LtG_R(0, idxe, i);
          tpl_ke.push_back(T(idxtiu, idxtiu, tmp_u));
          tpl_ke_atten.push_back(T(idxtiu, idxtiu, tmp_u_a));
        }
      }
      // surface gravity boundary term
      int idxpb = this->LtG_R(1, _mesh.NE() - 1, NQ - 1);
      double rpb = _mesh.NodeRadius(_mesh.NE() - 1, NQ - 1);
      tpl_ke.push_back(T(idxpb, idxpb, rpb / (4.0 * pi_db * bigg_db)));

      // coupling stiffness terms
      for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
        double elem_width = _mesh.EW(idxe);
        auto e2 = 0.5 * elem_width;
        double d_val = 2.0 / elem_width;
        stdvvec mat_d(q.N(), stdvec(q.N(), 0.0));
        for (int i = 0; i < q.N(); ++i)
          for (int j = 0; j < q.N(); ++j)
            mat_d[i][j] = d_val * vec_lag_deriv[j][i];

        std::vector<double> FV(q.N()), FVA(q.N());
        for (int i = 0; i < q.N(); ++i) {
          FV[i] = _mesh_model.F(idxe, i);
          FVA[i] = _mesh_model.F_atten(idxe, i);
        }

        for (int i = 0; i < q.N(); ++i) {
          for (int j = 0; j < q.N(); ++j) {
            auto rj = _mesh.NodeRadius(idxe, j);
            auto rhoj = _mesh_model.Density(idxe, j);
            auto tmp_uu = elem_width * q.W(j) * FV[j] * rj * mat_d[i][j];
            auto tmp_pdu = q.W(j) * rhoj * rj * rj * vec_lag_deriv[j][i];
            auto tmp_uu_a = elem_width * q.W(j) * FVA[j] * rj * mat_d[i][j];
            auto idx_u_i = this->LtG_R(0, idxe, i);
            auto idx_u_j = this->LtG_R(0, idxe, j);
            auto idx_p_i = this->LtG_R(1, idxe, i);
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
              auto rk = _mesh.NodeRadius(idxe, k);
              auto ddrr = rk * rk * mat_d[i][k] * mat_d[j][k];
              tmp_uu += q.W(k) * _mesh_model.C(idxe, k) * ddrr;
              tmp_pp += q.W(k) * ddrr;
              tmp_uu_a += q.W(k) * _mesh_model.C_atten(idxe, k) * ddrr;
            }
            tmp_uu *= e2;
            tmp_pp *= e2 / (4.0 * pi_db * bigg_db);
            tmp_uu_a *= e2;
            auto idx_u_i = this->LtG_R(0, idxe, i);
            auto idx_u_j = this->LtG_R(0, idxe, j);
            auto idx_p_i = this->LtG_R(1, idxe, i);
            auto idx_p_j = this->LtG_R(1, idxe, j);
            tpl_ke.push_back(T(idx_u_i, idx_u_j, tmp_uu));
            tpl_ke.push_back(T(idx_p_i, idx_p_j, tmp_pp));
            tpl_ke_atten.push_back(T(idx_u_i, idx_u_j, tmp_uu_a));
          }
        }
      }

      mat_ke_0.resize(num_radial_dof, num_radial_dof);
      mat_ke_0.setFromTriplets(tpl_ke.begin(), tpl_ke.end());
      mat_ke_0.makeCompressed();
      mat_ke_0_atten.resize(num_radial_dof, num_radial_dof);
      mat_ke_0_atten.setFromTriplets(tpl_ke_atten.begin(), tpl_ke_atten.end());
      mat_ke_0_atten.makeCompressed();
    }
  }

  // -----------------------------------------------------------------------
  // toroidal base matrices
  {
    using T = Eigen::Triplet<double>;
    std::vector<T> tpl_in_0, tpl_ke_1, tpl_ke_2, tpl_ke_a1, tpl_ke_a2;
    auto ntdof = this->LtG_T(_mesh.NE() - 1, NQ - 1) + 1;

    for (int idxe = _el; idxe < _eu; ++idxe) {
      double elem_width = _mesh.EW(idxe);
      for (int i = 0; i < q.N(); ++i) {
        double xrad = _mesh.NodeRadius(idxe, i);
        double tmp = elem_width / 2.0 * q.W(i) * _mesh_model.Density(idxe, i) *
                     xrad * xrad;
        auto idx_ww = this->LtG_T(idxe, i);
        tpl_in_0.push_back(T(idx_ww, idx_ww, tmp));
      }
    }
    mat_in_t_base.resize(ntdof, ntdof);
    mat_in_t_base.setFromTriplets(tpl_in_0.begin(), tpl_in_0.end());
    mat_in_t_base.makeCompressed();

    for (int idxe = _el; idxe < _eu; ++idxe) {
      double elem_width = _mesh.EW(idxe);
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
        LV[i] = _mesh_model.L(idxe, i);
        NV[i] = _mesh_model.N(idxe, i);
        RR[i] = _mesh.NodeRadius(idxe, i);
        LVA[i] = _mesh_model.L_atten(idxe, i);
        NVA[i] = _mesh_model.N_atten(idxe, i);
        IDXW[i] = this->LtG_T(idxe, i);
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
    vec_ke_t_base.push_back(make_smat(tpl_ke_1));
    vec_ke_t_base.push_back(make_smat(tpl_ke_2));
    vec_ke_t_atten.push_back(make_smat(tpl_ke_a1));
    vec_ke_t_atten.push_back(make_smat(tpl_ke_a2));
  }

  // -----------------------------------------------------------------------
  // spheroidal base matrices
  {
    vec_ke_s_base.reserve(3);
    vec_ke_s_atten.reserve(3);
    vec_in_s_base.reserve(2);
    auto totlen_S = this->LtG_S(2, _mesh.NE() - 1, NQ - 1) + 1;
    using T = Eigen::Triplet<double>;
    std::vector<T> tpl_in_0, tpl_in_1;
    std::vector<T> tpl_ke_0, tpl_ke_1, tpl_ke_2;
    std::vector<T> tpl_ke_a0, tpl_ke_a1, tpl_ke_a2;

    // inertia
    for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
      for (int i = 0; i < q.N(); ++i) {
        double xrad = _mesh.NodeRadius(idxe, i);
        double tmp = _mesh.EW(idxe) / 2.0 * q.W(i) *
                     _mesh_model.Density(idxe, i) * xrad * xrad;
        auto idx_uu = this->LtG_S(0, idxe, i);
        auto idx_vv = this->LtG_S(1, idxe, i);
        tpl_in_0.push_back(T(idx_uu, idx_uu, tmp));
        tpl_in_1.push_back(T(idx_vv, idx_vv, tmp));
      }
    }

    // diagonal stiffness
    for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
      double elem_width = _mesh.EW(idxe);
      for (int i = 0; i < q.N(); ++i) {
        double tmp0 = elem_width / 2.0 * q.W(i);
        auto Li = _mesh_model.L(idxe, i);
        auto Ai = _mesh_model.A(idxe, i);
        auto Ni = _mesh_model.N(idxe, i);
        auto gi = _mesh_model.Gravity(idxe, i);
        auto crad = _mesh.NodeRadius(idxe, i);
        auto crho = _mesh_model.Density(idxe, i);
        auto Li_a = _mesh_model.L_atten(idxe, i);
        auto Ai_a = _mesh_model.A_atten(idxe, i);
        auto Ni_a = _mesh_model.N_atten(idxe, i);
        auto idxtiu = this->LtG_S(0, idxe, i);
        auto idxtiv = this->LtG_S(1, idxe, i);
        auto idxtip = this->LtG_S(2, idxe, i);

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
    for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
      double elem_width = _mesh.EW(idxe);
      double e2 = 0.5 * elem_width;
      double d_val = 2.0 / elem_width;
      stdvvec mat_d(q.N(), stdvec(q.N(), 0.0));
      for (int i = 0; i < q.N(); ++i)
        for (int j = 0; j < q.N(); ++j)
          mat_d[i][j] = d_val * vec_lag_deriv[j][i];

      std::vector<double> LV(q.N()), NV(q.N()), AV(q.N()), FV(q.N()), RR(q.N());
      std::vector<double> LVA(q.N()), NVA(q.N()), AVA(q.N()), FVA(q.N());
      for (int i = 0; i < q.N(); ++i) {
        LV[i] = _mesh_model.L(idxe, i);
        NV[i] = _mesh_model.N(idxe, i);
        AV[i] = _mesh_model.A(idxe, i);
        FV[i] = _mesh_model.F(idxe, i);
        RR[i] = _mesh.NodeRadius(idxe, i);
        LVA[i] = _mesh_model.L_atten(idxe, i);
        NVA[i] = _mesh_model.N_atten(idxe, i);
        AVA[i] = _mesh_model.A_atten(idxe, i);
        FVA[i] = _mesh_model.F_atten(idxe, i);
      }

      for (int i = 0; i < q.N(); ++i) {
        auto ri = _mesh.NodeRadius(idxe, i);
        auto rhoi = _mesh_model.Density(idxe, i);
        for (int j = 0; j < q.N(); ++j) {
          auto rj = _mesh.NodeRadius(idxe, j);
          auto rhoj = _mesh_model.Density(idxe, j);
          auto tmp_uu_0 = elem_width * q.W(j) * FV[j] * rj * mat_d[i][j];
          auto tmp_vv_1 = -e2 * q.W(j) * LV[j] * rj * mat_d[i][j];
          auto tmp_uvd_1 = e2 * q.W(i) * LV[i] * ri * mat_d[j][i];
          auto tmp_udv_1 = -e2 * q.W(j) * FV[j] * rj * mat_d[i][j];
          auto tmp_pdu_0 = q.W(j) * rhoj * rj * rj * vec_lag_deriv[j][i];
          auto tmp_uu_0a = elem_width * q.W(j) * FVA[j] * rj * mat_d[i][j];
          auto tmp_vv_1a = -e2 * q.W(j) * LVA[j] * rj * mat_d[i][j];
          auto tmp_uvd_1a = e2 * q.W(i) * LVA[i] * ri * mat_d[j][i];
          auto tmp_udv_1a = -e2 * q.W(j) * FVA[j] * rj * mat_d[i][j];
          auto idx_u_i = this->LtG_S(0, idxe, i),
               idx_u_j = this->LtG_S(0, idxe, j);
          auto idx_v_i = this->LtG_S(1, idxe, i),
               idx_v_j = this->LtG_S(1, idxe, j);
          auto idx_p_i = this->LtG_S(2, idxe, i);
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
            auto rk = _mesh.NodeRadius(idxe, k);
            auto ddrr = rk * rk * mat_d[i][k] * mat_d[j][k];
            tmp_uu_0 += q.W(k) * _mesh_model.C(idxe, k) * ddrr;
            tmp_vv_1 += q.W(k) * _mesh_model.L(idxe, k) * ddrr;
            tmp_pp_0 += q.W(k) * ddrr;
            tmp_uu_0a += q.W(k) * _mesh_model.C_atten(idxe, k) * ddrr;
            tmp_vv_1a += q.W(k) * _mesh_model.L_atten(idxe, k) * ddrr;
          }
          auto idx_u_i = this->LtG_S(0, idxe, i),
               idx_u_j = this->LtG_S(0, idxe, j);
          auto idx_v_i = this->LtG_S(1, idxe, i),
               idx_v_j = this->LtG_S(1, idxe, j);
          auto idx_p_i = this->LtG_S(2, idxe, i),
               idx_p_j = this->LtG_S(2, idxe, j);
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
    vec_in_s_base.push_back(make_in(tpl_in_0));
    vec_in_s_base.push_back(make_in(tpl_in_1));
    vec_ke_s_base.push_back(make_smat(tpl_ke_0));
    vec_ke_s_base.push_back(make_smat(tpl_ke_1));
    vec_ke_s_base.push_back(make_smat(tpl_ke_2));
    vec_ke_s_atten.push_back(make_smat(tpl_ke_a0));
    vec_ke_s_atten.push_back(make_smat(tpl_ke_a1));
    vec_ke_s_atten.push_back(make_smat(tpl_ke_a2));
  }
};

}   // namespace Full1D

#endif   // SPECSEM_CONSTRUCTOR_H
