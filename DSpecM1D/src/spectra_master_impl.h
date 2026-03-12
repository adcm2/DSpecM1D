#ifndef SPECTRA_MASTER_IMPL_H
#define SPECTRA_MASTER_IMPL_H

#include "spectra_master.h"

namespace Full1D {

auto
specsem::LtG_S(int neig, int idx_e, int idx_n) const {
  assert((neig >= 0) && (neig <= 2) &&
         "Error: neig must be 0, 1 or 2 in LtG_S");
  assert((idx_e >= 0) && (idx_e < _mesh.NE()) &&
         "Error: idx_e out of range in LtG_S");
  assert((idx_n >= 0) && (idx_n < _mesh.NN()) &&
         "Error: idx_n out of range in LtG_S");
  std::size_t retval = 0;

  int offset_val = vec_offset[idx_e];
  if (idx_n == 0) {
    if ((std::find(fsb.begin(), fsb.end(), idx_e - 1) != fsb.end())) {
      if (neig == 1) {
        offset_val += 1;
      } else {
        offset_val -= 1;
      }
    }
  }
  retval = (3 * idx_e * (_mesh.NN() - 1) + idx_n * 3 + neig) + offset_val;
  return retval;
};

auto
specsem::LtG_R(int neig, int idx_e, int idx_n) const {
  assert((neig >= 0) && (neig < 2) && "Error: neig must be 0 or 1 in LtG_R");
  assert((idx_e >= 0) && (idx_e < _mesh.NE()) &&
         "Error: idx_e out of range in LtG_R");
  assert((idx_n >= 0) && (idx_n < _mesh.NN()) &&
         "Error: idx_n out of range in LtG_R");

  // first index in element
  std::size_t retval = 2 * idx_e * (_mesh.NN() - 1);

  // then node
  retval += idx_n * 2;

  // then function
  retval += neig;

  // return
  return retval;
};

auto
specsem::LtG_T(int idx_e, int idx_n) const {
  assert((idx_e >= _el) && (idx_e < _eu) && "idxe out of range in LtG_T");
  assert((idx_n >= 0) && (idx_n < _mesh.NN()) &&
         "Error: idx_n out of range in LtG_T");
  std::size_t retval = (idx_e - _el) * (_mesh.NN() - 1) + idx_n;
  return retval;
};

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

  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  // getting boundary and offset information for the local to global map
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

  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  // getting matrices of derivative values for Lagrange polynomials
  auto q = _mesh.GLL();
  auto pleg =
      Interpolation::LagrangePolynomial(q.Points().begin(), q.Points().end());
  {

    // storing h_i'(x_k)
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

  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  {

    // size of matrix for full problem
    totlen = this->LtG_S(2, _mesh.NE() - 1, NQ - 1) + 1;

    // triplet list for matrices
    using T = Eigen::Triplet<double>;

    // make the vector of matrices for the different ls
    {
      using T = Eigen::Triplet<double>;
      std::size_t matlen = totlen;

      ///////////////////////////////////////////////////////////////////////
      // toroidals
      // vector to store whether an element is fluid or not
      _vec_fluid = std::vector<int>(_mesh.NE(), 0);
      _vec_dof = std::vector<bool>(_mesh.NE(), false);
      bool _has_fluid = false;
      for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
        if (inp_model.IsFluid(_mesh.LayerNumber(idxe))) {
          _vec_fluid[idxe] = 1;
          _has_fluid = true;
        }

        // if there is a change from fluid to solid or vice versa, we set the
        // dof to true at the lower element
        if (idxe > 0) {
          if ((std::abs(_vec_fluid[idxe] - _vec_fluid[idxe - 1]) == 1)) {
            _vec_dof[idxe - 1] = true;
          }
        }
      }
      {
        bool found = false;
        bool prev_fluid = false;
        while (!found) {
          ++_el;
          auto laynum = _mesh.LayerNumber(_el);
          if (prev_fluid) {
            if (inp_model.IsSolid(laynum)) {
              found = true;
              break;
            }
          }
          if (inp_model.IsFluid(laynum)) {
            prev_fluid = true;
          }
        }
      }

      // check if there is an ocean layer
      _eu = _el;

      {
        bool found = false;
        while (!found) {
          auto laynum = _mesh.LayerNumber(_eu);
          if (inp_model.IsFluid(laynum) || (++_eu == _mesh.NE())) {
            found = true;
            break;
          }
        }
      }

      // radial inertia matrix
      {
        auto num_radial_dof = this->LtG_R(1, _mesh.NE() - 1, NQ - 1) + 1;
        std::vector<T> tpl_in, tpl_ke, tpl_ke_atten;
        {
          for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
            // int imin = (idxe == 0);
            double elem_width = _mesh.EW(idxe);
            int laynum = _mesh.LayerNumber(idxe);
            for (int i = 0; i < q.N(); ++i) {
              // pre-factor
              double xrad = _mesh.NodeRadius(idxe, i);
              double tmp = elem_width / 2.0 * q.W(i) *
                           _mesh_model.Density(idxe, i) * xrad * xrad;

              // indices
              auto idx_uu = this->LtG_R(0, idxe, i);

              // U'U term
              tpl_in.push_back(T(idx_uu, idx_uu, tmp));
            }
          }
          // Eigen::SparseMatrix<double> mat_tmp(totlen, totlen);
          mat_inertia_0.resize(num_radial_dof, num_radial_dof);
          mat_inertia_0.setFromTriplets(tpl_in.begin(), tpl_in.end());
          mat_inertia_0.makeCompressed();
        }
        // std::cout << "Post radial inertia matrix setting.\n";

        //  radial kinetic energy matrix
        {

          // diagonal terms
          {
            for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
              // int imin = (idxe == 0);
              double elem_width = _mesh.EW(idxe);

              int laynum = _mesh.LayerNumber(idxe);
              double d_val = 2.0 / elem_width;
              for (int i = 0; i < q.N(); ++i) {
                auto crad = _mesh.NodeRadius(idxe, i);
                auto crho = _mesh_model.Density(idxe, i);
                auto gi = _mesh_model.Gravity(idxe, i);
                // universal factors
                double tmp0 = elem_width / 2.0 * q.W(i);

                // material parameters
                auto Li = _mesh_model.L(idxe, i);
                auto Fi = _mesh_model.F(idxe, i);
                auto Ai = _mesh_model.A(idxe, i);
                auto Ni = _mesh_model.N(idxe, i);
                auto Ai_atten = _mesh_model.A_atten(idxe, i);
                auto Ni_atten = _mesh_model.N_atten(idxe, i);
                //////////////////////////////
                // U'U term
                double tmp_u =
                    tmp0 *
                    (4.0 * crho * (pi_db * bigg_db * crho * crad - gi) * crad +
                     4 * (Ai - Ni));
                double tmp_u_atten = 4 * tmp0 * (Ai_atten - Ni_atten);

                //////////////////////////////

                // indices
                std::size_t idxtiu = this->LtG_R(0, idxe, i);

                //////////////////////////////
                // pushback
                tpl_ke.push_back(T(idxtiu, idxtiu, tmp_u));   // U'U term
                tpl_ke_atten.push_back(
                    T(idxtiu, idxtiu, tmp_u_atten));   // U'U term
                //////////////////////////////
              }
            }

            // add in boundary term for gravity
            int idxpb = this->LtG_R(1, _mesh.NE() - 1, NQ - 1);
            double rpb = _mesh.NodeRadius(_mesh.NE() - 1, NQ - 1);
            tpl_ke.push_back(T(idxpb, idxpb, rpb / (4.0 * pi_db * bigg_db)));
          }
          // std::cout << "Post radial KE diagonal terms.\n";
          // coupling terms
          {
            for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
              // int imin = (idxe == 0);
              double elem_width = _mesh.EW(idxe);
              auto e2 = 0.5 * _mesh.EW(idxe);
              int laynum = _mesh.LayerNumber(idxe);
              double d_val = 2.0 / _mesh.EW(idxe);

              // calculate M and N
              stdvvec mat_d(q.N(), stdvec(q.N(), 0.0));
              // stdvvec mat_n = mat_m;

              for (int i = 0; i < q.N(); ++i) {
                for (int j = 0; j < q.N(); ++j) {
                  mat_d[i][j] = d_val * vec_lag_deriv[j][i];
                }
              }

              std::vector<double> LV(q.N(), 0.0), NV(q.N(), 0.0),
                  AV(q.N(), 0.0), FV(q.N(), 0.0), RR(q.N(), 0.0);
              std::vector<double> LVA(q.N(), 0.0), NVA(q.N(), 0.0),
                  AVA(q.N(), 0.0), FVA(q.N(), 0.0);
              for (int i = 0; i < q.N(); ++i) {
                LV[i] = _mesh_model.L(idxe, i);
                NV[i] = _mesh_model.N(idxe, i);
                AV[i] = _mesh_model.A(idxe, i);
                FV[i] = _mesh_model.F(idxe, i);
                RR[i] = _mesh.NodeRadius(idxe, i);

                // attenuation terms can be added here later
                LVA[i] = _mesh_model.L_atten(idxe, i);
                NVA[i] = _mesh_model.N_atten(idxe, i);
                AVA[i] = _mesh_model.A_atten(idxe, i);
                FVA[i] = _mesh_model.F_atten(idxe, i);
              }

              for (int i = 0; i < q.N(); ++i) {
                auto ri = _mesh.NodeRadius(idxe, i);
                auto rhoi = _mesh_model.Density(idxe, i);
                auto Ai = _mesh_model.A(idxe, i);
                auto Ni = _mesh_model.N(idxe, i);
                auto Fi = _mesh_model.F(idxe, i);
                auto Li = _mesh_model.L(idxe, i);
                for (int j = 0; j < q.N(); ++j) {
                  auto rj = _mesh.NodeRadius(idxe, j);
                  auto rhoj = _mesh_model.Density(idxe, j);
                  // auto Aj = _mesh_model.A(idxe, j);
                  // auto Nj = _mesh_model.N(idxe, j);
                  // auto Fj = _mesh_model.F(idxe, j);
                  // auto Lj = _mesh_model.L(idxe, j);

                  // finding values to push back
                  auto tmp_uu =
                      elem_width * q.W(j) * FV[j] * RR[j] * mat_d[i][j];

                  auto tmp_pdu =
                      q.W(j) * rhoj * RR[j] * RR[j] * vec_lag_deriv[j][i];

                  // attenuation terms
                  auto tmp_uu_atten =
                      elem_width * q.W(j) * FVA[j] * RR[j] * mat_d[i][j];

                  // indices
                  auto idx_u_i = this->LtG_R(0, idxe, i);
                  auto idx_u_j = this->LtG_R(0, idxe, j);
                  auto idx_p_i = this->LtG_R(1, idxe, i);
                  auto idx_p_j = this->LtG_R(1, idxe, j);

                  // pushbacks
                  tpl_ke.push_back(T(idx_u_i, idx_u_j, tmp_uu));    // U'U
                  tpl_ke.push_back(T(idx_u_j, idx_u_i, tmp_uu));    // UU'
                  tpl_ke.push_back(T(idx_p_i, idx_u_j, tmp_pdu));   // P'U
                  tpl_ke.push_back(T(idx_u_j, idx_p_i, tmp_pdu));   // UP'

                  // attenuation pushbacks
                  tpl_ke_atten.push_back(
                      T(idx_u_i, idx_u_j, tmp_uu_atten));   // U'U
                  tpl_ke_atten.push_back(
                      T(idx_u_j, idx_u_i, tmp_uu_atten));   // UU'
                }
              }

              for (int i = 0; i < q.N(); ++i) {
                for (int j = 0; j < q.N(); ++j) {
                  auto tmp_uu = 0.0;
                  auto tmp_pp = 0.0;
                  auto tmp_uu_atten = 0.0;

                  // go over element
                  for (int k = 0; k < q.N(); ++k) {
                    auto rk = _mesh.NodeRadius(idxe, k);
                    auto tmp_ddrr = rk * rk * mat_d[i][k] * mat_d[j][k];
                    tmp_uu += q.W(k) * _mesh_model.C(idxe, k) * tmp_ddrr;

                    tmp_pp += q.W(k) * tmp_ddrr;

                    tmp_uu_atten +=
                        q.W(k) * _mesh_model.C_atten(idxe, k) * tmp_ddrr;
                  }

                  // multiply
                  tmp_uu *= e2;
                  tmp_pp *= e2 / (4.0 * pi_db * bigg_db);
                  tmp_uu_atten *= e2;

                  // indices
                  auto idx_u_i = this->LtG_R(0, idxe, i);
                  auto idx_u_j = this->LtG_R(0, idxe, j);
                  auto idx_p_i = this->LtG_R(1, idxe, i);
                  auto idx_p_j = this->LtG_R(1, idxe, j);

                  // pushbacks
                  tpl_ke.push_back(T(idx_u_i, idx_u_j, tmp_uu));   // U'U
                  tpl_ke.push_back(T(idx_p_i, idx_p_j, tmp_pp));   // P'P

                  // attenuation pushbacks
                  tpl_ke_atten.push_back(
                      T(idx_u_i, idx_u_j, tmp_uu_atten));   // U'U
                }
              }
            }
          }

          {
            mat_ke_0.resize(num_radial_dof, num_radial_dof);
            mat_ke_0.setFromTriplets(tpl_ke.begin(), tpl_ke.end());
            mat_ke_0.makeCompressed();

            // attenuation matrix
            mat_ke_0_atten.resize(num_radial_dof, num_radial_dof);
            mat_ke_0_atten.setFromTriplets(tpl_ke_atten.begin(),
                                           tpl_ke_atten.end());
            mat_ke_0_atten.makeCompressed();
          }
          // std::cout << "Post radial KE FINAL.\n";
          // std::cout << "Completed L = " << idxl << " matrices\n";
          // }
        }
      }
    }
  }
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  // toroidals base matrices
  {
    // vec_in_t_base.reserve(1);
    using T = Eigen::Triplet<double>;
    std::vector<T> tpl_in_0, tpl_ke_1, tpl_ke_2;
    std::vector<T> tpl_ke_atten_1, tpl_ke_atten_2;
    auto ntdof = this->LtG_T(_mesh.NE() - 1, NQ - 1) + 1;
    {
      for (int idxe = _el; idxe < _eu; ++idxe) {
        double elem_width = _mesh.EW(idxe);
        for (int i = 0; i < q.N(); ++i) {
          // pre-factor
          double xrad = _mesh.NodeRadius(idxe, i);
          double tmp = _mesh.EW(idxe) / 2.0;
          tmp *= q.W(i) * _mesh_model.Density(idxe, i);
          tmp *= xrad * xrad;

          // indices and pushback
          auto idx_ww = this->LtG_T(idxe, i);
          tpl_in_0.push_back(T(idx_ww, idx_ww, tmp));
        }
      }
      mat_in_t_base.resize(ntdof, ntdof);
      mat_in_t_base.setFromTriplets(tpl_in_0.begin(), tpl_in_0.end());
      mat_in_t_base.makeCompressed();
    }

    {

      for (int idxe = _el; idxe < _eu; ++idxe) {
        // int imin = (idxe == 0);
        stdvvec mat_d(q.N(), stdvec(q.N(), 0.0));
        // stdvvec mat_n = mat_m;
        double elem_width = _mesh.EW(idxe);
        double e2 = 0.5 * _mesh.EW(idxe);
        int laynum = _mesh.LayerNumber(idxe);
        double d_val = 2.0 / elem_width;

        for (int i = 0; i < q.N(); ++i) {
          for (int j = 0; j < q.N(); ++j) {
            mat_d[i][j] = d_val * vec_lag_deriv[j][i];
          }
        }
        std::vector<double> LV(q.N(), 0.0), NV(q.N(), 0.0), RR(q.N(), 0.0);
        std::vector<double> LVA(q.N(), 0.0), NVA(q.N(), 0.0);
        std::vector<int> IDXW(q.N(), 0);
        for (int i = 0; i < q.N(); ++i) {
          LV[i] = _mesh_model.L(idxe, i);
          NV[i] = _mesh_model.N(idxe, i);
          RR[i] = _mesh.NodeRadius(idxe, i);
          IDXW[i] = this->LtG_T(idxe, i);

          // attenuation terms can be added here later
          LVA[i] = _mesh_model.L_atten(idxe, i);
          NVA[i] = _mesh_model.N_atten(idxe, i);
        }

        for (int i = 0; i < q.N(); ++i) {
          double tmp0 = elem_width / 2.0 * q.W(i);
          auto tmp_ww_1 = tmp0 * (LV[i] - 2.0 * NV[i]);
          auto tmp_ww_2 = tmp0 * NV[i];

          tpl_ke_1.push_back(T(IDXW[i], IDXW[i], tmp_ww_1));   // W'W term
          tpl_ke_2.push_back(T(IDXW[i], IDXW[i], tmp_ww_2));   // W'W term

          // attenuation terms
          auto tmp_w_atten_1 = tmp0 * (LVA[i] - 2.0 * NVA[i]);
          auto tmp_w_atten_2 = tmp0 * NVA[i];
          tpl_ke_atten_1.push_back(T(IDXW[i], IDXW[i], tmp_w_atten_1));
          tpl_ke_atten_2.push_back(T(IDXW[i], IDXW[i], tmp_w_atten_2));
        }

        for (int i = 0; i < q.N(); ++i) {
          for (int j = 0; j < q.N(); ++j) {
            auto tmp_wdw = -e2 * q.W(j) * LV[j] * RR[j] * mat_d[i][j];
            tpl_ke_1.push_back(T(IDXW[i], IDXW[j], tmp_wdw));
            tpl_ke_1.push_back(T(IDXW[j], IDXW[i], tmp_wdw));

            // attenuation terms
            auto tmp_wdw_atten = -e2 * q.W(j) * LVA[j] * RR[j] * mat_d[i][j];
            tpl_ke_atten_1.push_back(T(IDXW[i], IDXW[j], tmp_wdw_atten));
            tpl_ke_atten_1.push_back(T(IDXW[j], IDXW[i], tmp_wdw_atten));
          }
        }

        for (int i = 0; i < q.N(); ++i) {
          for (int j = 0; j < q.N(); ++j) {
            auto tmp_ww = 0.0;
            auto tmp_ww_atten = 0.0;

            // go over element
            for (int k = 0; k < q.N(); ++k) {
              auto tmp_ddrr = RR[k] * RR[k] * mat_d[i][k] * mat_d[j][k];
              tmp_ww += q.W(k) * LV[k] * tmp_ddrr;
              tmp_ww_atten += q.W(k) * LVA[k] * tmp_ddrr;
            }
            tmp_ww *= e2;
            tmp_ww_atten *= e2;

            tpl_ke_1.push_back(T(IDXW[i], IDXW[j], tmp_ww));
            tpl_ke_atten_1.push_back(T(IDXW[i], IDXW[j], tmp_ww_atten));
          }
        }
      }

      // setting matrices
      Eigen::SparseMatrix<double> mat_tmp_0(ntdof, ntdof);
      mat_tmp_0.setFromTriplets(tpl_ke_1.begin(), tpl_ke_1.end());
      mat_tmp_0.makeCompressed();
      vec_ke_t_base.push_back(mat_tmp_0);

      Eigen::SparseMatrix<double> mat_tmp_1(ntdof, ntdof);
      mat_tmp_1.setFromTriplets(tpl_ke_2.begin(), tpl_ke_2.end());
      mat_tmp_1.makeCompressed();
      vec_ke_t_base.push_back(mat_tmp_1);

      // attenuation matrices
      Eigen::SparseMatrix<double> mat_tmp_atten_0(ntdof, ntdof);
      mat_tmp_atten_0.setFromTriplets(tpl_ke_atten_1.begin(),
                                      tpl_ke_atten_1.end());
      mat_tmp_atten_0.makeCompressed();
      vec_ke_t_atten.push_back(mat_tmp_atten_0);

      Eigen::SparseMatrix<double> mat_tmp_atten_1(ntdof, ntdof);
      mat_tmp_atten_1.setFromTriplets(tpl_ke_atten_2.begin(),
                                      tpl_ke_atten_2.end());
      mat_tmp_atten_1.makeCompressed();
      vec_ke_t_atten.push_back(mat_tmp_atten_1);
    }
  }

  // spheroidals base matrices
  {
    vec_ke_s_base.reserve(3);
    vec_ke_s_atten.reserve(3);
    vec_in_s_base.reserve(2);
    auto totlen_S = this->LtG_S(2, _mesh.NE() - 1, NQ - 1) + 1;
    using T = Eigen::Triplet<double>;
    std::vector<T> tpl_in_0, tpl_in_1, tpl_ke_0, tpl_ke_1, tpl_ke_2;
    std::vector<T> tpl_ke_atten_0, tpl_ke_atten_1, tpl_ke_atten_2;
    {
      for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
        for (int i = 0; i < q.N(); ++i) {
          // pre-factor
          double xrad = _mesh.NodeRadius(idxe, i);
          double tmp = _mesh.EW(idxe) / 2.0 * q.W(i) *
                       _mesh_model.Density(idxe, i) * xrad * xrad;

          // indices
          auto idx_uu = this->LtG_S(0, idxe, i);
          auto idx_vv = this->LtG_S(1, idxe, i);

          // push back
          tpl_in_0.push_back(T(idx_uu, idx_uu, tmp));   // U'U term
          tpl_in_1.push_back(T(idx_vv, idx_vv, tmp));   // V'V term
        }
      }

      // setting matrices
      Eigen::SparseMatrix<double> mat_tmp_0(totlen_S, totlen_S);
      mat_tmp_0.setFromTriplets(tpl_in_0.begin(), tpl_in_0.end());
      mat_tmp_0.makeCompressed();
      vec_in_s_base.push_back(mat_tmp_0);

      Eigen::SparseMatrix<double> mat_tmp_1(totlen_S, totlen_S);
      mat_tmp_1.setFromTriplets(tpl_in_1.begin(), tpl_in_1.end());
      mat_tmp_1.makeCompressed();
      vec_in_s_base.push_back(mat_tmp_1);
    }

    //  kinetic energy matrix

    {

      // diagonal terms
      {
        for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
          // int imin = (idxe == 0);
          double elem_width = _mesh.EW(idxe);

          for (int i = 0; i < q.N(); ++i) {
            // universal factors
            double tmp0 = elem_width / 2.0 * q.W(i);

            // material parameters
            auto Li = _mesh_model.L(idxe, i);
            auto Ai = _mesh_model.A(idxe, i);
            auto Ni = _mesh_model.N(idxe, i);
            auto Fi = _mesh_model.F(idxe, i);

            auto gi = _mesh_model.Gravity(idxe, i);
            auto crad = _mesh.NodeRadius(idxe, i);
            auto crho = _mesh_model.Density(idxe, i);

            // attenuation
            auto Li_atten = _mesh_model.L_atten(idxe, i);
            auto Ai_atten = _mesh_model.A_atten(idxe, i);
            auto Ni_atten = _mesh_model.N_atten(idxe, i);
            auto Fi_atten = _mesh_model.F_atten(idxe, i);

            //////////////////////////////
            // U'U term
            double tmp_u_0 =
                4.0 * tmp0 *
                (crho * (pi_db * bigg_db * crho * crad - gi) * crad + Ai - Ni);
            double tmp_u_1 = tmp0 * Li;
            double tmp_v_1 = tmp0 * (Li - 2 * Ni);
            double tmp_v_2 = tmp0 * Ai;
            double tmp_pp_1 = 1.0 / (4.0 * pi_db * bigg_db) * tmp0;
            double tmp_pv_1 = tmp0 * crho * crad;
            double tmp_uv_1 = tmp0 * (crho * gi * crad - Li - 2 * (Ai - Ni));

            // attenuation terms
            double tmp_u_0_atten = 4.0 * tmp0 * (Ai_atten - Ni_atten);
            double tmp_u_1_atten = tmp0 * Li_atten;
            double tmp_v_1_atten = tmp0 * (Li_atten - 2 * Ni_atten);
            double tmp_v_2_atten = tmp0 * Ai_atten;
            double tmp_uv_1_atten =
                -tmp0 * (Li_atten + 2 * (Ai_atten - Ni_atten));

            //////////////////////////////

            // indices
            std::size_t idxtiu = this->LtG_S(0, idxe, i);
            std::size_t idxtiv = this->LtG_S(1, idxe, i);
            std::size_t idxtip = this->LtG_S(2, idxe, i);

            //////////////////////////////
            // pushback
            tpl_ke_0.push_back(T(idxtiu, idxtiu, tmp_u_0));    // U'U term
            tpl_ke_1.push_back(T(idxtiu, idxtiu, tmp_u_1));    // U'U term
            tpl_ke_1.push_back(T(idxtiv, idxtiv, tmp_v_1));    // V'V term
            tpl_ke_2.push_back(T(idxtiv, idxtiv, tmp_v_2));    // V'V term
            tpl_ke_1.push_back(T(idxtip, idxtip, tmp_pp_1));   // P'P term
            tpl_ke_1.push_back(T(idxtip, idxtiv, tmp_pv_1));   // P'V term
            tpl_ke_1.push_back(T(idxtiv, idxtip, tmp_pv_1));   // V'P term
            tpl_ke_1.push_back(T(idxtiu, idxtiv, tmp_uv_1));   // U'V term
            tpl_ke_1.push_back(T(idxtiv, idxtiu, tmp_uv_1));   // V'U term
            //////////////////////////////

            // attenuation pushbacks
            tpl_ke_atten_0.push_back(T(idxtiu, idxtiu, tmp_u_0_atten));
            tpl_ke_atten_1.push_back(T(idxtiu, idxtiu, tmp_u_1_atten));
            tpl_ke_atten_1.push_back(T(idxtiv, idxtiv, tmp_v_1_atten));
            tpl_ke_atten_2.push_back(T(idxtiv, idxtiv, tmp_v_2_atten));
            tpl_ke_atten_1.push_back(T(idxtiu, idxtiv, tmp_uv_1_atten));
            tpl_ke_atten_1.push_back(T(idxtiv, idxtiu, tmp_uv_1_atten));
          }
        }
      }

      // coupling terms
      {
        for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
          // int imin = (idxe == 0);
          double elem_width = _mesh.EW(idxe);
          auto e2 = 0.5 * _mesh.EW(idxe);
          int laynum = _mesh.LayerNumber(idxe);
          double d_val = 2.0 / _mesh.EW(idxe);

          // calculate M and N
          stdvvec mat_d(q.N(), stdvec(q.N(), 0.0));

          for (int i = 0; i < q.N(); ++i) {
            for (int j = 0; j < q.N(); ++j) {
              mat_d[i][j] = d_val * vec_lag_deriv[j][i];
            }
          }

          std::vector<double> LV(q.N(), 0.0), NV(q.N(), 0.0), AV(q.N(), 0.0),
              FV(q.N(), 0.0), RR(q.N(), 0.0);
          std::vector<double> LVA(q.N(), 0.0), NVA(q.N(), 0.0), AVA(q.N(), 0.0),
              FVA(q.N(), 0.0);
          for (int i = 0; i < q.N(); ++i) {
            LV[i] = _mesh_model.L(idxe, i);
            NV[i] = _mesh_model.N(idxe, i);
            AV[i] = _mesh_model.A(idxe, i);
            FV[i] = _mesh_model.F(idxe, i);
            RR[i] = _mesh.NodeRadius(idxe, i);

            // attenuation terms can be added here later
            LVA[i] = _mesh_model.L_atten(idxe, i);
            NVA[i] = _mesh_model.N_atten(idxe, i);
            AVA[i] = _mesh_model.A_atten(idxe, i);
            FVA[i] = _mesh_model.F_atten(idxe, i);
          }

          for (int i = 0; i < q.N(); ++i) {
            auto ri = _mesh.NodeRadius(idxe, i);
            auto rhoi = _mesh_model.Density(idxe, i);
            // auto Ai = _mesh_model.A(idxe, i);
            // auto Ni = _mesh_model.N(idxe, i);
            // auto Fi = _mesh_model.F(idxe, i);
            // auto Li = _mesh_model.L(idxe, i);
            for (int j = 0; j < q.N(); ++j) {
              auto rj = _mesh.NodeRadius(idxe, j);
              auto rhoj = _mesh_model.Density(idxe, j);
              // auto Aj = _mesh_model.A(idxe, j);
              // auto Nj = _mesh_model.N(idxe, j);
              // auto Fj = _mesh_model.F(idxe, j);
              // auto Lj = _mesh_model.L(idxe, j);

              // finding values to push back
              auto tmp_uu_0 = elem_width * q.W(j) * FV[j] * rj * mat_d[i][j];
              auto tmp_vv_1 = -e2 * q.W(j) * LV[j] * rj * mat_d[i][j];
              auto tmp_uvd_1 = e2 * q.W(i) * LV[i] * ri * mat_d[j][i];
              auto tmp_udv_1 = -e2 * q.W(j) * FV[j] * rj * mat_d[i][j];
              auto tmp_pdu_0 = q.W(j) * rhoj * rj * rj * vec_lag_deriv[j][i];

              // attenuation terms
              auto tmp_uu_0_atten =
                  elem_width * q.W(j) * FVA[j] * rj * mat_d[i][j];
              auto tmp_vv_1_atten = -e2 * q.W(j) * LVA[j] * rj * mat_d[i][j];
              auto tmp_uvd_1_atten = e2 * q.W(i) * LVA[i] * ri * mat_d[j][i];
              auto tmp_udv_1_atten = -e2 * q.W(j) * FVA[j] * rj * mat_d[i][j];

              // indices
              auto idx_u_i = this->LtG_S(0, idxe, i);
              auto idx_u_j = this->LtG_S(0, idxe, j);
              auto idx_v_i = this->LtG_S(1, idxe, i);
              auto idx_v_j = this->LtG_S(1, idxe, j);
              auto idx_p_i = this->LtG_S(2, idxe, i);
              auto idx_p_j = this->LtG_S(2, idxe, j);

              // push backs
              tpl_ke_0.push_back(T(idx_u_i, idx_u_j, tmp_uu_0));    // U'U
              tpl_ke_0.push_back(T(idx_u_j, idx_u_i, tmp_uu_0));    // U'U
              tpl_ke_1.push_back(T(idx_v_i, idx_v_j, tmp_vv_1));    // V'V
              tpl_ke_1.push_back(T(idx_v_j, idx_v_i, tmp_vv_1));    // V'V
              tpl_ke_1.push_back(T(idx_u_i, idx_v_j, tmp_uvd_1));   // U'V
              tpl_ke_1.push_back(T(idx_v_j, idx_u_i, tmp_uvd_1));   // VU'
              tpl_ke_1.push_back(T(idx_u_i, idx_v_j, tmp_udv_1));   // UV'
              tpl_ke_1.push_back(T(idx_v_j, idx_u_i, tmp_udv_1));   // V'U
              tpl_ke_0.push_back(T(idx_p_i, idx_u_j, tmp_pdu_0));   // P'U
              tpl_ke_0.push_back(T(idx_u_j, idx_p_i, tmp_pdu_0));   // UP'

              // attenuation push backs
              tpl_ke_atten_0.push_back(T(idx_u_i, idx_u_j, tmp_uu_0_atten));
              tpl_ke_atten_0.push_back(T(idx_u_j, idx_u_i, tmp_uu_0_atten));
              tpl_ke_atten_1.push_back(T(idx_v_i, idx_v_j, tmp_vv_1_atten));
              tpl_ke_atten_1.push_back(T(idx_v_j, idx_v_i, tmp_vv_1_atten));
              tpl_ke_atten_1.push_back(T(idx_u_i, idx_v_j, tmp_uvd_1_atten));
              tpl_ke_atten_1.push_back(T(idx_v_j, idx_u_i, tmp_uvd_1_atten));
              tpl_ke_atten_1.push_back(T(idx_u_i, idx_v_j, tmp_udv_1_atten));
              tpl_ke_atten_1.push_back(T(idx_v_j, idx_u_i, tmp_udv_1_atten));
            }
          }

          for (int i = 0; i < q.N(); ++i) {
            for (int j = 0; j < q.N(); ++j) {
              auto tmp_uu_0 = 0.0;
              auto tmp_vv_1 = 0.0;
              auto tmp_pp_0 = 0.0;

              // attenuation terms
              auto tmp_uu_0_atten = 0.0;
              auto tmp_vv_1_atten = 0.0;

              // go over element
              for (int k = 0; k < q.N(); ++k) {
                auto rk = _mesh.NodeRadius(idxe, k);
                auto tmp_ddrr = rk * rk * mat_d[i][k] * mat_d[j][k];
                tmp_uu_0 += q.W(k) * _mesh_model.C(idxe, k) * tmp_ddrr;
                tmp_vv_1 += q.W(k) * _mesh_model.L(idxe, k) * tmp_ddrr;
                tmp_pp_0 += q.W(k) * tmp_ddrr;

                // attenuation terms
                tmp_uu_0_atten +=
                    q.W(k) * _mesh_model.C_atten(idxe, k) * tmp_ddrr;
                tmp_vv_1_atten +=
                    q.W(k) * _mesh_model.L_atten(idxe, k) * tmp_ddrr;
              }

              // multiply
              tmp_uu_0 *= e2;
              tmp_vv_1 *= e2;
              tmp_pp_0 *= e2 / (4.0 * pi_db * bigg_db);

              // attenuation multiply
              tmp_uu_0_atten *= e2;
              tmp_vv_1_atten *= e2;

              // indices
              auto idx_u_i = this->LtG_S(0, idxe, i);
              auto idx_u_j = this->LtG_S(0, idxe, j);
              auto idx_v_i = this->LtG_S(1, idxe, i);
              auto idx_v_j = this->LtG_S(1, idxe, j);
              auto idx_p_i = this->LtG_S(2, idxe, i);
              auto idx_p_j = this->LtG_S(2, idxe, j);

              // pushbacks
              tpl_ke_0.push_back(T(idx_u_i, idx_u_j, tmp_uu_0));   // U'U
              tpl_ke_1.push_back(T(idx_v_i, idx_v_j, tmp_vv_1));   // V'V
              tpl_ke_0.push_back(T(idx_p_i, idx_p_j, tmp_pp_0));   // P'P

              // attenuation pushbacks
              tpl_ke_atten_0.push_back(T(idx_u_i, idx_u_j, tmp_uu_0_atten));
              tpl_ke_atten_1.push_back(T(idx_v_i, idx_v_j, tmp_vv_1_atten));
            }
          }
        }
      }
      // setting matrices
      Eigen::SparseMatrix<double> mat_tmp_0(totlen_S, totlen_S);
      mat_tmp_0.setFromTriplets(tpl_ke_0.begin(), tpl_ke_0.end());
      mat_tmp_0.makeCompressed();
      vec_ke_s_base.push_back(mat_tmp_0);

      Eigen::SparseMatrix<double> mat_tmp_1(totlen_S, totlen_S);
      mat_tmp_1.setFromTriplets(tpl_ke_1.begin(), tpl_ke_1.end());
      mat_tmp_1.makeCompressed();
      vec_ke_s_base.push_back(mat_tmp_1);

      Eigen::SparseMatrix<double> mat_tmp_2(totlen_S, totlen_S);
      mat_tmp_2.setFromTriplets(tpl_ke_2.begin(), tpl_ke_2.end());
      mat_tmp_2.makeCompressed();
      vec_ke_s_base.push_back(mat_tmp_2);

      // attenuation matrices
      Eigen::SparseMatrix<double> mat_tmp_atten_0(totlen_S, totlen_S);
      mat_tmp_atten_0.setFromTriplets(tpl_ke_atten_0.begin(),
                                      tpl_ke_atten_0.end());
      mat_tmp_atten_0.makeCompressed();
      vec_ke_s_atten.push_back(mat_tmp_atten_0);

      Eigen::SparseMatrix<double> mat_tmp_atten_1(totlen_S, totlen_S);
      mat_tmp_atten_1.setFromTriplets(tpl_ke_atten_1.begin(),
                                      tpl_ke_atten_1.end());
      mat_tmp_atten_1.makeCompressed();
      vec_ke_s_atten.push_back(mat_tmp_atten_1);

      Eigen::SparseMatrix<double> mat_tmp_atten_2(totlen_S, totlen_S);
      mat_tmp_atten_2.setFromTriplets(tpl_ke_atten_2.begin(),
                                      tpl_ke_atten_2.end());
      mat_tmp_atten_2.makeCompressed();
      vec_ke_s_atten.push_back(mat_tmp_atten_2);
    }
  }

  // std::cout << "Finished building base matrices.\n";
};

Eigen::SparseMatrix<double>
specsem::H_S(int idxl) const {
  auto k2 = idxl * (idxl + 1);
  Eigen::SparseMatrix<double> tmp = vec_ke_s_base[0];
  tmp += vec_ke_s_base[1] * k2;
  tmp += vec_ke_s_base[2] * k2 * k2;

  // add in P'P boundary
  auto idx_pp = this->LtG_S(2, _mesh.NE() - 1, _NQ - 1);
  double rpb = _mesh.NodeRadius(_mesh.NE() - 1, _NQ - 1);
  tmp.coeffRef(idx_pp, idx_pp) += rpb * (idxl + 1) / (4.0 * pi_db * bigg_db);

  // return
  return tmp;
};

Eigen::SparseMatrix<double>
specsem::H_SA(int idxl) const {
  auto k2 = idxl * (idxl + 1);
  Eigen::SparseMatrix<double> tmp = vec_ke_s_atten[0];
  tmp += vec_ke_s_atten[1] * k2;
  tmp += vec_ke_s_atten[2] * k2 * k2;

  // return
  return tmp;
};

Eigen::SparseMatrix<double>
specsem::P_S(int idxl) const {
  auto k2 = idxl * (idxl + 1);
  Eigen::SparseMatrix<double> tmp = vec_in_s_base[0];
  tmp += vec_in_s_base[1] * k2;
  return tmp;
};

Eigen::SparseMatrix<double>
specsem::H_TK(int idxl) const {
  auto k2 = idxl * (idxl + 1);
  Eigen::SparseMatrix<double> tmp = vec_ke_t_base[0] * k2;
  tmp += vec_ke_t_base[1] * k2 * k2;
  return tmp;
};

Eigen::SparseMatrix<double>
specsem::H_TA(int idxl) const {
  auto k2 = idxl * (idxl + 1);
  Eigen::SparseMatrix<double> tmp = vec_ke_t_atten[0] * k2;
  tmp += vec_ke_t_atten[1] * k2 * k2;
  return tmp;
};

Eigen::SparseMatrix<double>
specsem::P_TK(int idxl) const {
  auto k2 = idxl * (idxl + 1);
  Eigen::SparseMatrix<double> tmp = mat_in_t_base * k2;
  return tmp;
};

Eigen::MatrixXcd
specsem::CalculateForce(SourceInfo::EarthquakeCMT &cmt, int idxl) {

  //////////////////////////////
  int NQ = _mesh.NN();
  totlen = this->LtG_S(2, _mesh.NE() - 1, NQ - 1) + 1;
  Eigen::MatrixXcd vec_lforce = Eigen::MatrixXcd::Zero(totlen, 2 * idxl + 1);
  double kval =
      std::sqrt(static_cast<double>(idxl) * (static_cast<double>(idxl) + 1.0));
  double kd2 = kval / std::sqrt(2.0);

  // find element within which the source sits
  double rad_source = _mesh.PR() - 1000.0 * cmt.Depth() / _length_norm;

  // source location in spherical coordinates
  double theta_s = (90.0 - cmt.Latitude()) * EIGEN_PI / (180.0);
  double phi_s = cmt.Longitude() * EIGEN_PI / (180.0);

  //////////////////////////////
  // wigner matrix for evaluation
  int maxn = 2;
  if (maxn > idxl) {
    maxn = idxl;
  }
  auto wigdmat =
      GSHTrans::Wigner<double, GSHTrans::Ortho, GSHTrans::All, GSHTrans::All,
                       GSHTrans::Single, GSHTrans::ColumnMajor>(idxl, idxl,
                                                                maxn, theta_s);
  // ylmn lambda
  auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
    auto dl = wigdmat[N];
    auto tmp = dl[l, m];
    auto ylm = tmp * std::exp(Complex(0.0, m * phi));
    return ylm;
  };

  //////////////////////////////
  // parameters
  double invsq2 = 1.0 / std::sqrt(2.0);
  Complex isq2 = Complex(0.0, 1.0 / std::sqrt(2.0));
  double omegal2 = std::sqrt((idxl + 2) * (idxl - 1) / 2.0);
  // double lprefac = std::exp(-2.0 * this->pi_db * (idxl + 1) / (1 + 0.5));
  double lprefac = 1.0;

  // loop through the elements
  for (int idx = 0; idx < _mesh.NE(); ++idx) {
    if ((_mesh.ELR(idx) <= rad_source) && (_mesh.EUR(idx) > rad_source)) {
      // std::cout << "\nSource located in element " << idx << "\n\n";
      //////////////////////////////
      stdvec vec_nodes(NQ, 0.0);
      for (int idxn = 0; idxn < NQ; ++idxn) {
        vec_nodes[idxn] = _mesh.NodeRadius(idx, idxn);
      }

      auto pleg =
          Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());
      //////////////////////////////

      for (int idxq = 0; idxq < NQ; ++idxq) {
        auto w_val = pleg(idxq, rad_source) / rad_source;
        auto w_deriv = pleg.Derivative(idxq, rad_source);

        //////////////////////////////
        // indices
        auto idx_u = this->LtG_S(0, idx, idxq);
        auto idx_v = this->LtG_S(1, idx, idxq);

        // loop over m
        for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
          // spherical harmonic values
          Complex y0c = std::conj(ylmn(idxl, idxm, 0, phi_s));
          Complex ymc = std::conj(ylmn(idxl, idxm, -1, phi_s));
          Complex ypc = std::conj(ylmn(idxl, idxm, 1, phi_s));
          Complex ymmc = 0.0, yppc = 0.0;
          if (idxl > 1) {
            ymmc = std::conj(ylmn(idxl, idxm, -2, phi_s));
            yppc = std::conj(ylmn(idxl, idxm, 2, phi_s));
          }

          //////////////////////////////
          // common values
          Complex tmp_0pm = cmt.MC0p() * ypc + cmt.MC0m() * ymc;
          Complex tmp_pm = 2.0 * cmt.MCmp() * y0c;
          Complex tmp_ppmm = cmt.MCpp() * yppc + cmt.MCmm() * ymmc;

          // u component
          Complex tmp_u = w_val * (kd2 * tmp_0pm - tmp_pm);
          tmp_u += cmt.MC00() * w_deriv * y0c;

          // v component
          Complex tmp_v = w_deriv * tmp_0pm;
          tmp_v += w_val * (kd2 * tmp_pm - tmp_0pm + omegal2 * tmp_ppmm);
          tmp_v *= kd2;

          //////////////////////////////
          // put in force
          vec_lforce(idx_u, idxm + idxl) = tmp_u * lprefac;
          vec_lforce(idx_v, idxm + idxl) = tmp_v * lprefac;

          //////////////////////////////
        };
      };
    };
  };

  vec_lforce *= (1.0 / _moment_norm);
  return vec_lforce;
};

Eigen::MatrixXcd
specsem::CalculateForce_All(SourceInfo::EarthquakeCMT &cmt, int idxl) {

  //////////////////////////////
  int NQ = _mesh.NN();
  totlen = this->LtG_S(2, _mesh.NE() - 1, NQ - 1) + 1;
  Eigen::MatrixXcd vec_lforce = Eigen::MatrixXcd::Zero(totlen, 4);
  double kval =
      std::sqrt(static_cast<double>(idxl) * (static_cast<double>(idxl) + 1.0));
  double kd2 = kval / std::sqrt(2.0);

  // find element within which the source sits
  double rad_source = _mesh.PR() - 1000.0 * cmt.Depth() / _length_norm;

  // loop through the elements
  for (int idx = 0; idx < _mesh.NE(); ++idx) {
    if ((_mesh.ELR(idx) <= rad_source) && (_mesh.EUR(idx) > rad_source)) {
      // std::cout << "\nSource located in element " << idx << "\n\n";
      //////////////////////////////
      stdvec vec_nodes(NQ, 0.0);
      for (int idxn = 0; idxn < NQ; ++idxn) {
        vec_nodes[idxn] = _mesh.NodeRadius(idx, idxn);
      }

      auto pleg =
          Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());
      //////////////////////////////

      for (int idxq = 0; idxq < NQ; ++idxq) {
        auto w_val = pleg(idxq, rad_source) / rad_source;
        auto w_deriv = pleg.Derivative(idxq, rad_source);

        //////////////////////////////
        // indices
        auto idx_u = this->LtG_S(0, idx, idxq);
        auto idx_v = this->LtG_S(1, idx, idxq);

        // put in force
        vec_lforce(idx_u, 0) = w_deriv;
        vec_lforce(idx_u, 1) = w_val;
        vec_lforce(idx_v, 2) = kd2 * w_deriv;
        vec_lforce(idx_v, 3) = kd2 * w_val;
      };
    };
  };

  return vec_lforce;
};

auto
specsem::Source_Element(SourceInfo::EarthquakeCMT &cmt)
    const {   // find element within which the source sits
  double rad_source = _mesh.PR() - 1000.0 * cmt.Depth() / _length_norm;

  // loop through the elements
  int idxout = _mesh.NE() - 1;
  for (int idx = 0; idx < _mesh.NE(); ++idx) {
    if ((_mesh.ELR(idx) <= rad_source) && (_mesh.EUR(idx) > rad_source)) {
      idxout = idx;
      break;
    };
  };
  return idxout;
};

Eigen::MatrixXcd
specsem::CalculateForce_All_T(SourceInfo::EarthquakeCMT &cmt, int idxl) {

  //////////////////////////////
  int NQ = _mesh.NN();
  totlen = this->LtG_T(_mesh.NE() - 1, NQ - 1) + 1;
  Eigen::MatrixXcd vec_lforce = Eigen::MatrixXcd::Zero(totlen, 2);
  double kval =
      std::sqrt(static_cast<double>(idxl) * (static_cast<double>(idxl) + 1.0));
  double kd2 = kval / std::sqrt(2.0);

  // find element within which the source sits
  double rad_source = _mesh.PR() - 1000.0 * cmt.Depth() / _length_norm;

  // loop through the elements
  for (int idx = 0; idx < _mesh.NE(); ++idx) {
    if ((_mesh.ELR(idx) <= rad_source) && (_mesh.EUR(idx) > rad_source)) {
      // std::cout << "\nSource located in element " << idx << "\n\n";
      //////////////////////////////
      stdvec vec_nodes(NQ, 0.0);
      for (int idxn = 0; idxn < NQ; ++idxn) {
        vec_nodes[idxn] = _mesh.NodeRadius(idx, idxn);
      }

      auto pleg =
          Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());
      //////////////////////////////

      for (int idxq = 0; idxq < NQ; ++idxq) {
        auto w_val = pleg(idxq, rad_source) / rad_source;
        auto w_deriv = pleg.Derivative(idxq, rad_source);

        //////////////////////////////
        // indices
        auto idx_v = this->LtG_T(idx, idxq);

        // put in force
        {
          vec_lforce(idx_v, 0) = kval * w_deriv;
          vec_lforce(idx_v, 1) = kval * w_val;
        };
      };
    };
  };

  return vec_lforce;
};

Eigen::MatrixXcd
specsem::CalculateForce_Coefficients(SourceInfo::EarthquakeCMT &cmt, int idxl) {

  //////////////////////////////
  int NQ = _mesh.NN();
  totlen = this->LtG_S(2, _mesh.NE() - 1, NQ - 1) + 1;
  Eigen::MatrixXcd vec_lforce = Eigen::MatrixXcd::Zero(2 * idxl + 1, 4);
  double kval =
      std::sqrt(static_cast<double>(idxl) * (static_cast<double>(idxl) + 1.0));
  double kd2 = kval / std::sqrt(2.0);

  // source location in spherical coordinates
  double theta_s = (90.0 - cmt.Latitude()) * EIGEN_PI / (180.0);
  double phi_s = cmt.Longitude() * EIGEN_PI / (180.0);

  //////////////////////////////
  // wigner matrix for evaluation
  int maxn = 2;
  if (maxn > idxl) {
    maxn = idxl;
  }
  auto wigdmat =
      GSHTrans::Wigner<double, GSHTrans::Ortho, GSHTrans::All, GSHTrans::All,
                       GSHTrans::Single, GSHTrans::ColumnMajor>(idxl, idxl,
                                                                maxn, theta_s);
  // ylmn lambda
  auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
    auto dl = wigdmat[N];
    auto tmp = dl[l, m];
    auto ylm = tmp * std::exp(Complex(0.0, m * phi));
    return ylm;
  };

  //////////////////////////////
  // parameters
  // double invsq2 = 1.0 / std::sqrt(2.0);
  // Complex isq2 = Complex(0.0, 1.0 / std::sqrt(2.0));
  double omegal2 = std::sqrt((idxl + 2) * (idxl - 1) / 2.0);
  // double lprefac = std::exp(-2.0 * this->pi_db * (idxl + 1) / (1 + 0.5));
  // double lprefac = 1.0;

  // loop through the elements

  // loop over m
  for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
    // spherical harmonic values
    Complex y0c = std::conj(ylmn(idxl, idxm, 0, phi_s));
    Complex ymc = std::conj(ylmn(idxl, idxm, -1, phi_s));
    Complex ypc = std::conj(ylmn(idxl, idxm, 1, phi_s));
    Complex ymmc = 0.0, yppc = 0.0;
    if (idxl > 1) {
      ymmc = std::conj(ylmn(idxl, idxm, -2, phi_s));
      yppc = std::conj(ylmn(idxl, idxm, 2, phi_s));
    }

    //////////////////////////////
    // common values
    Complex tmp_0pm = cmt.MC0p() * ypc + cmt.MC0m() * ymc;
    Complex tmp_pm = 2.0 * cmt.MCmp() * y0c;
    Complex tmp_ppmm = cmt.MCpp() * yppc + cmt.MCmm() * ymmc;

    //////////////////////////////
    // put in force
    vec_lforce(idxm + idxl, 0) = cmt.MC00() * y0c;
    vec_lforce(idxm + idxl, 1) = kd2 * tmp_0pm - tmp_pm;
    vec_lforce(idxm + idxl, 2) = tmp_0pm;
    vec_lforce(idxm + idxl, 3) = kd2 * tmp_pm - tmp_0pm + omegal2 * tmp_ppmm;

    //////////////////////////////
  };

  vec_lforce *= (1.0 / _moment_norm);
  return vec_lforce;
};

Eigen::MatrixXcd
specsem::CalculateForce_RED_Coefficients(SourceInfo::EarthquakeCMT &cmt,
                                         int idxl, double az) {

  //////////////////////////////
  // force vector
  int maxn = 2;
  (maxn > idxl) ? maxn = idxl : maxn = maxn;
  Eigen::MatrixXcd vec_force = Eigen::MatrixXcd::Zero(2 * maxn + 1, 4);

  //////////////////////////////
  // parameters
  double lv = static_cast<double>(idxl);
  auto mfact = std::sqrt((2.0 * lv + 1.0) / (4.0 * EIGEN_PI));
  auto omegal2 = std::sqrt((lv + 2) * (lv - 1) / 2.0);
  auto kval = std::sqrt(lv * (lv + 1.0));
  auto kd2 = kval / std::sqrt(2.0);
  auto isq2 = Complex(0.0, 1.0 / std::sqrt(2.0));

  // exponential factors for azimuthal dependence
  Complex expm2 = std::exp(Complex(0.0, -2.0 * az));
  Complex expm1 = std::exp(Complex(0.0, -1.0 * az));
  Complex expp1 = std::exp(Complex(0.0, 1.0 * az));
  Complex expp2 = std::exp(Complex(0.0, 2.0 * az));
  expm2 = 1.0;
  expm1 = 1.0;
  expp1 = 1.0;
  expp2 = 1.0;

  // compute the combinations of CMT coefficients and spherical harmonics
  Complex tmp_mm = cmt.MCmm() * expm2;
  Complex tmp_0m = cmt.MC0m() * expm1;
  Complex tmp_0p = cmt.MC0p() * expp1;
  Complex tmp_pp = cmt.MCpp() * expp2;

  if (idxl == 1) {
    // partial_r U terms for l=1
    vec_force(1, 0) = cmt.MC00();

    // r^{-1} U terms for l=1
    vec_force(0, 1) = kd2 * tmp_0m;
    vec_force(1, 1) = -2.0 * cmt.MCmp();
    vec_force(2, 1) = kd2 * tmp_0p;

    // partial_r V terms for l=1
    vec_force(0, 2) = tmp_0m;
    vec_force(2, 2) = tmp_0p;

    // r^{-1} V terms for l=1
    vec_force(0, 3) = -tmp_0m;
    vec_force(1, 3) = kd2 * 2.0 * cmt.MCmp();
    vec_force(2, 3) = -tmp_0p;

  } else if (idxl > 1) {

    // partial_r U terms
    vec_force(2, 0) = cmt.MC00();

    // r^{-1} U terms
    vec_force(1, 1) = kd2 * tmp_0m;
    vec_force(2, 1) = -2.0 * cmt.MCmp();
    vec_force(3, 1) = kd2 * tmp_0p;

    // partial_r V terms
    vec_force(1, 2) = tmp_0m;
    vec_force(3, 2) = tmp_0p;

    // r^{-1} V terms
    vec_force(0, 3) = omegal2 * tmp_mm;
    vec_force(1, 3) = -tmp_0m;
    vec_force(2, 3) = kd2 * 2.0 * cmt.MCmp();
    vec_force(3, 3) = -tmp_0p;
    vec_force(4, 3) = omegal2 * tmp_pp;
  }
  vec_force *= mfact;
  // vec_force *= Complex(0, 1.0);

  vec_force *= (1.0 / _moment_norm);
  return vec_force;
};

Eigen::MatrixXcd
specsem::CalculateForce_Coefficients_T(SourceInfo::EarthquakeCMT &cmt,
                                       int idxl) {

  //////////////////////////////
  int NQ = _mesh.NN();
  totlen = this->LtG_T(_mesh.NE() - 1, NQ - 1) + 1;
  Eigen::MatrixXcd vec_lforce = Eigen::MatrixXcd::Zero(2 * idxl + 1, 2);
  double kval =
      std::sqrt(static_cast<double>(idxl) * (static_cast<double>(idxl) + 1.0));
  double kd2 = kval / std::sqrt(2.0);

  // source location in spherical coordinates
  double theta_s = (90.0 - cmt.Latitude()) * EIGEN_PI / (180.0);
  double phi_s = cmt.Longitude() * EIGEN_PI / (180.0);

  //////////////////////////////
  // wigner matrix for evaluation
  int maxn = 2;
  if (maxn > idxl) {
    maxn = idxl;
  }
  auto wigdmat =
      GSHTrans::Wigner<double, GSHTrans::Ortho, GSHTrans::All, GSHTrans::All,
                       GSHTrans::Single, GSHTrans::ColumnMajor>(idxl, idxl,
                                                                maxn, theta_s);
  // ylmn lambda
  auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
    auto dl = wigdmat[N];
    auto tmp = dl[l, m];
    auto ylm = tmp * std::exp(Complex(0.0, m * phi));
    return ylm;
  };

  //////////////////////////////
  // parameters
  Complex isq2 = Complex(0.0, 1.0 / std::sqrt(2.0));
  double omegal2 = std::sqrt((idxl + 2) * (idxl - 1) / 2.0);

  // loop over m
  for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
    // spherical harmonic values
    Complex ymc = std::conj(ylmn(idxl, idxm, -1, phi_s));
    Complex ypc = std::conj(ylmn(idxl, idxm, 1, phi_s));
    Complex ymmc = 0.0, yppc = 0.0;
    if (idxl > 1) {
      ymmc = std::conj(ylmn(idxl, idxm, -2, phi_s));
      yppc = std::conj(ylmn(idxl, idxm, 2, phi_s));
    }

    auto tmp1 = cmt.MC0m() * ymc - cmt.MC0p() * ypc;
    auto tmp2 = cmt.MCmm() * ymmc - cmt.MCpp() * yppc;
    //////////////////////////////
    // put in force
    vec_lforce(idxm + idxl, 0) = isq2 * tmp1;
    vec_lforce(idxm + idxl, 1) = isq2 * (omegal2 * tmp2 - tmp1);

    //////////////////////////////
  };

  vec_lforce *= (1.0 / _moment_norm);
  return vec_lforce;
};

Eigen::MatrixXcd
specsem::CalculateForce_RED_Coefficients_T(SourceInfo::EarthquakeCMT &cmt,
                                           int idxl, double az) {

  //////////////////////////////
  // force vector
  int maxn = 2;
  (maxn > idxl) ? maxn = idxl : maxn = maxn;
  Eigen::MatrixXcd vec_force = Eigen::MatrixXcd::Zero(2 * maxn + 1, 2);

  //////////////////////////////
  // parameters
  double lv = static_cast<double>(idxl);
  auto mfact = std::sqrt((2.0 * lv + 1.0) / (4.0 * EIGEN_PI));
  // mfact = 1.0;
  auto omegal2 = std::sqrt((lv + 2) * (lv - 1) / 2.0);
  auto kval = std::sqrt(lv * (lv + 1.0));
  auto kd2 = kval / std::sqrt(2.0);
  auto isq2 = Complex(0.0, 1.0 / std::sqrt(2.0));

  // exponential factors for azimuthal dependence
  Complex expm2 = std::exp(Complex(0.0, -2.0 * az));
  Complex expm1 = std::exp(Complex(0.0, -1.0 * az));
  Complex expp1 = std::exp(Complex(0.0, 1.0 * az));
  Complex expp2 = std::exp(Complex(0.0, 2.0 * az));
  expm2 = 1.0;
  expm1 = 1.0;
  expp1 = 1.0;
  expp2 = 1.0;

  // compute the combinations of CMT coefficients and spherical harmonics
  Complex tmp_mm = cmt.MCmm() * expm2;
  Complex tmp_0m = cmt.MC0m() * expm1;
  Complex tmp_0p = cmt.MC0p() * expp1;
  Complex tmp_pp = cmt.MCpp() * expp2;

  // different if l = 1 or l > 1

  if (idxl == 1) {
    // partial_r W terms for l=1
    vec_force(0, 0) = tmp_0m;
    vec_force(2, 0) = -tmp_0p;

    // r^{-1} W terms for l=1
    vec_force(0, 1) = -tmp_0m;
    vec_force(2, 1) = tmp_0p;

  } else if (idxl > 1) {
    // partial_r W terms for l>1
    vec_force(1, 0) = tmp_0m;
    vec_force(3, 0) = -tmp_0p;

    // r^{-1} W terms for l>1
    vec_force(1, 1) = -tmp_0m;
    vec_force(3, 1) = tmp_0p;

    // r^{-1} W terms for l>1
    vec_force(0, 1) = omegal2 * tmp_mm;
    vec_force(4, 1) = -omegal2 * tmp_pp;
  }
  vec_force *= mfact;
  vec_force *= isq2;

  vec_force *= (1.0 / _moment_norm);
  return vec_force;
};

Eigen::MatrixXcd
specsem::CalculateForce_T(SourceInfo::EarthquakeCMT &cmt, int idxl) {

  //////////////////////////////
  int NQ = _mesh.NN();
  totlen = this->LtG_T(_eu - 1, NQ - 1) + 1;
  Eigen::MatrixXcd vec_lforce = Eigen::MatrixXcd::Zero(totlen, 2 * idxl + 1);
  double kval =
      std::sqrt(static_cast<double>(idxl) * (static_cast<double>(idxl) + 1.0));
  double kd2 = kval / std::sqrt(2.0);

  // find element within which the source sits
  double rad_source = _mesh.PR() - 1000.0 * cmt.Depth() / _length_norm;

  // source location in spherical coordinates
  double theta_s = (90.0 - cmt.Latitude()) * EIGEN_PI / (180.0);
  double phi_s = cmt.Longitude() * EIGEN_PI / (180.0);

  //////////////////////////////
  // wigner matrix for evaluation
  int maxn = 2;
  if (maxn > _lmax) {
    maxn = _lmax;
  }
  auto wigdmat =
      GSHTrans::Wigner<double, GSHTrans::Ortho, GSHTrans::All, GSHTrans::All,
                       GSHTrans::Single, GSHTrans::ColumnMajor>(_lmax, _lmax,
                                                                maxn, theta_s);
  // ylmn lambda
  auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
    auto dl = wigdmat[N];
    auto tmp = dl[l, m];
    auto ylm = tmp * std::exp(Complex(0.0, m * phi));
    return ylm;
  };

  //////////////////////////////
  // parameters
  double invsq2 = 1.0 / std::sqrt(2.0);
  Complex isq2 = Complex(0.0, 1.0 / std::sqrt(2.0));
  double omegal2 = std::sqrt((idxl + 2) * (idxl - 1) / 2.0);
  // double lprefac = std::exp(-2.0 * this->pi_db * (idxl + 1) / (1 + 0.5));
  double lprefac = 1.0;

  // loop through the elements
  //   to find the element that contains the source
  for (int idx = _el; idx < _eu; ++idx) {
    // std::cout << _mesh.ELR(idx) << " " << _mesh.EUR(idx) << " " <<
    // rad_source
    //           << "\n";
    if ((_mesh.ELR(idx) <= rad_source) && (_mesh.EUR(idx) > rad_source)) {
      // std::cout << idx << " " << rad_source << " " << _mesh.ELR(idx) << " "
      //           << _mesh.EUR(idx) << "\n";
      std::vector<double> vec_nodes(NQ, 0.0);
      for (int idxn = 0; idxn < NQ; ++idxn) {
        vec_nodes[idxn] = _mesh.NodeRadius(idx, idxn);
      }

      auto pleg =
          Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());

      for (int idxq = 0; idxq < NQ; ++idxq) {
        auto w_val = pleg(idxq, rad_source) / rad_source;
        auto w_prefactor = pleg.Derivative(idxq, rad_source) -
                           w_val;   // prefactor for first term
        // for (int idxl = 1; idxl < _lmax + 1; ++idxl) {

        // get the index of the spherical harmonic
        std::size_t ridx = this->LtG_T(idx, idxq);
        for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
          // spherical harmonic
          Complex ymc = std::conj(ylmn(idxl, idxm, -1, phi_s));
          Complex ypc = std::conj(ylmn(idxl, idxm, 1, phi_s));
          Complex ymmc = 0.0, yppc = 0.0;

          if (idxl > 1) {
            ymmc = std::conj(ylmn(idxl, idxm, -2, phi_s));
            yppc = std::conj(ylmn(idxl, idxm, 2, phi_s));
          }

          // std::cout << "Size: " << vec_force.rows() << ", ovidx: " << ovidx
          //           << "\n";
          Complex tmp = w_prefactor * (cmt.MC0m() * ymc - cmt.MC0p() * ypc);
          tmp += w_val * omegal2 * (cmt.MCmm() * ymmc - cmt.MCpp() * yppc);
          tmp *= isq2;

          vec_lforce(ridx, idxm + idxl) = tmp;
        };
        // };
      };
    };
  };
  vec_lforce *= (1.0 / _moment_norm);
  return vec_lforce;
};

Eigen::MatrixXcd
specsem::CalculateForce_R(SourceInfo::EarthquakeCMT &cmt) {

  //////////////////////////////
  int NQ = _mesh.NN();
  totlen = this->LtG_R(1, _mesh.NE() - 1, NQ - 1) + 1;
  Eigen::MatrixXcd vec_lforce = Eigen::MatrixXcd::Zero(totlen, 1);
  // double kval =
  //     std::sqrt(static_cast<double>(idxl) * ( 1.0));
  // double kd2 = kval / std::sqrt(2.0);

  // find element within which the source sits
  double rad_source = _mesh.PR() - 1000.0 * cmt.Depth() / _length_norm;

  // source location in spherical coordinates
  double theta_s = (90.0 - cmt.Latitude()) * EIGEN_PI / (180.0);
  double phi_s = cmt.Longitude() * EIGEN_PI / (180.0);

  //////////////////////////////
  // wigner matrix for evaluation
  int maxn = 2;
  if (maxn > _lmax) {
    maxn = _lmax;
  }
  auto wigdmat =
      GSHTrans::Wigner<double, GSHTrans::Ortho, GSHTrans::All, GSHTrans::All,
                       GSHTrans::Single, GSHTrans::ColumnMajor>(_lmax, _lmax,
                                                                maxn, theta_s);
  auto wigdmat2 =
      GSHTrans::Wigner<double, GSHTrans::Ortho, GSHTrans::All, GSHTrans::All,
                       GSHTrans::Single, GSHTrans::ColumnMajor>(0, 0, 0, 0.0);
  // ylmn lambda
  auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
    auto dl = wigdmat[N];
    auto tmp = dl[l, m];
    auto ylm = tmp * std::exp(Complex(0.0, m * phi));
    return ylm;
  };

  //////////////////////////////
  // parameters
  double invsq2 = 1.0 / std::sqrt(2.0);
  Complex isq2 = Complex(0.0, 1.0 / std::sqrt(2.0));
  // double omegal2 = std::sqrt((idxl + 2) * (idxl - 1) / 2.0);
  // double lprefac = std::exp(-2.0 * this->pi_db * (idxl + 1) / (1 + 0.5));
  double lprefac = 1.0;

  // loop through the elements
  for (int idx = 0; idx < _mesh.NE(); ++idx) {
    if ((_mesh.ELR(idx) <= rad_source) && (_mesh.EUR(idx) > rad_source)) {
      //////////////////////////////
      stdvec vec_nodes(NQ, 0.0);
      for (int idxn = 0; idxn < NQ; ++idxn) {
        vec_nodes[idxn] = _mesh.NodeRadius(idx, idxn);
      }

      auto pleg =
          Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());
      //////////////////////////////

      for (int idxq = 0; idxq < NQ; ++idxq) {
        auto w_val = pleg(idxq, rad_source) / rad_source;
        auto w_deriv = pleg.Derivative(idxq, rad_source);

        //////////////////////////////
        // indices
        auto idx_u = this->LtG_R(0, idx, idxq);

        // values to use
        // Complex y0c = std::conj(ylmn(0, 0, 0, phi_s));
        Complex y0c = wigdmat2[0][0, 0];
        Complex tmp_pm = 2.0 * cmt.MCmp() * y0c;

        // u component
        Complex tmp_u = -w_val * tmp_pm;
        tmp_u += cmt.MC00() * w_deriv * y0c;

        //////////////////////////////
        // put in force
        vec_lforce(idx_u, 0) = tmp_u * lprefac;

        //////////////////////////////
      };
    };
  };

  vec_lforce *= (1.0 / _moment_norm);
  return vec_lforce;
};

Eigen::MatrixXcd
specsem::CalculateForce_Red_R(SourceInfo::EarthquakeCMT &cmt) {

  //////////////////////////////
  int NQ = _mesh.NN();
  totlen = this->LtG_R(1, _mesh.NE() - 1, NQ - 1) + 1;
  Eigen::MatrixXcd vec_lforce = Eigen::MatrixXcd::Zero(totlen, 1);
  // find element within which the source sits
  double rad_source = _mesh.PR() - 1000.0 * cmt.Depth() / _length_norm;

  auto wigdmat =
      GSHTrans::Wigner<double, GSHTrans::Ortho, GSHTrans::All, GSHTrans::All,
                       GSHTrans::Single, GSHTrans::ColumnMajor>(0, 0, 0, 0.0);
  Complex y0c = wigdmat[0][0, 0];
  // std::cout << "y0c: " << y0c << "\n";
  //////////////////////////////
  // loop through the elements
  for (int idx = 0; idx < _mesh.NE(); ++idx) {
    if ((_mesh.ELR(idx) <= rad_source) && (_mesh.EUR(idx) > rad_source)) {
      //////////////////////////////
      stdvec vec_nodes(NQ, 0.0);
      for (int idxn = 0; idxn < NQ; ++idxn) {
        vec_nodes[idxn] = _mesh.NodeRadius(idx, idxn);
      }

      auto pleg =
          Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());
      //////////////////////////////

      for (int idxq = 0; idxq < NQ; ++idxq) {
        auto w_val = pleg(idxq, rad_source) / rad_source;
        auto w_deriv = pleg.Derivative(idxq, rad_source);

        //////////////////////////////
        // indices
        auto idx_u = this->LtG_R(0, idx, idxq);

        // values to use

        Complex tmp_pm = 2.0 * cmt.MCmp() * y0c;

        // u component
        Complex tmp_u = -w_val * tmp_pm;
        tmp_u += cmt.MC00() * w_deriv * y0c;

        //////////////////////////////
        // put in force
        // std::cout << "Non-zero at: " << idx_u << "\n";
        vec_lforce(idx_u, 0) = tmp_u;

        //////////////////////////////
      };
    };
  };

  vec_lforce *= (1.0 / _moment_norm);
  return vec_lforce;
};

Eigen::MatrixXcd
specsem::RV_BASE_Z(InputParameters &param, int idxl, int idxr) {
  double rad_r = _mesh.PR() - param.receiver_depth() * 1000.0 / _length_norm;

  // create the receiver vector
  std::size_t flen = this->LtG_S(2, _mesh.NE() - 1, _mesh.NN() - 1) + 1;
  Eigen::MatrixXcd vec_receiver = Eigen::MatrixXcd::Zero(1, flen);

  // int elem_idx = 0;
  bool evaluated = false;
  for (int idx = 0; idx < _mesh.NE(); ++idx) {
    if ((_mesh.ELR(idx) <= rad_r) && (_mesh.EUR(idx) >= rad_r) &&
        (!evaluated)) {
      evaluated = true;

      // get the Lagrange polynomial for this element
      std::vector<double> vec_nodes(_mesh.NN(), 0.0);
      for (int idxn = 0; idxn < _mesh.NN(); ++idxn) {
        vec_nodes[idxn] = _mesh.NodeRadius(idx, idxn);
      }
      auto pleg =
          Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());

      // loop through the nodes and evaluate
      for (int idxq = 0; idxq < _mesh.NN(); ++idxq) {
        auto idx_u = this->LtG_S(0, idx, idxq);
        // receiver vector
        vec_receiver(0, idx_u) = pleg(idxq, rad_r);
      }
    }
  }

  return vec_receiver;
};

auto
specsem::Receiver_Elements(InputParameters &param) const {
  std::vector<int> receiver_elems;
  double rad_receiver =
      _mesh.PR() - 1000.0 * param.receiver_depth() / _length_norm;

  for (int idx = 0; idx < _mesh.NE(); ++idx) {
    if ((_mesh.ELR(idx) <= rad_receiver) && (_mesh.EUR(idx) >= rad_receiver)) {
      receiver_elems.push_back(idx);
    }
  }

  return receiver_elems;
};

Eigen::MatrixXcd
specsem::RV_FULL(InputParameters &param, int idxl) {
  auto nrec = param.num_receivers();
  using MATRIX = Eigen::MatrixXcd;
  MATRIX vec_receiver = MATRIX::Zero(3 * nrec, 2 * idxl + 1);

  using namespace GSHTrans;
  int maxn = 1;
  auto i1 = std::complex<double>(0.0, 1.0);
  auto deg2rad = EIGEN_PI / 180.0;

  // vector of theta values:
  std::vector<double> vec_theta_r(nrec, 0.0);
  std::vector<double> vec_phi_r(nrec, 0.0);
  for (int idx = 0; idx < nrec; ++idx) {
    auto rec = param.receivers()[idx];
    vec_theta_r[idx] = (90.0 - rec.first) * deg2rad;
    vec_phi_r[idx] = rec.second * deg2rad;
  }
  auto wigdmat_all = Wigner<double, Ortho, All, All, Multiple, ColumnMajor>(
      idxl, idxl, maxn, vec_theta_r);
  auto ylmn_all = [&wigdmat_all, &vec_phi_r](int l, int m, int N, int idxr) {
    auto dl = wigdmat_all[N, idxr];
    auto tmp = dl[l, m];
    auto ylm = tmp * std::exp(std::complex<double>(0.0, m * vec_phi_r[idxr]));
    return ylm;
  };

  // loop through receivers
  for (int idxr = 0; idxr < nrec; ++idxr) {
    for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
      // spherical harmonic values
      Complex yl0 = ylmn_all(idxl, idxm, 0, idxr);
      Complex ylm = ylmn_all(idxl, idxm, -1, idxr);
      Complex ylp = ylmn_all(idxl, idxm, 1, idxr);

      // adding to receiver vector
      vec_receiver(3 * idxr, idxm + idxl) = yl0;
      vec_receiver(3 * idxr + 1, idxm + idxl) = ylp - ylm;
      vec_receiver(3 * idxr + 2, idxm + idxl) = -i1 * (ylm + ylp);
    }
  }

  return vec_receiver;
};

Eigen::MatrixXcd
specsem::RV_FULL_T(InputParameters &param, int idxl) {
  auto nrec = param.num_receivers();
  using MATRIX = Eigen::MatrixXcd;
  MATRIX vec_receiver = MATRIX::Zero(3 * nrec, 2 * idxl + 1);

  using namespace GSHTrans;
  auto i1 = std::complex<double>(0.0, 1.0);
  auto deg2rad = EIGEN_PI / 180.0;

  // vector of theta values:
  std::vector<double> vec_theta_r(nrec, 0.0);
  std::vector<double> vec_phi_r(nrec, 0.0);
  for (int idx = 0; idx < nrec; ++idx) {
    auto rec = param.receivers()[idx];
    vec_theta_r[idx] = (90.0 - rec.first) * deg2rad;
    vec_phi_r[idx] = rec.second * deg2rad;
  }
  auto wigdmat_all = Wigner<double, Ortho, All, All, Multiple, ColumnMajor>(
      idxl, idxl, 1, vec_theta_r);
  auto ylmn_all = [&wigdmat_all, &vec_phi_r](int l, int m, int N, int idxr) {
    auto dl = wigdmat_all[N, idxr];
    auto tmp = dl[l, m];
    auto ylm = tmp * std::exp(std::complex<double>(0.0, m * vec_phi_r[idxr]));
    return ylm;
  };

  // loop through receivers
  for (int idxr = 0; idxr < nrec; ++idxr) {
    for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
      // spherical harmonic values
      Complex ylm = ylmn_all(idxl, idxm, -1, idxr);
      Complex ylp = ylmn_all(idxl, idxm, 1, idxr);

      // adding to receiver vector
      vec_receiver(3 * idxr + 1, idxm + idxl) = i1 * (ylm + ylp);
      vec_receiver(3 * idxr + 2, idxm + idxl) = ylp - ylm;
    }
  }

  return vec_receiver;
};

Eigen::MatrixXcd
specsem::RV_VAL_Z(InputParameters &param, int idxl, int idxr) {
  auto rec = param.receivers()[idxr];
  double theta_r = (90.0 - rec.first) * EIGEN_PI / 180.0;
  double phi_r = rec.second * EIGEN_PI / 180.0;

  // create the receiver vector
  std::size_t fcols = 2 * idxl + 1;
  Eigen::MatrixXcd vec_receiver = Eigen::MatrixXcd::Zero(fcols, 1);

  //////////////////////////////
  // wigner matrix for evaluation
  using namespace GSHTrans;
  int maxn = 2;
  if (maxn > _lmax) {
    maxn = _lmax;
  }
  auto wigdmat = Wigner<double, Ortho, All, All, Single, ColumnMajor>(
      _lmax, _lmax, maxn, theta_r);
  auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
    auto dl = wigdmat[N];
    auto tmp = dl[l, m];
    auto ylm = tmp * std::exp(std::complex<double>(0.0, m * phi));
    return ylm;
  };

  for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
    // spherical harmonic values
    Complex yl0 = ylmn(idxl, idxm, 0, phi_r);
    vec_receiver(idxm + idxl, 0) = yl0;
  }
  // Complex yl0 = ylmn(idxl, 0, 0, phi_r);
  // vec_receiver(idxl, 0) = yl0;

  return vec_receiver;
};

Eigen::MatrixXcd
specsem::RV_Z_R(InputParameters &param, int idxr) {
  auto rec = param.receivers()[idxr];
  double theta_r = (90.0 - rec.first) * EIGEN_PI / 180.0;
  double phi_r = rec.second * EIGEN_PI / 180.0;
  double rad_r = _mesh.PR() - param.receiver_depth() * 1000.0 / _length_norm;

  // create the receiver vector
  std::size_t flen = this->LtG_R(1, _mesh.NE() - 1, _mesh.NN() - 1) + 1;
  // std::size_t fcols = 2 * idxl + 1;
  Eigen::MatrixXcd vec_receiver = Eigen::MatrixXcd::Zero(flen, 1);

  //////////////////////////////
  // wigner matrix for evaluation
  using namespace GSHTrans;
  int maxn = 2;
  if (maxn > _lmax) {
    maxn = _lmax;
  }
  auto wigdmat = Wigner<double, Ortho, All, All, Single, ColumnMajor>(
      _lmax, _lmax, maxn, theta_r);
  auto wigdmat2 =
      Wigner<double, Ortho, All, All, Single, ColumnMajor>(0, 0, 0, 0.0);
  auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
    auto dl = wigdmat[N];
    auto tmp = dl[l, m];
    auto ylm = tmp * std::exp(std::complex<double>(0.0, m * phi));
    return ylm;
  };

  // int elem_idx = 0;
  bool evaluated = false;
  for (int idx = 0; idx < _mesh.NE(); ++idx) {
    if ((_mesh.ELR(idx) <= rad_r) && (_mesh.EUR(idx) >= rad_r) &&
        (!evaluated)) {
      evaluated = true;

      // get the Lagrange polynomial for this element
      std::vector<double> vec_nodes(_mesh.NN(), 0.0);
      for (int idxn = 0; idxn < _mesh.NN(); ++idxn) {
        vec_nodes[idxn] = _mesh.NodeRadius(idx, idxn);
      }
      auto pleg =
          Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());

      // loop through the nodes and evaluate
      for (int idxq = 0; idxq < _mesh.NN(); ++idxq) {
        auto idx_u = this->LtG_R(0, idx, idxq);
        // for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
        // spherical harmonic values
        // Complex yl0 = ylmn(0, 0, 0, phi_r);
        Complex yl0 = wigdmat2[0][0, 0];

        // receiver vector
        vec_receiver(idx_u, 0) = yl0 * pleg(idxq, rad_r);
      }
    }
  }

  return vec_receiver;
};

Eigen::MatrixXcd
specsem::RV_RED_Z_R(InputParameters &param) {
  double rad_r = _mesh.PR() - param.receiver_depth() * 1000.0 / _length_norm;
  auto rec_elems = this->Receiver_Elements(param);
  auto lowidx = this->LtG_R(0, rec_elems[0], 0);
  auto upidx = this->LtG_R(1, rec_elems.back(), _mesh.NN() - 1);
  int lenidx = upidx - lowidx + 1;
  // std::cout << "RV_RED_Z_R: lowidx: " << lowidx << ", upidx: " << upidx
  //           << ", lenidx: " << lenidx << "\n";
  Eigen::MatrixXcd vec_receiver = Eigen::MatrixXcd::Zero(lenidx, 1);
  //////////////////////////////
  // wigner matrix for evaluation
  using namespace GSHTrans;
  auto wigdmat =
      Wigner<double, Ortho, All, All, Single, ColumnMajor>(0, 0, 0, 0.0);

  bool evaluated = false;
  for (int idx = 0; idx < _mesh.NE(); ++idx) {
    if ((_mesh.ELR(idx) <= rad_r) && (_mesh.EUR(idx) >= rad_r) &&
        (!evaluated)) {
      evaluated = true;

      // get the Lagrange polynomial for this element
      std::vector<double> vec_nodes(_mesh.NN(), 0.0);
      for (int idxn = 0; idxn < _mesh.NN(); ++idxn) {
        vec_nodes[idxn] = _mesh.NodeRadius(idx, idxn);
      }
      auto pleg =
          Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());
      Complex yl0 = wigdmat[0][0, 0];
      // loop through the nodes and evaluate
      for (int idxq = 0; idxq < _mesh.NN(); ++idxq) {
        auto idx_u = this->LtG_R(0, idx, idxq) - lowidx;

        // receiver vector
        vec_receiver(idx_u, 0) = yl0 * pleg(idxq, rad_r);
      }
    }
  }

  return vec_receiver;
};

Eigen::MatrixXcd
specsem::RV_BASE_THETA(InputParameters &param, int idxl, int idxr) {

  double rad_r = _mesh.PR() - param.receiver_depth() * 1000.0 / _length_norm;

  auto k =
      std::sqrt(static_cast<double>(idxl) * (static_cast<double>(idxl) + 1.0));
  // create the receiver vector
  std::size_t flen = this->LtG_S(2, _mesh.NE() - 1, _mesh.NN() - 1) + 1;
  Eigen::MatrixXcd vec_receiver = Eigen::MatrixXcd::Zero(1, flen);

  // int elem_idx = 0;
  bool evaluated = false;
  for (int idx = 0; idx < _mesh.NE(); ++idx) {
    if ((_mesh.ELR(idx) <= rad_r) && (_mesh.EUR(idx) >= rad_r) &&
        (!evaluated)) {
      evaluated = true;

      // get the Lagrange polynomial for this element
      std::vector<double> vec_nodes(_mesh.NN(), 0.0);
      for (int idxn = 0; idxn < _mesh.NN(); ++idxn) {
        vec_nodes[idxn] = _mesh.NodeRadius(idx, idxn);
      }
      auto pleg =
          Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());

      // loop through the nodes and evaluate
      for (int idxq = 0; idxq < _mesh.NN(); ++idxq) {
        auto idx_v = this->LtG_S(1, idx, idxq);

        // receiver vector
        vec_receiver(0, idx_v) = k / 2.0 * pleg(idxq, rad_r);
      }
    }
  }

  return vec_receiver;
};

Eigen::MatrixXcd
specsem::RV_BASE_THETA_T(InputParameters &param, int idxl, int idxr) {

  double rad_r = _mesh.PR() - param.receiver_depth() * 1000.0 / _length_norm;

  auto k =
      std::sqrt(static_cast<double>(idxl) * (static_cast<double>(idxl) + 1.0));
  // create the receiver vector
  std::size_t flen = this->LtG_T(_mesh.NE() - 1, _mesh.NN() - 1) + 1;
  Eigen::MatrixXcd vec_receiver = Eigen::MatrixXcd::Zero(1, flen);
  // std::complex<double> i1 = std::complex<double>(0.0, 1.0);
  // int elem_idx = 0;
  bool evaluated = false;
  for (int idx = 0; idx < _mesh.NE(); ++idx) {
    if ((_mesh.ELR(idx) <= rad_r) && (_mesh.EUR(idx) >= rad_r) &&
        (!evaluated)) {
      evaluated = true;

      // get the Lagrange polynomial for this element
      std::vector<double> vec_nodes(_mesh.NN(), 0.0);
      for (int idxn = 0; idxn < _mesh.NN(); ++idxn) {
        vec_nodes[idxn] = _mesh.NodeRadius(idx, idxn);
      }
      auto pleg =
          Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());

      // loop through the nodes and evaluate
      for (int idxq = 0; idxq < _mesh.NN(); ++idxq) {
        auto idx_v = this->LtG_T(idx, idxq);

        // receiver vector
        vec_receiver(0, idx_v) = 0.5 * pleg(idxq, rad_r);
      }
    }
  }

  return vec_receiver;
};

Eigen::MatrixXcd
specsem::RV_VAL_THETA(InputParameters &param, int idxl, int idxr) {
  auto rec = param.receivers()[idxr];
  double theta_r = (90.0 - rec.first) * EIGEN_PI / 180.0;
  double phi_r = rec.second * EIGEN_PI / 180.0;
  std::size_t fcols = 2 * idxl + 1;
  Eigen::MatrixXcd vec_receiver = Eigen::MatrixXcd::Zero(fcols, 1);

  //////////////////////////////
  // wigner matrix for evaluation
  using namespace GSHTrans;
  int maxn = 2;
  if (maxn > _lmax) {
    maxn = _lmax;
  }
  auto wigdmat = Wigner<double, Ortho, All, All, Single, ColumnMajor>(
      _lmax, _lmax, maxn, theta_r);
  auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
    auto dl = wigdmat[N];
    auto tmp = dl[l, m];
    auto ylm = tmp * std::exp(std::complex<double>(0.0, m * phi));
    return ylm;
  };

  for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
    // spherical harmonic values
    Complex ylm = ylmn(idxl, idxm, -1, phi_r);
    Complex ylp = ylmn(idxl, idxm, 1, phi_r);
    vec_receiver(idxm + idxl, 0) = ylm - ylp;
  }

  return vec_receiver;
};

Eigen::MatrixXcd
specsem::RV_VAL_THETA_T(InputParameters &param, int idxl, int idxr) {
  auto rec = param.receivers()[idxr];
  double theta_r = (90.0 - rec.first) * EIGEN_PI / 180.0;
  double phi_r = rec.second * EIGEN_PI / 180.0;
  std::size_t fcols = 2 * idxl + 1;
  Eigen::MatrixXcd vec_receiver = Eigen::MatrixXcd::Zero(fcols, 1);

  //////////////////////////////
  // wigner matrix for evaluation
  using namespace GSHTrans;
  int maxn = 2;
  if (maxn > _lmax) {
    maxn = _lmax;
  }
  auto wigdmat = Wigner<double, Ortho, All, All, Single, ColumnMajor>(
      _lmax, _lmax, maxn, theta_r);
  auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
    auto dl = wigdmat[N];
    auto tmp = dl[l, m];
    auto ylm = tmp * std::exp(std::complex<double>(0.0, m * phi));
    return ylm;
  };
  auto i1 = std::complex<double>(0.0, 1.0);

  for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
    // spherical harmonic values
    Complex ylm = ylmn(idxl, idxm, -1, phi_r);
    Complex ylp = ylmn(idxl, idxm, 1, phi_r);
    vec_receiver(idxm + idxl, 0) = -i1 * (ylm + ylp);
  }

  return vec_receiver;
};

Eigen::MatrixXcd
specsem::RV_THETA_T(InputParameters &param, int idxl, int idxr) {
  auto rec = param.receivers()[idxr];
  double theta_r = (90.0 - rec.first) * EIGEN_PI / 180.0;
  double phi_r = rec.second * EIGEN_PI / 180.0;
  double rad_r = _mesh.PR() - param.receiver_depth() * 1000.0 / _length_norm;

  auto k =
      std::sqrt(static_cast<double>(idxl) * (static_cast<double>(idxl) + 1.0));
  // create the receiver vector
  std::size_t flen = this->LtG_T(_eu - 1, _mesh.NN() - 1) + 1;
  std::size_t fcols = 2 * idxl + 1;
  Eigen::MatrixXcd vec_receiver = Eigen::MatrixXcd::Zero(flen, fcols);

  //////////////////////////////
  // wigner matrix for evaluation
  using namespace GSHTrans;
  int maxn = 2;
  if (maxn > _lmax) {
    maxn = _lmax;
  }
  auto wigdmat = Wigner<double, Ortho, All, All, Single, ColumnMajor>(
      _lmax, _lmax, maxn, theta_r);
  auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
    auto dl = wigdmat[N];
    auto tmp = dl[l, m];
    std::complex<double> ylm =
        tmp * std::exp(std::complex<double>(0.0, m * phi));
    return ylm;
  };
  Complex i1 = std::complex<double>(0.0, 1.0);
  // int elem_idx = 0;
  bool evaluated = false;
  for (int idx = 0; idx < _mesh.NE(); ++idx) {
    if ((_mesh.ELR(idx) <= rad_r) && (_mesh.EUR(idx) >= rad_r) &&
        (!evaluated)) {
      evaluated = true;

      // get the Lagrange polynomial for this element
      std::vector<double> vec_nodes(_mesh.NN(), 0.0);
      for (int idxn = 0; idxn < _mesh.NN(); ++idxn) {
        vec_nodes[idxn] = _mesh.NodeRadius(idx, idxn);
      }
      auto pleg =
          Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());

      // loop through the nodes and evaluate
      for (int idxq = 0; idxq < _mesh.NN(); ++idxq) {
        auto idx_v = this->LtG_T(idx, idxq);
        for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
          // spherical harmonic values
          Complex ylm = ylmn(idxl, idxm, -1, phi_r);
          Complex ylp = ylmn(idxl, idxm, 1, phi_r);

          // receiver vector
          vec_receiver(idx_v, idxm + idxl) =
              -i1 * pleg(idxq, rad_r) * (ylp + ylm) / 2.0;
        }
      }
    }
  }

  return vec_receiver;
};

Eigen::MatrixXcd
specsem::RV_BASE_PHI_T(InputParameters &param, int idxl, int idxr) {
  return this->RV_BASE_THETA_T(param, idxl, idxr);
};

Eigen::MatrixXcd
specsem::RV_VAL_PHI(InputParameters &param, int idxl, int idxr) {
  auto rec = param.receivers()[idxr];
  double theta_r = (90.0 - rec.first) * EIGEN_PI / 180.0;
  double phi_r = rec.second * EIGEN_PI / 180.0;
  std::size_t fcols = 2 * idxl + 1;
  Eigen::MatrixXcd vec_receiver = Eigen::MatrixXcd::Zero(fcols, 1);

  //////////////////////////////
  // wigner matrix for evaluation
  using namespace GSHTrans;
  int maxn = 2;
  if (maxn > _lmax) {
    maxn = _lmax;
  }
  auto wigdmat = Wigner<double, Ortho, All, All, Single, ColumnMajor>(
      _lmax, _lmax, maxn, theta_r);
  auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
    auto dl = wigdmat[N];
    auto tmp = dl[l, m];
    auto ylm = tmp * std::exp(std::complex<double>(0.0, m * phi));
    return ylm;
  };
  Complex i1 = std::complex<double>(0.0, 1.0);
  for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
    // spherical harmonic values
    Complex ylm = ylmn(idxl, idxm, -1, phi_r);
    Complex ylp = ylmn(idxl, idxm, 1, phi_r);
    vec_receiver(idxm + idxl, 0) = -i1 * (ylm + ylp);
  }

  return vec_receiver;
};

Eigen::MatrixXcd
specsem::RV_VAL_PHI_T(InputParameters &param, int idxl, int idxr) {
  auto rec = param.receivers()[idxr];
  double theta_r = (90.0 - rec.first) * EIGEN_PI / 180.0;
  double phi_r = rec.second * EIGEN_PI / 180.0;
  std::size_t fcols = 2 * idxl + 1;
  Eigen::MatrixXcd vec_receiver = Eigen::MatrixXcd::Zero(fcols, 1);

  //////////////////////////////
  // wigner matrix for evaluation
  using namespace GSHTrans;
  int maxn = 2;
  if (maxn > _lmax) {
    maxn = _lmax;
  }
  auto wigdmat = Wigner<double, Ortho, All, All, Single, ColumnMajor>(
      _lmax, _lmax, maxn, theta_r);
  auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
    auto dl = wigdmat[N];
    auto tmp = dl[l, m];
    auto ylm = tmp * std::exp(std::complex<double>(0.0, m * phi));
    return ylm;
  };

  for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
    // spherical harmonic values
    Complex ylm = ylmn(idxl, idxm, -1, phi_r);
    Complex ylp = ylmn(idxl, idxm, 1, phi_r);
    vec_receiver(idxm + idxl, 0) = (ylp - ylm);
  }

  return vec_receiver;
};

Eigen::MatrixXcd
specsem::RV_PHI_T(InputParameters &param, int idxl, int idxr) {
  auto rec = param.receivers()[idxr];
  double theta_r = (90.0 - rec.first) * EIGEN_PI / 180.0;
  double phi_r = rec.second * EIGEN_PI / 180.0;
  double rad_r = _mesh.PR() - param.receiver_depth() * 1000.0 / _length_norm;

  auto k =
      std::sqrt(static_cast<double>(idxl) * (static_cast<double>(idxl) + 1.0));
  // create the receiver vector
  std::size_t flen = this->LtG_T(_eu - 1, _mesh.NN() - 1) + 1;
  std::size_t fcols = 2 * idxl + 1;
  Eigen::MatrixXcd vec_receiver = Eigen::MatrixXcd::Zero(flen, fcols);

  //////////////////////////////
  // wigner matrix for evaluation
  using namespace GSHTrans;
  int maxn = 2;
  if (maxn > _lmax) {
    maxn = _lmax;
  }
  auto wigdmat = Wigner<double, Ortho, All, All, Single, ColumnMajor>(
      _lmax, _lmax, maxn, theta_r);
  auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
    auto dl = wigdmat[N];
    auto tmp = dl[l, m];
    std::complex<double> ylm =
        tmp * std::exp(std::complex<double>(0.0, m * phi));
    return ylm;
  };
  Complex i1 = std::complex<double>(0.0, 1.0);
  // int elem_idx = 0;
  bool evaluated = false;
  for (int idx = 0; idx < _mesh.NE(); ++idx) {
    if ((_mesh.ELR(idx) <= rad_r) && (_mesh.EUR(idx) >= rad_r) &&
        (!evaluated)) {
      evaluated = true;

      // get the Lagrange polynomial for this element
      std::vector<double> vec_nodes(_mesh.NN(), 0.0);
      for (int idxn = 0; idxn < _mesh.NN(); ++idxn) {
        vec_nodes[idxn] = _mesh.NodeRadius(idx, idxn);
      }
      auto pleg =
          Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());

      // loop through the nodes and evaluate
      for (int idxq = 0; idxq < _mesh.NN(); ++idxq) {
        auto idx_v = this->LtG_T(idx, idxq);
        for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
          // spherical harmonic values
          Complex ylm = ylmn(idxl, idxm, -1, phi_r);
          Complex ylp = ylmn(idxl, idxm, 1, phi_r);

          // receiver vector
          vec_receiver(idx_v, idxm + idxl) =
              pleg(idxq, rad_r) * (ylp - ylm) / 2.0;
        }
      }
    }
  }

  return vec_receiver;
};

Eigen::MatrixXcd
specsem::RV_BASE_FULL(InputParameters &param, int idxl) {
  // Implementation of RV_BASE_FULL
  auto rec_elems = this->Receiver_Elements(param);
  auto lowidx = this->LtG_S(0, rec_elems[0], 0);
  auto upidx = this->LtG_S(1, rec_elems.back(), this->_NQ - 1);
  int lenidx = upidx - lowidx + 1;
  auto nrec = param.num_receivers();
  Eigen::MatrixXcd mat_base = Eigen::MatrixXcd::Zero(3 * nrec, lenidx);
  double rad_r =
      _mesh.PR() - param.receiver_depth() * 1000.0 / this->_length_norm;
  double k = std::sqrt(1.0 * idxl * (idxl + 1.0));

  // loop through elements containing source
  for (int idx = rec_elems[0]; idx < rec_elems.back() + 1; ++idx) {
    // get the Lagrange polynomial for this element
    std::vector<double> vec_nodes(_mesh.NN(), 0.0);
    for (int idxn = 0; idxn < _mesh.NN(); ++idxn) {
      vec_nodes[idxn] = _mesh.NodeRadius(idx, idxn);
    }
    auto pleg =
        Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());

    // loop through the nodes and evaluate
    for (int idxq = 0; idxq < _mesh.NN(); ++idxq) {
      // indices
      auto idx_u = this->LtG_S(0, idx, idxq) - lowidx;
      auto idx_v = this->LtG_S(1, idx, idxq) - lowidx;

      // values
      auto zv = pleg(idxq, rad_r);
      auto tv = k / 2.0 * zv;

      // loop through receivers
      for (int idxr = 0; idxr < nrec; ++idxr) {
        mat_base(3 * idxr, idx_u) = zv;
        mat_base(3 * idxr + 1, idx_v) = tv;
        mat_base(3 * idxr + 2, idx_v) = tv;
      }
    }
  }

  return mat_base;
};

Eigen::MatrixXcd
specsem::RV_BASE_FULL_T(InputParameters &param, int idxl) {
  // Implementation of RV_BASE_FULL
  auto rec_elems = this->Receiver_Elements(param);
  auto lowidx = this->LtG_T(rec_elems[0], 0);
  auto upidx = this->LtG_T(rec_elems.back(), this->_NQ - 1);
  int lenidx = upidx - lowidx + 1;
  auto nrec = param.num_receivers();
  Eigen::MatrixXcd mat_base = Eigen::MatrixXcd::Zero(3 * nrec, lenidx);
  double rad_r =
      _mesh.PR() - param.receiver_depth() * 1000.0 / this->_length_norm;
  double k = std::sqrt(1.0 * idxl * (idxl + 1.0));

  // loop through elements containing source
  for (int idx = rec_elems[0]; idx < rec_elems.back() + 1; ++idx) {
    // get the Lagrange polynomial for this element
    std::vector<double> vec_nodes(_mesh.NN(), 0.0);
    for (int idxn = 0; idxn < _mesh.NN(); ++idxn) {
      vec_nodes[idxn] = _mesh.NodeRadius(idx, idxn);
    }
    auto pleg =
        Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());

    // loop through the nodes and evaluate
    for (int idxq = 0; idxq < _mesh.NN(); ++idxq) {
      // indices
      auto idx_v = this->LtG_T(idx, idxq) - lowidx;

      // values
      auto zv = pleg(idxq, rad_r);
      auto tv = k / 2.0 * zv;

      // loop through receivers
      for (int idxr = 0; idxr < nrec; ++idxr) {
        mat_base(3 * idxr + 1, idx_v) = tv;
        mat_base(3 * idxr + 2, idx_v) = tv;
      }
    }
  }

  return mat_base;
};

}   // namespace Full1D

#endif   // SPECTRA_MASTER_IMPL_H
