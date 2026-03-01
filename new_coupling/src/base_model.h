#ifndef NEW_COUPLING_SRC_BASE_H
#define NEW_COUPLING_SRC_BASE_H

// Header files
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <chrono>
#include <iomanip>
#include <Interpolation/All>
#include <EarthMesh/All>
// #include <Gravitational_Field/All>
#include <PlanetaryModel/All>
// #include <Gravitational_Field/Test>
// #include <Gravitational_Field/Timer>

namespace ModeCoupling {
class basetest {
public:
  basetest() {};
  basetest(double inpd) { _bvec = std::vector<double>(1, inpd); };
  auto veceval() { return _bvec[0]; };

private:
  std::vector<double> _bvec;
};

class mineos_eigenfunction {
public:
  mineos_eigenfunction() {};
  mineos_eigenfunction(int n, int l, double q, double gv, double w,
                       std::vector<double> &u, std::vector<double> &up,
                       std::vector<double> &v, std::vector<double> &vp,
                       std::vector<double> &p, std::vector<double> &pp)
      : _n(n), _l(l), _q(q), _gv(gv), _w(w), _u(u), _up(up), _v(v), _vp(vp),
        _p(p), _pp(pp) {};
  mineos_eigenfunction(int n, int l, double q, double gv, double w,
                       std::vector<double> &vec_longinput)
      : _n(n), _l(l), _q(q), _gv(gv), _w(w) {
    int lnum = vec_longinput.size() / 6;
    auto myit = vec_longinput.begin();
    _u = std::vector<double>(myit, myit + lnum);
    _up = std::vector<double>(myit + lnum, myit + 2 * lnum);
    _v = std::vector<double>(myit + 2 * lnum, myit + 3 * lnum);
    _vp = std::vector<double>(myit + 3 * lnum, myit + 4 * lnum);
    _p = std::vector<double>(myit + 4 * lnum, myit + 5 * lnum);
    _pp = std::vector<double>(myit + 5 * lnum, vec_longinput.end());
  };

  auto n() { return _n; };
  auto l() { return _l; };
  auto q() { return _q; };
  auto gv() { return _gv; };
  auto w() { return _w; };

  auto u() { return _u; };
  auto up() { return _up; };
  auto v() { return _v; };
  auto vp() { return _vp; };
  auto p() { return _p; };
  auto pp() { return _pp; };

  //  auto ku() { return _ku; };
  //  auto kr() { return _kr; };
  //  auto kd() { return _kd; };

private:
  int _n, _l;
  double _q, _gv, _w;
  std::vector<double> _u, _up, _v, _vp, _p, _pp;   //, _ku, _kr, _kd;
  //  std::vector<double> u = std::vector<double>(185);
  //  std::vector<double> up = std::vector<double>(185);
  //  std::vector<double> v = std::vector<double>(185);
  //  std::vector<double> vp = std::vector<double>(185);
  //  std::vector<double> p = std::vector<double>(185);
  //  std::vector<double> pp = std::vector<double>(185);
  //  std::vector<double> ku = std::vector<double>(185);
  //  std::vector<double> kr = std::vector<double>(185);
  //  std::vector<double> kd = std::vector<double>(185);

  // std::vector<double> u, up, v, vp, p, pp;
};

class mineos_eigenfunction_continuous {
public:
  using myvector = std::vector<double>;
  using myiter = myvector::iterator;
  using InterpA = Interpolation::CubicSpline<myiter, myiter>;
  // using InterpA = Interpolation::Akima<myiter, myiter>;

  mineos_eigenfunction_continuous() {};
  mineos_eigenfunction_continuous(
      int n, int l, double q, double gv, double w, std::vector<double> &u,
      std::vector<double> &up, std::vector<double> &v, std::vector<double> &vp,
      std::vector<double> &p, std::vector<double> &pp)
      : _n(n), _l(l), _q(q), _gv(gv), _w(w), _u(u), _up(up), _v(v), _vp(vp),
        _p(p), _pp(pp) {};

  mineos_eigenfunction_continuous(
      int n, int l, double q, double gv, double w,
      std::vector<double> &vec_longinput,
      EarthModels::ModelInput<double, int> &inp_model,
      std::vector<std::vector<double>> &vec_radii)
      : _n(n), _l(l), _q(q), _gv(gv), _w(w), _earthmodel(inp_model) {
    int lnum = vec_longinput.size() / 6;
    auto myit = vec_longinput.begin();
    _u = std::vector<double>(myit, myit + lnum);
    double mfact = 2.0 * 3.14159265358979 * w;
    // for (auto &idx : _u) {
    //    idx *= w * 2.0 * twopi;
    // }

    _up = std::vector<double>(myit + lnum, myit + 2 * lnum);
    _v = std::vector<double>(myit + 2 * lnum, myit + 3 * lnum);
    _vp = std::vector<double>(myit + 3 * lnum, myit + 4 * lnum);
    _p = std::vector<double>(myit + 4 * lnum, myit + 5 * lnum);
    _pp = std::vector<double>(myit + 5 * lnum, vec_longinput.end());
    // std::transform(_u.begin(), _u.end(), _u.begin(),
    //                [mfact](double c) { return c * mfact; });
    // std::transform(_up.begin(), _up.end(), _up.begin(),
    //                [mfact](double c) { return c * mfact; });
    // std::transform(_v.begin(), _v.end(), _v.begin(),
    //                [mfact](double c) { return c * mfact; });
    // std::transform(_vp.begin(), _vp.end(), _vp.begin(),
    //                [mfact](double c) { return c * mfact; });
    // std::transform(_p.begin(), _p.end(), _p.begin(),
    // [mfact](double c) { return c * mfact; });
    // std::transform(_pp.begin(), _pp.end(), _pp.begin(),
    // [mfact](double c) { return c * mfact; });
    // std::cout << _u.size() << "\n";

    for (int idx = 0; idx < inp_model.NumberOfLayers(); ++idx) {
      auto itr1 = vec_radii[idx].begin();
      auto itr2 = vec_radii[idx].end();

      auto itu1 = _u.begin() + inp_model.LayerLowerIndex(idx);
      auto itu2 = _u.begin() + inp_model.LayerUpperIndex(idx) + 1;
      _ul.push_back(std::vector<double>(itu1, itu2));
      auto tmpu = InterpA(itr1, itr2, _ul[idx].begin());
      _ucont.push_back(tmpu);

      auto itup1 = _up.begin() + inp_model.LayerLowerIndex(idx);
      auto itup2 = _up.begin() + inp_model.LayerUpperIndex(idx) + 1;
      _upl.push_back(std::vector<double>(itup1, itup2));
      auto tmpup = InterpA(itr1, itr2, _upl[idx].begin());
      _upcont.push_back(tmpup);

      auto itv1 = _v.begin() + inp_model.LayerLowerIndex(idx);
      auto itv2 = _v.begin() + inp_model.LayerUpperIndex(idx) + 1;
      _vl.push_back(std::vector<double>(itv1, itv2));
      auto tmpv = InterpA(itr1, itr2, _vl[idx].begin());
      _vcont.push_back(tmpv);

      auto itvp1 = _vp.begin() + inp_model.LayerLowerIndex(idx);
      auto itvp2 = _vp.begin() + inp_model.LayerUpperIndex(idx) + 1;
      _vpl.push_back(std::vector<double>(itvp1, itvp2));
      auto tmpvp = InterpA(itr1, itr2, _vpl[idx].begin());
      _vpcont.push_back(tmpvp);

      auto itp1 = _p.begin() + inp_model.LayerLowerIndex(idx);
      auto itp2 = _p.begin() + inp_model.LayerUpperIndex(idx) + 1;
      _pl.push_back(std::vector<double>(itp1, itp2));
      auto tmpp = InterpA(itr1, itr2, _pl[idx].begin());
      _pcont.push_back(tmpp);

      auto itpp1 = _pp.begin() + inp_model.LayerLowerIndex(idx);
      auto itpp2 = _pp.begin() + inp_model.LayerUpperIndex(idx) + 1;
      _ppl.push_back(std::vector<double>(itpp1, itpp2));
      auto tmppp = InterpA(itr1, itr2, _ppl[idx].begin());
      _ppcont.push_back(tmppp);
    };
  }

  auto n() { return _n; };
  auto l() { return _l; };
  auto q() { return _q; };
  auto gv() { return _gv; };
  auto w() { return _w; };

  // auto radii() { return _radii; };
  auto u() { return _u; };
  auto up() { return _up; };
  auto v() { return _v; };
  auto vp() { return _vp; };
  auto p() { return _p; };
  auto pp() { return _pp; };

  // return continuous u
  auto u(int i) { return _ucont[i]; };
  auto up(int i) { return _upcont[i]; };
  auto v(int i) { return _vcont[i]; };
  auto vp(int i) { return _vpcont[i]; };
  auto p(int i) { return _pcont[i]; };
  auto pp(int i) { return _ppcont[i]; };

  // auto u_pointer(int i) { return _u[i].begin(); };

  //  auto ku() { return _ku; };
  //  auto kr() { return _kr; };
  //  auto kd() { return _kd; };

private:
  int _n, _l;
  double _q, _gv, _w;
  std::vector<double> _u, _up, _v, _vp, _p, _pp;   //, _ku, _kr, _kd;
  std::vector<std::vector<double>> _ul, _upl, _vl, _vpl, _pl,
      _ppl;   //, _ku, _kr, _kd;
  std::vector<InterpA> _ucont, _upcont, _vcont, _vpcont, _pcont, _ppcont;
  EarthModels::ModelInput<double, int> _earthmodel;

  //  std::vector<double> u = std::vector<double>(185);
  //  std::vector<double> up = std::vector<double>(185);
  //  std::vector<double> v = std::vector<double>(185);
  //  std::vector<double> vp = std::vector<double>(185);
  //  std::vector<double> p = std::vector<double>(185);
  //  std::vector<double> pp = std::vector<double>(185);
  //  std::vector<double> ku = std::vector<double>(185);
  //  std::vector<double> kr = std::vector<double>(185);
  //  std::vector<double> kd = std::vector<double>(185);

  // std::vector<double> u, up, v, vp, p, pp;
};

class eigfunc_complete {

  std::vector<mineos_eigenfunction> normal_modes;
};

class modecataloguecontinuous {

public:
  modecataloguecontinuous() {};
  template <class MODEL> modecataloguecontinuous(std::string &, MODEL &);

  auto allmodes() { return _mode_storage; };
  auto isincatalogue(int n, int l) {
    bool boolout = false;
    for (int idx = 0; idx < _vec_nl.size(); ++idx) {
      //  int idxoutput;

      if (n == _vec_nl[idx][0]) {
        if (l == _vec_nl[idx][1]) {
          boolout = true;
          break;
        }
      }
    }
    return boolout;
  };
  auto singlemode(int i) {
    assert(i < _mode_storage.size() && "Not that many modes!");
    return _mode_storage[i];
  };
  mineos_eigenfunction_continuous &singlemodep(int i) {
    assert(i < _mode_storage.size() && "Not that many modes!");
    return _mode_storage[i];
  };
  auto NumberOfModes() { return _mode_storage.size(); };

private:
  std::vector<mineos_eigenfunction_continuous> _mode_storage;
  std::vector<std::vector<int>> _vec_nl;
  std::vector<std::vector<double>> _vec_radii;
};

template <class MODEL>
modecataloguecontinuous::modecataloguecontinuous(std::string &pathtofile,
                                                 MODEL &inp_model) {
  _vec_radii = inp_model.LayerRadii();
  // open the file for input, in binary form
  std::ifstream infile(pathtofile, std::ifstream::binary);

  // check opened correctly
  if (!infile) {
    std::cout << "Cannot open file!" << std::endl;
  }

  // find the length of the file
  infile.seekg(0, std::ios::end);
  std::size_t file_size = infile.tellg();
  std::cout << "Size of the file is " << file_size << " bytes" << std::endl;
  infile.seekg(0, infile.beg);   // reset back to beginning

  // Declarations
  int idx, k, i;                 // indices at start and end
  double q, gv, freq, myfloat;   // attenuation, group velocity and frequency
  int n, l, myint;               // n and l of mode
  int sizen, sizeq, sizeu;

  // find the number of eigenfunctions
  int npts, nsidx, ntotidx;
  npts = 185;
  npts = inp_model.LayerUpperIndex(inp_model.NumberOfLayers() - 1) + 1;
  nsidx = 2 * sizeof(n) + (6 * npts + 3) * sizeof(q);   // accounting for data
                                                        //  output
  nsidx = nsidx + 4 * sizeof(i);   // accounting for fortran line data
  ntotidx = file_size / nsidx;
  sizen = sizeof(n);
  sizeq = sizeof(q);

  /* ///////////////////////////////////////////////////////////////////////
      //////Extraction of the data from the binary MINEOS output
      files///////
      ///////////////////////////////////////////////////////////////////////
  */
  // bytes in one eigenfunction
  sizeu = npts * 8;

  //  std::vector<double> uvec(npts);
  std::vector<double> longvec(npts * 6);

  _mode_storage.resize(ntotidx);
  for (int idxl = 0; idxl < ntotidx; idxl++) {
    // separator 4 byte sequence
    ///////////////////////////////////////////////////////////////////////
    infile.read(reinterpret_cast<char *>(&i), sizeof(i));
    ///////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    // get mode information straight into vector
    infile.read(reinterpret_cast<char *>(&n), sizen);
    infile.read(reinterpret_cast<char *>(&l), sizen);
    infile.read(reinterpret_cast<char *>(&freq), sizeq);
    infile.read(reinterpret_cast<char *>(&q), sizeq);
    infile.read(reinterpret_cast<char *>(&gv), sizeq);

    // fill out vectors from arrays
    infile.read(reinterpret_cast<char *>(longvec.data()), sizeu * 6);

    // separator 4 byte sequence
    ///////////////////////////////////////////////////////////////////////
    infile.read(reinterpret_cast<char *>(&i), sizen);
    infile.read(reinterpret_cast<char *>(&i), sizen);
    infile.read(reinterpret_cast<char *>(&i), sizen);
    ///////////////////////////////////////////////////////////////////////

    _mode_storage[idxl] = mineos_eigenfunction_continuous(
        n, l, q, gv, freq / (2.0 * 3.1415926535), longvec, inp_model,
        _vec_radii);
    // }
    _vec_nl.push_back(std::vector<int>{n, l});
  }

  infile.close();
  if (!infile.good()) {
    std::cout << "Error occurred at reading time!" << std::endl;
  }
};

class eigenfunction_toroidal {
public:
  using myvector = std::vector<double>;
  using myiter = myvector::iterator;
  using InterpA = Interpolation::CubicSpline<myiter, myiter>;
  // using InterpA = Interpolation::Akima<myiter, myiter>;

  eigenfunction_toroidal() {};
  eigenfunction_toroidal(int n, int l, double q, double gv, double w,
                         std::vector<double> &u, std::vector<double> &up)
      : _n(n), _l(l), _q(q), _gv(gv), _w(w), _u(u), _up(up) {};

  eigenfunction_toroidal(int n, int l, double q, double gv, double w,
                         std::vector<double> &vec_longinput,
                         EarthModels::ModelInput<double, int> &inp_model,
                         std::vector<std::vector<double>> &vec_radii)
      : _n(n), _l(l), _q(q), _gv(gv), _w(w), _earthmodel(inp_model) {
    int lnum = vec_longinput.size() / 2;
    auto myit = vec_longinput.begin();
    _u = std::vector<double>(myit, myit + lnum);
    _up = std::vector<double>(myit + lnum, myit + 2 * lnum);

    double mfact = 2.0 * 3.14159265358979 * w;

    for (int idx = 0; idx < inp_model.NumberOfLayers(); ++idx) {
      auto itr1 = vec_radii[idx].begin();
      auto itr2 = vec_radii[idx].end();

      auto itu1 = _u.begin() + inp_model.LayerLowerIndex(idx);
      auto itu2 = _u.begin() + inp_model.LayerUpperIndex(idx) + 1;
      _ul.push_back(std::vector<double>(itu1, itu2));
      auto tmpu = InterpA(itr1, itr2, _ul[idx].begin());
      _ucont.push_back(tmpu);

      auto itup1 = _up.begin() + inp_model.LayerLowerIndex(idx);
      auto itup2 = _up.begin() + inp_model.LayerUpperIndex(idx) + 1;
      _upl.push_back(std::vector<double>(itup1, itup2));
      auto tmpup = InterpA(itr1, itr2, _upl[idx].begin());
      _upcont.push_back(tmpup);
    };
  }

  // returns
  auto n() { return _n; };
  auto l() { return _l; };
  auto q() { return _q; };
  auto gv() { return _gv; };
  auto w() { return _w; };

  auto u() { return _u; };
  auto up() { return _up; };
  auto ul() { return _ul; };

  // return continuous u
  auto u(int i) { return _ucont[i]; };
  auto up(int i) { return _upcont[i]; };

private:
  int _n, _l;
  double _q, _gv, _w;
  std::vector<double> _u, _up;                  //, _ku, _kr, _kd;
  std::vector<std::vector<double>> _ul, _upl;   //, _ku, _kr, _kd;
  std::vector<InterpA> _ucont, _upcont;
  EarthModels::ModelInput<double, int> _earthmodel;
};

class modecattoroidal {
public:
  modecattoroidal() {};
  modecattoroidal(std::string &, EarthModels::ModelInput<double, int> &);

  auto allmodes() { return _mode_storage; };
  auto isincatalogue(int n, int l) {
    bool boolout = false;
    for (int idx = 0; idx < _vec_nl.size(); ++idx) {
      //  int idxoutput;

      if (n == _vec_nl[idx][0]) {
        if (l == _vec_nl[idx][1]) {
          boolout = true;
          break;
        }
      }
    }
    return boolout;
  };
  auto singlemode(int i) {
    assert(i < _mode_storage.size() && "Not that many modes!");
    return _mode_storage[i];
  };
  eigenfunction_toroidal &singlemodep(int i) {
    assert(i < _mode_storage.size() && "Not that many modes!");
    return _mode_storage[i];
  };
  auto NumberOfModes() { return _mode_storage.size(); };

  auto AllModes(const EarthMesh::RadialMesh &radmesh) {
    std::vector<std::vector<std::vector<double>>> vec_minval;
    for (int idxe = 0; idxe < radmesh.NE(); ++idxe) {
      std::vector<double> vec_x;
      std::vector<std::vector<double>> vec_out1, vec_out2;

      int laynum = radmesh.LayerNumber(idxe);
      for (int idxq = 0; idxq < radmesh.NN(); ++idxq) {
        std::vector<double> vec_x2;
        for (int idx = 0; idx < this->NumberOfModes(); ++idx) {
          auto mode = this->singlemodep(idx);
          vec_x2.push_back(mode.u(laynum)(radmesh.NodeRadius(idxe, idxq)));
        }
        vec_out2.push_back(vec_x2);
      }
      vec_minval.push_back(vec_out2);
    }
    return vec_minval;
  };

  auto AllModesDeriv(const EarthMesh::RadialMesh &radmesh) {
    std::vector<std::vector<std::vector<double>>> vec_minval;
    for (int idxe = 0; idxe < radmesh.NE(); ++idxe) {
      std::vector<double> vec_x;
      std::vector<std::vector<double>> vec_out1, vec_out2;

      int laynum = radmesh.LayerNumber(idxe);
      for (int idxq = 0; idxq < radmesh.NN(); ++idxq) {
        std::vector<double> vec_x2;
        for (int idx = 0; idx < this->NumberOfModes(); ++idx) {
          auto mode = this->singlemodep(idx);
          vec_x2.push_back(mode.up(laynum)(radmesh.NodeRadius(idxe, idxq)));
        }
        vec_out2.push_back(vec_x2);
      }
      vec_minval.push_back(vec_out2);
    }
    return vec_minval;
  };

private:
  std::vector<eigenfunction_toroidal> _mode_storage;
  std::vector<std::vector<int>> _vec_nl;
  std::vector<std::vector<double>> _vec_radii;
  std::size_t numnodes;
};

modecattoroidal::modecattoroidal(
    std::string &pathtofile, EarthModels::ModelInput<double, int> &inp_model) {
  _vec_radii = inp_model.LayerRadii();
  // open the file for input, in binary form
  std::ifstream infile(pathtofile, std::ifstream::binary);

  numnodes = inp_model.LayerUpperIndex(inp_model.NumberOfLayers() - 1) + 1;
  // check opened correctly
  if (!infile) {
    std::cout << "Cannot open file!" << std::endl;
  }

  // find the length of the file
  infile.seekg(0, std::ios::end);
  std::size_t file_size = infile.tellg();
  std::cout << "Size of the file is " << file_size << " bytes" << std::endl;
  infile.seekg(0, infile.beg);   // reset back to beginning

  // Declarations
  int idx, k, i;                 // indices at start and end
  double q, gv, freq, myfloat;   // attenuation, group velocity and frequency
  int n, l, myint;               // n and l of mode
  int sizen, sizeq, sizeu;

  // find the number of eigenfunctions
  int npts, nsidx, ntotidx;
  npts = numnodes;
  npts = inp_model.LayerUpperIndex(inp_model.NumberOfLayers() - 1) + 1;
  nsidx = 2 * sizeof(n) + (2 * npts + 3) * sizeof(q);   // accounting for data
                                                        //  output
  nsidx = nsidx + 4 * sizeof(i);   // accounting for fortran line data
  ntotidx = file_size / nsidx;
  sizen = sizeof(n);
  sizeq = sizeof(q);

  /* ///////////////////////////////////////////////////////////////////////
      //////Extraction of the data from the binary MINEOS output
      files///////
      ///////////////////////////////////////////////////////////////////////
  */
  // bytes in one eigenfunction
  sizeu = npts * 8;

  //  std::vector<double> uvec(npts);
  std::vector<double> longvec(npts * 2);

  _mode_storage.resize(ntotidx);
  for (int idxl = 0; idxl < ntotidx; idxl++) {
    // separator 4 byte sequence
    ///////////////////////////////////////////////////////////////////////
    infile.read(reinterpret_cast<char *>(&i), sizeof(i));
    ///////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    // get mode information straight into vector
    infile.read(reinterpret_cast<char *>(&n), sizen);
    infile.read(reinterpret_cast<char *>(&l), sizen);
    infile.read(reinterpret_cast<char *>(&freq), sizeq);
    infile.read(reinterpret_cast<char *>(&q), sizeq);
    infile.read(reinterpret_cast<char *>(&gv), sizeq);

    // fill out vectors from arrays
    infile.read(reinterpret_cast<char *>(longvec.data()), sizeu * 2);

    // separator 4 byte sequence
    ///////////////////////////////////////////////////////////////////////
    infile.read(reinterpret_cast<char *>(&i), sizen);
    infile.read(reinterpret_cast<char *>(&i), sizen);
    infile.read(reinterpret_cast<char *>(&i), sizen);
    ///////////////////////////////////////////////////////////////////////

    _mode_storage[idxl] =
        eigenfunction_toroidal(n, l, q, gv, freq / (2.0 * 3.1415926535),
                               longvec, inp_model, _vec_radii);
    // }
    _vec_nl.push_back(std::vector<int>{n, l});
  }

  infile.close();
  if (!infile.good()) {
    std::cout << "Error occurred at reading time!" << std::endl;
  }
};

class modecatalogue {

public:
  modecatalogue() {};
  modecatalogue(std::string &, EarthModels::ModelInput<double, int> &);

  auto allmodes() { return _mode_storage; };
  auto isincatalogue(int n, int l) {
    bool boolout = false;
    for (int idx = 0; idx < _vec_nl.size(); ++idx) {
      //  int idxoutput;

      if (n == _vec_nl[idx][0]) {
        if (l == _vec_nl[idx][1]) {
          boolout = true;
          break;
        }
      }
    }
    return boolout;
  };
  auto singlemodeoutput(int i) {
    assert(i < _mode_storage.size() && "Not that many modes!");
    return _mode_storage[i];
  };
  auto NumberOfModes() { return _mode_storage.size(); };

private:
  std::vector<mineos_eigenfunction> _mode_storage;
  std::vector<std::vector<int>> _vec_nl;
  std::vector<std::vector<double>> _vec_radii;
};

modecatalogue::modecatalogue(std::string &pathtofile,
                             EarthModels::ModelInput<double, int> &inp_model) {

  // radii vector
  _vec_radii = inp_model.LayerRadii();
  // open the file for input, in binary form
  std::ifstream infile(pathtofile, std::ifstream::binary);

  // check opened correctly
  if (!infile) {
    std::cout << "Cannot open file!" << std::endl;
  }

  // find the length of the file
  infile.seekg(0, std::ios::end);
  std::size_t file_size = infile.tellg();
  std::cout << "Size of the file is " << file_size << " bytes" << std::endl;
  infile.seekg(0, infile.beg);   // reset back to beginning

  // Declarations
  int idx, k, i;                 // indices at start and end
  double q, gv, freq, myfloat;   // attenuation, group velocity and frequency
  int n, l, myint;               // n and l of mode
  int sizen, sizeq, sizeu;

  // find the number of eigenfunctions
  int npts, nsidx, ntotidx;
  npts = 185;
  nsidx = 2 * sizeof(n) + (6 * npts + 3) * sizeof(q);   // accounting for data
                                                        //  output
  nsidx = nsidx + 4 * sizeof(i);   // accounting for fortran line data
  ntotidx = file_size / nsidx;
  sizen = sizeof(n);
  sizeq = sizeof(q);

  /* ///////////////////////////////////////////////////////////////////////
      //////Extraction of the data from the binary MINEOS output
      files///////
      ///////////////////////////////////////////////////////////////////////
  */
  // bytes in one eigenfunction
  sizeu = npts * 8;

  //  std::vector<double> uvec(npts);
  std::vector<double> longvec(npts * 6);

  _mode_storage.resize(ntotidx);
  for (int idxl = 0; idxl < ntotidx; idxl++) {
    // separator 4 byte sequence
    ///////////////////////////////////////////////////////////////////////
    infile.read(reinterpret_cast<char *>(&i), sizeof(i));
    ///////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    // get mode information straight into vector
    infile.read(reinterpret_cast<char *>(&n), sizen);
    infile.read(reinterpret_cast<char *>(&l), sizen);
    infile.read(reinterpret_cast<char *>(&freq), sizeq);
    infile.read(reinterpret_cast<char *>(&q), sizeq);
    infile.read(reinterpret_cast<char *>(&gv), sizeq);

    // fill out vectors from arrays
    infile.read(reinterpret_cast<char *>(longvec.data()), sizeu * 6);

    // separator 4 byte sequence
    ///////////////////////////////////////////////////////////////////////
    infile.read(reinterpret_cast<char *>(&i), sizen);
    infile.read(reinterpret_cast<char *>(&i), sizen);
    infile.read(reinterpret_cast<char *>(&i), sizen);
    ///////////////////////////////////////////////////////////////////////

    _mode_storage[idxl] =
        mineos_eigenfunction(n, l, q, gv, freq / (2.0 * 3.1415926535), longvec);
    if (idxl == 0) {
      auto tmpeig = mineos_eigenfunction_continuous(
          n, l, q, gv, freq / (2.0 * 3.1415926535), longvec, inp_model,
          _vec_radii);
    }
    _vec_nl.push_back(std::vector<int>{n, l});
  }

  infile.close();
  if (!infile.good()) {
    std::cout << "Error occurred at reading time!" << std::endl;
  }
};

class modesummary {

public:
  modesummary(ModeCoupling::modecataloguecontinuous &,
              EarthModels::ModelInput<double, int> &, int);

  // returns
  std::vector<std::vector<std::vector<std::vector<double>>>> &uvec() {
    return vec_uval;
  };
  double uvec(int idx1, int idx2, int idx3, int idx4) {
    return vec_uval[idx1][idx2][idx3][idx4];
  };
  std::vector<std::vector<std::vector<std::vector<double>>>> &upvec() {
    return vec_upval;
  };
  double upvec(int idx1, int idx2, int idx3, int idx4) {
    return vec_upval[idx1][idx2][idx3][idx4];
  };
  std::vector<std::vector<std::vector<std::vector<double>>>> &vvec() {
    return vec_vval;
  };
  double vvec(int idx1, int idx2, int idx3, int idx4) {
    return vec_vval[idx1][idx2][idx3][idx4];
  };
  std::vector<std::vector<std::vector<std::vector<double>>>> &vpvec() {
    return vec_vpval;
  };
  double vpvec(int idx1, int idx2, int idx3, int idx4) {
    return vec_vpval[idx1][idx2][idx3][idx4];
  };
  std::vector<std::vector<std::vector<double>>> &multvec() { return vec_wdrr; };
  double multvec(int idx1, int idx2, int idx3) {
    return vec_wdrr[idx1][idx2][idx3];
  };
  int npoints() { return _q.N(); };
  double w(int idx) { return vec_w[idx]; };
  int NumberOfModes() { return vec_w.size(); };
  int NumberOfLayers() { return numlayers; };
  auto q() { return _q; };

private:
  std::vector<std::vector<std::vector<std::vector<double>>>> vec_uval, vec_vval,
      vec_upval, vec_vpval;
  std::vector<std::vector<std::vector<double>>> vec_wdrr;
  std::vector<double> vec_w;
  GaussQuad::Quadrature1D<double> _q;
  int numlayers;
};

modesummary::modesummary(ModeCoupling::modecataloguecontinuous &mode_cat,
                         EarthModels::ModelInput<double, int> &inp_model,
                         int npoints)
    : numlayers{inp_model.NumberOfLayers()} {
  _q = GaussQuad::GaussLobattoLegendreQuadrature1D<double>(npoints);
  std::vector<double> qx = _q.Points();
  std::vector<double> qw = _q.Weights();
  int qsize = qx.size();

  int maxnum = mode_cat.NumberOfModes();

  double densitynorm = 5515.0;
  double radiusnorm = inp_model.OuterRadius() * inp_model.LengthNorm();
  //    std::cout << radiusnorm << "\n\n";
  double pi_db = 3.14159265358979;
  double bigg_db = 6.6723 * std::pow(10.0, -11.0);
  double velocitynorm = radiusnorm / std::sqrt(pi_db * bigg_db * densitynorm);
  double frequencynorm = velocitynorm / radiusnorm;

  for (int idxlayer = 0; idxlayer < inp_model.NumberOfLayers(); ++idxlayer) {
    auto vec_radii = inp_model.LayerRadii(idxlayer);
    std::vector<std::vector<double>> vec_tmp1;
    // std::cout << idxlayer << "\n";
    for (int idxpl = 0; idxpl < vec_radii.size() - 1; ++idxpl) {

      double x1 = vec_radii[idxpl];
      double x2 = vec_radii[idxpl + 1];
      double delta = (x2 - x1) / 2.0;
      double xpl = (x1 + x2) / 2.0;
      double tmpsum = 0.0;
      std::vector<double> vec_tmp2(qsize, 0.0);
      for (int idxq = 0; idxq < qsize; ++idxq) {
        double xscale = delta * qx[idxq] + xpl;
        vec_tmp2[idxq] =
            qw[idxq] * inp_model.Density(idxlayer)(xscale) * xscale * xscale;
      }
      vec_tmp1.push_back(vec_tmp2);
    };
    vec_wdrr.push_back(vec_tmp1);
  };
  for (int idx1 = 0; idx1 < maxnum; ++idx1) {
    std::vector<std::vector<std::vector<double>>> vec_tmp1, vec_tmp11,
        vec_tmpup1, vec_tmpvp1;
    vec_w.push_back(mode_cat.singlemodep(idx1).w());
    for (int idxlayer = 0; idxlayer < inp_model.NumberOfLayers(); ++idxlayer) {
      auto vec_radii = inp_model.LayerRadii(idxlayer);
      std::vector<std::vector<double>> vec_tmp2, vec_tmp21, vec_tmpup2,
          vec_tmpvp2;
      for (int idxpl = 0; idxpl < vec_radii.size() - 1; ++idxpl) {

        double x1 = vec_radii[idxpl];
        double x2 = vec_radii[idxpl + 1];
        double delta = (x2 - x1) / 2.0;
        double xpl = (x1 + x2) / 2.0;
        // double tmpsum = 0.0;
        std::vector<double> vec_tmp3(qsize, 0.0), vec_tmp31(qsize, 0.0),
            vec_tmpup3(qsize, 0.0), vec_tmpvp3(qsize, 0.0);
        for (int idxq = 0; idxq < qsize; ++idxq) {
          double xscale = delta * qx[idxq] + xpl;
          vec_tmp3[idxq] = mode_cat.singlemodep(idx1).u(idxlayer)(xscale);
          vec_tmp31[idxq] = mode_cat.singlemodep(idx1).v(idxlayer)(xscale);
          vec_tmpup3[idxq] = mode_cat.singlemodep(idx1).up(idxlayer)(xscale);
          vec_tmpvp3[idxq] = mode_cat.singlemodep(idx1).vp(idxlayer)(xscale);
        }
        vec_tmp2.push_back(vec_tmp3);
        vec_tmp21.push_back(vec_tmp31);
        vec_tmpup2.push_back(vec_tmpup3);
        vec_tmpvp2.push_back(vec_tmpvp3);
      };
      vec_tmp1.push_back(vec_tmp2);
      vec_tmp11.push_back(vec_tmp21);
      vec_tmpup1.push_back(vec_tmpup2);
      vec_tmpvp1.push_back(vec_tmpvp2);
    };
    vec_uval.push_back(vec_tmp1);
    vec_vval.push_back(vec_tmp11);
    vec_upval.push_back(vec_tmpup1);
    vec_vpval.push_back(vec_tmpvp1);
  }
};

class modesummary2 {

public:
  modesummary2(ModeCoupling::modecataloguecontinuous &,
               EarthModels::ModelInput<double, int> &, int);

  // returns
  std::vector<std::vector<std::vector<double>>> &uvec() { return vec_uval2; };
  double uvec(int idx1, int idx2, int idx3) {
    return vec_uval2[idx1][idx2][idx3];
  };
  std::vector<std::vector<std::vector<double>>> &upvec() { return vec_upval2; };
  double upvec(int idx1, int idx2, int idx3) {
    return vec_upval2[idx1][idx2][idx3];
  };
  std::vector<std::vector<std::vector<double>>> &vvec() { return vec_vval2; };
  double vvec(int idx1, int idx2, int idx3) {
    return vec_vval2[idx1][idx2][idx3];
  };
  std::vector<std::vector<std::vector<double>>> &vpvec() { return vec_vpval2; };
  double vpvec(int idx1, int idx2, int idx3) {
    return vec_vpval2[idx1][idx2][idx3];
  };
  std::vector<std::vector<double>> &multvec() { return vec_wdrr2; };
  double multvec(int idx1, int idx2) { return vec_wdrr2[idx1][idx2]; };
  int npoints() { return _q.N(); };
  double w(int idx) { return vec_w[idx]; };
  int NumberOfModes() { return vec_w.size(); };
  int NumberOfLayers() { return numlayers; };
  auto q() { return _q; };
  EarthMesh::RadialMesh &fullmesh() { return _radial_mesh; };

private:
  std::vector<std::vector<std::vector<std::vector<double>>>> vec_uval, vec_vval,
      vec_upval, vec_vpval;
  std::vector<std::vector<std::vector<double>>> vec_wdrr;

  std::vector<std::vector<std::vector<double>>> vec_uval2, vec_vval2,
      vec_upval2, vec_vpval2;
  std::vector<std::vector<double>> vec_wdrr2;
  std::vector<double> vec_w;
  GaussQuad::Quadrature1D<double> _q;
  int numlayers;
  EarthMesh::RadialMesh _radial_mesh;
};

modesummary2::modesummary2(ModeCoupling::modecataloguecontinuous &mode_cat,
                           EarthModels::ModelInput<double, int> &inp_model,
                           int npoints)
    : numlayers{inp_model.NumberOfLayers()} {
  _q = GaussQuad::GaussLobattoLegendreQuadrature1D<double>(npoints);
  std::vector<double> qx = _q.Points();
  std::vector<double> qw = _q.Weights();
  int qsize = qx.size();

  // declare radial mesh
  _radial_mesh = EarthMesh::RadialMesh(inp_model, qsize, 0.01, 1.0, false);

  int maxnum = mode_cat.NumberOfModes();

  double densitynorm = 5515.0;
  double radiusnorm = inp_model.OuterRadius() * inp_model.LengthNorm();
  //    std::cout << radiusnorm << "\n\n";
  double pi_db = 3.14159265358979;
  double bigg_db = 6.6723 * std::pow(10.0, -11.0);
  double velocitynorm = radiusnorm / std::sqrt(pi_db * bigg_db * densitynorm);
  double frequencynorm = velocitynorm / radiusnorm;

  // evaluate on the mesh
  for (int idxe = 0; idxe < _radial_mesh.NE(); ++idxe) {
    // auto vec_radii = inp_model.LayerRadii(idxlayer);
    // std::vector<std::vector<double>> vec_tmp1;
    // std::cout << idxlayer << "\n";
    int laynum = _radial_mesh.LayerNumber(idxe);
    std::vector<double> vec_tmp2(qsize, 0.0);
    for (int idxq = 0; idxq < qsize; ++idxq) {
      double xscale = _radial_mesh.NodeRadius(idxe, idxq);
      vec_tmp2[idxq] =
          qw[idxq] * inp_model.Density(laynum)(xscale) * xscale * xscale;
    }
    vec_wdrr2.push_back(vec_tmp2);
  };

  for (int idx1 = 0; idx1 < maxnum; ++idx1) {
    std::vector<std::vector<double>> vec_tmp2, vec_tmp21, vec_tmpup2,
        vec_tmpvp2;
    vec_w.push_back(mode_cat.singlemodep(idx1).w());
    for (int idxe = 0; idxe < _radial_mesh.NE(); ++idxe) {
      // auto vec_radii = inp_model.LayerRadii(idxlayer);
      // std::vector<std::vector<double>> vec_tmp1;
      // std::cout << idxlayer << "\n";
      std::vector<double> vec_tmp3(qsize, 0.0), vec_tmp31(qsize, 0.0),
          vec_tmpup3(qsize, 0.0), vec_tmpvp3(qsize, 0.0);

      int idxlayer = _radial_mesh.LayerNumber(idxe);
      for (int idxq = 0; idxq < qsize; ++idxq) {
        double xscale = _radial_mesh.NodeRadius(idxe, idxq);

        vec_tmp3[idxq] = mode_cat.singlemodep(idx1).u(idxlayer)(xscale);
        vec_tmp31[idxq] = mode_cat.singlemodep(idx1).v(idxlayer)(xscale);
        vec_tmpup3[idxq] = mode_cat.singlemodep(idx1).up(idxlayer)(xscale);
        vec_tmpvp3[idxq] = mode_cat.singlemodep(idx1).vp(idxlayer)(xscale);
      }
      vec_tmp2.push_back(vec_tmp3);
      vec_tmp21.push_back(vec_tmp31);
      vec_tmpup2.push_back(vec_tmpup3);
      vec_tmpvp2.push_back(vec_tmpvp3);
    };
    vec_uval2.push_back(vec_tmp2);
    vec_vval2.push_back(vec_tmp21);
    vec_upval2.push_back(vec_tmpup2);
    vec_vpval2.push_back(vec_tmpvp2);
  }
};

}   // namespace ModeCoupling

#endif