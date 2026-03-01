#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <chrono>
#include <iomanip>

using std::cout;
using std::endl;
using std::ifstream;
using std::ios;
using std::setprecision;
using std::string;
using std::vector;

// using namespace std;
using namespace std::chrono;

struct mineos_eigenfunction {
  int n, l;
  double q, group_velocity, w;
  vector<double> u = vector<double>(185);
  vector<double> up = vector<double>(185);
  vector<double> v = vector<double>(185);
  vector<double> vp = vector<double>(185);
  vector<double> p = vector<double>(185);
  vector<double> pp = vector<double>(185);
  vector<double> ku = vector<double>(185);
  vector<double> kr = vector<double>(185);
  vector<double> kd = vector<double>(185);

  // vector<double> u, up, v, vp, p, pp;
};

struct eigenfunction_toroidal {
  int n, l;
  double q, group_velocity, w;

  vector<double> u;
  vector<double> up;
  // vector<double> ku = vector<double>(185);
  // vector<double> kr = vector<double>(185);
  // vector<double> kd = vector<double>(185);

  eigenfunction_toroidal(int n) {
    u = vector<double>(n);
    up = vector<double>(n);
  };
  // vector<double> u, up, v, vp, p, pp;
};

struct eigfunc_complete {
  std::vector<mineos_eigenfunction> normal_modes;
};

int
main() {
  // open the file for input, in binary form
  ifstream infile("/space/adcm2/mineos-1.0.2/OUTPUT/eprem_noocean_S_IC",
                  ifstream::binary);

  // check opened correctly
  if (!infile) {
    cout << "Cannot open file!" << endl;
    return 1;
  }

  // find the length of the file
  infile.seekg(0, ios::end);
  int file_size = infile.tellg();
  std::cout << "Size of the file is " << file_size << " bytes" << endl;
  infile.seekg(0, infile.beg);   // reset back to beginning

  // get number of lines

  // Declarations
  int idx, k, i;   // indices at start and end
  // long int i;
  bool mydo, mydo2;              // while loop booleans
  double q, gv, freq, myfloat;   // attenuation, group velocity and frequency
  int n, l, myint;               // n and l of mode
  int sizen, sizeq, sizeu, sizearray;

  // find the number of eigenfunctions
  int npts, nsidx, ntotidx;
  npts = 218;
  nsidx = 2 * sizeof(n) +
          (2 * npts + 3) * sizeof(q);   // accounting for data output
  nsidx = nsidx + 4 * sizeof(i);        // accounting for fortran line data
  ntotidx = file_size / nsidx;
  sizen = sizeof(n);
  sizeq = sizeof(q);
  // std::cout << file_size << std::endl;
  // std::cout << nsidx << std::endl;
  // std::cout << "Number of sets is " << ntotidx << endl;

  /*  ///////////////////////////////////////////////////////////////////////
      //////Extraction of the data from the binary MINEOS output files///////
      ///////////////////////////////////////////////////////////////////////
  */

  // float freq;
  double u[npts], up[npts], v[npts], vp[npts], p[npts],
      pp[npts];   // arrays to store extracted data
  double readarray[6 * npts];
  sizeu = sizeof(u);
  sizearray = sizeof(readarray);

  // timing
  auto start = high_resolution_clock::now();

  // setting booleans
  mydo2 = 1;
  idx = 0;

  // use class
  mineos_eigenfunction myeig;
  vector<eigenfunction_toroidal> myeigvec(ntotidx,
                                          eigenfunction_toroidal(npts));
  // eigfunc_complete myeigvec;
  // myeigvec.normal_modes(ntotidx);

  // looping and extracting all data
  //  std::cout << setprecision(5);

  for (idx == 0; idx < ntotidx; idx++) {
    // separator 4 byte sequence
    ///////////////////////////////////////////////////////////////////////
    infile.read(reinterpret_cast<char *>(&i), sizeof(i));
    ///////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    // get mode information straight into vector
    infile.read(reinterpret_cast<char *>(&myeigvec[idx].n), sizen);
    infile.read(reinterpret_cast<char *>(&myeigvec[idx].l), sizen);
    infile.read(reinterpret_cast<char *>(&myeigvec[idx].w), sizeq);
    infile.read(reinterpret_cast<char *>(&myeigvec[idx].q), sizeq);
    infile.read(reinterpret_cast<char *>(&myeigvec[idx].group_velocity), sizeq);
    myeigvec[idx].w = myeigvec[idx].w / (2.0 * 3.1415926535);

    // fill out vectors from arrays
    infile.read(reinterpret_cast<char *>(myeigvec[idx].u.data()), sizeu);
    infile.read(reinterpret_cast<char *>(myeigvec[idx].up.data()), sizeu);
    // infile.read(reinterpret_cast<char *>(myeigvec[idx].v.data()), sizeu);
    // infile.read(reinterpret_cast<char *>(myeigvec[idx].vp.data()), sizeu);
    // infile.read(reinterpret_cast<char *>(myeigvec[idx].p.data()), sizeu);
    // infile.read(reinterpret_cast<char *>(myeigvec[idx].pp.data()), sizeu);

    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////

    // separator 4 byte sequence
    ///////////////////////////////////////////////////////////////////////
    infile.read(reinterpret_cast<char *>(&i), sizen);
    infile.read(reinterpret_cast<char *>(&i), sizen);
    infile.read(reinterpret_cast<char *>(&i), sizen);
    ///////////////////////////////////////////////////////////////////////
  }

  // close and check closed
  infile.close();
  // std::cout << "Hello 3" << std::endl;
  if (!infile.good()) {
    cout << "Error occurred at reading time!" << endl;
    return 1;
  }
  // std::cout << "Hello 4" << std::endl;
  // assigning to appropriate class structure:
  // eigfunc_complete myeigvec2;
  // myeigvec2.normal_modes = myeigvec;
  // k = 500;
  // std::cout << "The size of the vector is " << myeigvec.size() << std::endl;
  idx = 0;
  std::cout << setprecision(12) << myeigvec[idx].n << endl;
  std::cout << setprecision(12) << myeigvec[idx].l << endl;
  std::cout << setprecision(12) << myeigvec[idx].q << endl;
  std::cout << setprecision(12) << myeigvec[idx].group_velocity << endl;
  std::cout << setprecision(12) << myeigvec[idx].w * 1000.0 << endl;
  std::cout << myeigvec[idx].u[0] << " " << myeigvec[idx].u[1] << std::endl;
  std::cout << "\n\n";
  for (auto idx : myeigvec[idx].u) {
    std::cout << idx << "\n";
  }
  std::cout << "\n\n";
  std::cout << std::endl;

  idx = 1;
  std::cout << setprecision(12) << myeigvec[idx].n << endl;
  std::cout << setprecision(12) << myeigvec[idx].l << endl;
  std::cout << setprecision(12) << myeigvec[idx].q << endl;
  std::cout << setprecision(12) << myeigvec[idx].group_velocity << endl;
  std::cout << setprecision(12) << myeigvec[idx].w * 1000.0 << endl;
  // std::cout << "\n\n";
  // for (int idx : myeigvec[idx].u) {
  //   std::cout << idx << "\n";
  // }
  // std::cout << "\n\n";
  std::cout << std::endl;

  // find the sensitivity kernels
  // for (idx = 0; idx < ntotidx; ++idx)
  // {
  //     //Ku which doesn't need integration
  //     myeigvec[idx].ku = (1/3.0) * (2.0 * )
  // }
  // finish timing
  auto stop = high_resolution_clock::now();

  // difference
  auto duration = duration_cast<microseconds>(stop - start);
  std::cout << "Time: " << duration.count() / 1000000.0 << "s" << endl;
  // std::cout << "Hello 5" << std::endl;
}