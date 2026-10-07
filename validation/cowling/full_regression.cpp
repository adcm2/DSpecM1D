#include <DSpecM1D/src/FullSpec.h>
#include <DSpecM1D/src/InputParametersNew.h>
#include "test_utils.h"
#include <fstream>
#include <string>
void dump(std::ofstream &out, const Eigen::MatrixXcd &m) {
  const long rows = m.rows(), cols = m.cols();
  out.write(reinterpret_cast<const char*>(&rows), sizeof(rows));
  out.write(reinterpret_cast<const char*>(&cols), sizeof(cols));
  for (long j = 0; j < cols; ++j) for (long i = 0; i < rows; ++i) {
    double r = m(i,j).real(), v = m(i,j).imag();
    out.write(reinterpret_cast<const char*>(&r), sizeof(r));
    out.write(reinterpret_cast<const char*>(&v), sizeof(v));
  }
}
int main(int argc, char **argv) {
  if (argc != 2) return 2;
  std::ofstream out(argv[1], std::ios::binary);
  for (bool atten : {false, true}) {
    DSpecMTest::TempDir tmp;
    DSpecMTest::ParameterOptions o;
    o.type = 3; o.attenuation = atten; o.lmin = 1; o.lmax = 4;
    o.f1 = o.f11 = o.f12 = 2.0; o.f2 = o.f21 = o.f22 = 8.0;
    o.tOutMinutes = 1.0; o.timeStepSec = 5.0;
    o.relativeError = 1.0;
    auto model = DSpecMTest::repoRoot()/"data"/"models"/(atten?"prem.200.no.txt":"prem.200.noatten.txt");
    auto p = DSpecMTest::writeFile(tmp.path()/"params.txt", DSpecMTest::makeParameterText(model.string(),o));
    InputParametersNew inp(p.string(),3,2,0.05,1.0,0.05,0.0,1);
    Full1D::SEM sem(inp);
    for (int l : {1,2,7}) {
      dump(out, Eigen::MatrixXcd(sem.hS(l).cast<std::complex<double>>()));
      dump(out, Eigen::MatrixXcd(sem.hSa(l).cast<std::complex<double>>()));
      dump(out, Eigen::MatrixXcd(sem.pS(l).cast<std::complex<double>>()));
      dump(out, sem.calculateForceAll(inp.cmt(),l));
      dump(out, sem.rvBaseFull(inp.inputParameters(),l));
    }
    SPARSESPEC::SparseFSpec solver;
    dump(out, solver.spectra(inp,sem));
    dump(out, solver.spectra(inp));
  }
  return out.good()?0:3;
}
