#include <DSpecM1D/SensitivityKernels>
#include <DSpecM1D/ModelInput>
#include <DSpecM1D/FrequencyTools>
#include <DSpecM1D/src/InputParser.h>
#include <DSpecM1D/src/NormClass.h>
#include <DSpecM1D/src/SourceInfo.h>

#include <filesystem>
#include <algorithm>
#include <complex>
#include <fstream>
#include <iomanip>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

namespace fs = std::filesystem;

constexpr int kLmin=0, kLmax=2, kNQ=4;
constexpr double kMaxstep=0.25;
constexpr double kF1=0.1, kF2=0.8, kF11=0.1, kF12=0.2, kF21=0.6, kF22=0.8;
constexpr double kTimeStepSeconds=20.0, kOutputMinutes=0.5;
constexpr double kSourceDepthKm=20, kSourceLatitude=33, kSourceLongitude=10;
constexpr double kTensor[6]={1e20,2e20,3e20,4e20,5e20,6e20};
constexpr int kFirstBin=2, kEndBin=4;
const std::complex<double> kQBin2{0.7,0.2}, kQBin3{-0.25,0.45};

static void requireEmpty(const fs::path &out) {
  if (fs::exists(out)) {
    if (!fs::is_directory(out) || fs::directory_iterator(out) != fs::directory_iterator{})
      throw std::runtime_error("output directory must be empty: " + out.string());
  } else fs::create_directories(out);
}

static std::string parameters(const fs::path &model) {
  std::ostringstream s;
  s << "\"unused.out\"\n\"" << fs::absolute(model).string() << "\"\n"
    << "4\n0\n2\n0\n0\n1e-5\n" << kLmin << '\n' << kLmax << '\n'
    << kF1 << '\n' << kF2 << '\n' << kOutputMinutes << '\n' << kTimeStepSeconds << '\n'
    << kF11 << '\n' << kF12 << '\n' << kF21 << '\n' << kF22 << '\n'
    << kSourceDepthKm << '\n' << kSourceLatitude << '\n' << kSourceLongitude << '\n'
    << kTensor[0] << '\n' << kTensor[1] << '\n' << kTensor[2] << '\n'
    << kTensor[3] << '\n' << kTensor[4] << '\n' << kTensor[5] << '\n'
    << "3\n1\n30 65\n";
  return s.str();
}

template <class Base, class Norm>
class PerturbedModel : public Base {
public:
  PerturbedModel(const std::string &path, const Norm &norm,
                 int layer, double lo, double hi, double amplitude,
                 double propertyAmplitude, bool density)
      : Base(path, norm), layer_(layer), lo_(lo), hi_(hi), amp_(amplitude),
        propertyAmplitude_(propertyAmplitude), density_(density) {}
  auto Density(int layer) const {
    auto f = Base::Density(layer);
    const bool active = density_ && layer == layer_;
    return [f, lo=lo_, hi=hi_, a=active?amp_:0.0, strength=propertyAmplitude_, scale=Base::DensityNorm()](double r) {
      const double x = (2*r-lo-hi)/(hi-lo);
      const double b = (r >= lo && r <= hi) ? std::max(0.0, 1-x*x) : 0.0;
      return f(r) + a * strength / scale * b;
    };
  }
  auto A(int layer) const { return modulus(Base::A(layer), layer); }
  auto C(int layer) const { return modulus(Base::C(layer), layer); }
  auto F(int layer) const { return modulus(Base::F(layer), layer); }
  auto L(int layer) const { return Base::L(layer); }
  auto N(int layer) const { return Base::N(layer); }
private:
  template<class F> auto modulus(F f, int layer) const {
    const bool active = !density_ && layer == layer_ && Base::IsFluid(layer);
    return [f, lo=lo_, hi=hi_, a=active?amp_:0.0, strength=propertyAmplitude_, scale=Base::StressNorm()](double r) {
      const double x = (2*r-lo-hi)/(hi-lo);
      const double b = (r >= lo && r <= hi) ? std::max(0.0, 1-x*x) : 0.0;
      return f(r) + a * strength / scale * b;
    };
  }
  int layer_; double lo_, hi_, amp_, propertyAmplitude_; bool density_;
};

int main(int argc, char **argv) try {
  if (argc != 3) throw std::runtime_error("usage: sensitivity_kernels <model.txt> <empty-output-dir>");
  const fs::path modelPath=argv[1], out=argv[2];
  if (!fs::is_regular_file(modelPath)) throw std::runtime_error("model file not found");
  requireEmpty(out); // guard before creating any output file
  const auto paramPath=out/".params";
  { std::ofstream f(paramPath); f << parameters(modelPath); }
  InputParameters p(paramPath.string());
  EarthModels::EarthConstants<double> norm;
  using Model=EarthModels::ModelInput<double>;
  Model base(modelPath.string(), norm);
  SourceInfo::EarthquakeCMT source(p);
  Full1D::SEM sem(base, kMaxstep, kNQ, p.lmax(), source);
  SpectraSolver::FreqFull freq(kF1,kF2,kF11,kF12,kF21,kF22,kTimeStepSeconds,kOutputMinutes,0.1,0.5,0.0,0.3,2,base.TimeNorm());
  Eigen::VectorXd window(freq.nt());
  for (int i=0;i<freq.nt();++i) window[i]=0.5*(1-std::cos(2*PI*i/(freq.nt()-1)));
  Eigen::MatrixXcd q=Eigen::MatrixXcd::Zero(3,freq.nt()/2+1);
  q(0,kFirstBin)=kQBin2; q(0,kFirstBin+1)=kQBin3;
  const auto result=SPARSESPEC::WindowedSensitivityEngine::evaluate(sem,source,p,freq,window,q,kFirstBin,kEndBin);

  int solid=-1, fluid=-1;
  double solidDistance=std::numeric_limits<double>::infinity();
  double fluidDistance=std::numeric_limits<double>::infinity();
  for(int e=0;e<sem.mesh().NE();++e) {
    const double mid=0.5*(sem.mesh().NodeRadius(e,0)+sem.mesh().NodeRadius(e,sem.mesh().NN()-1));
    if(!sem.elementIsFluid(e) && e>=sem.el() && e<sem.eu() && std::abs(mid-0.7*sem.mesh().OR())<solidDistance) {
      solid=e; solidDistance=std::abs(mid-0.7*sem.mesh().OR());
    }
    if(sem.elementIsFluid(e) && std::abs(mid-0.3*sem.mesh().OR())<fluidDistance) {
      fluid=e; fluidDistance=std::abs(mid-0.3*sem.mesh().OR());
    }
  }
  if(solid<0 || fluid<0) throw std::runtime_error("example requires at least one solid and one fluid element");
  const auto &mesh=sem.mesh();
  auto direction=[&](const std::string &name, bool density, int e, double strength, bool selected) {
    const int layer=mesh.LayerNumber(e);
    const double lo=mesh.NodeRadius(e,0), hi=mesh.NodeRadius(e,mesh.NN()-1);
    double predicted=0;
    for(int i=0;i<mesh.NN();++i) {
      const double x=(2*mesh.NodeRadius(e,i)-lo-hi)/(hi-lo), bump=std::max(0.0,1-x*x);
      const double delta=strength*bump;
      const auto &k=density?result.materialKernels.densityKgPerM3:result.materialKernels.bulkModulusPa;
      if(density || sem.elementIsFluid(e)) predicted += k(e,i)*result.materialKernels.drWeightsMetres(e,i)*delta;
    }
    using PM=PerturbedModel<Model, EarthModels::EarthConstants<double>>;
    auto objective=[&](double a) {
      PM altered(modelPath.string(),norm,layer,lo,hi,a,strength,density);
      Full1D::SEM changed(altered,kMaxstep,kNQ,p.lmax(),source);
      return SPARSESPEC::WindowedSensitivityEngine::evaluate(changed,source,p,freq,window,q,kFirstBin,kEndBin).objective;
    };
    std::array<double,3> steps{0.02,0.01,0.005};
    std::array<double,3> measured{};
    for(int j=0;j<3;++j) measured[j]=(objective(steps[j])-objective(-steps[j]))/(2*steps[j]);
    struct DirectionResult { std::string name; double predicted; std::array<double,3> measured; int layer; double strength; bool selected; };
    return DirectionResult{name,predicted,measured,layer,strength,selected};
  };
  const auto d=direction("solid_density",true,solid,100.0,true);
  const auto b10=direction("fluid_bulk_10GPa",false,fluid,1e10,false);
  const auto b100=direction("fluid_bulk_100GPa",false,fluid,1e11,true);
  const auto receiver=p.receivers().at(0);
  const double GUnitSI=base.GravitationalConstant();
  const double GHat=sem.gravityCoefficient();
  const double physicalG=GHat*GUnitSI;

  std::ofstream samples(out/"samples.csv");
  samples << "radius_m,solid,density_kg_m3,A_Pa,C_Pa,F_Pa,L_Pa,N_Pa,fluid_bulk_pa,density_kernel,bulk_kernel,C_kernel,dr_weight_m,layer\n" << std::setprecision(17);
  for(int e=0;e<mesh.NE();++e) for(int i=0;i<mesh.NN();++i) {
    const int layer=mesh.LayerNumber(e); const double r=mesh.NodeRadius(e,i);
    const auto rho=base.Density(layer);
    const auto kap=base.Kappa(layer);
    const auto c=base.C(layer);
    const auto a=base.A(layer);
    const auto f=base.F(layer);
    const auto l=base.L(layer);
    const auto n=base.N(layer);
    samples << mesh.NodeRadius(e,i)*sem.lengthScaleMetres() << ',' << (!sem.elementIsFluid(e)) << ','
      << rho(r)*base.DensityNorm() << ',' << a(r)*base.StressNorm() << ',' << c(r)*base.StressNorm() << ','
      << f(r)*base.StressNorm() << ',' << l(r)*base.StressNorm() << ',' << n(r)*base.StressNorm() << ','
      << (sem.elementIsFluid(e)?kap(r)*base.StressNorm():0.0) << ','
      << result.materialKernels.densityKgPerM3(e,i) << ',' << result.materialKernels.bulkModulusPa(e,i) << ','
      << result.kernels.C(e,i) << ','
      << result.materialKernels.drWeightsMetres(e,i) << ',' << layer << '\n';
  }
  std::ofstream dirs(out/"directions.csv");
  dirs << "direction,element,layer,r_lower_m,r_upper_m,delta_amplitude,selected_for_gate,predicted,fd_0.02,fd_0.01,fd_0.005\n" << std::setprecision(17);
  auto writeDirection=[&](const auto &x,int e) {
    dirs << x.name<<','<<e<<','<<x.layer<<','<<mesh.NodeRadius(e,0)*sem.lengthScaleMetres()<<','
         <<mesh.NodeRadius(e,mesh.NN()-1)*sem.lengthScaleMetres()<<','<<x.strength<<','<<x.selected<<','
         <<x.predicted<<','<<x.measured[0]<<','<<x.measured[1]<<','<<x.measured[2]<<'\n';
  };
  writeDirection(d,solid); writeDirection(b10,fluid); writeDirection(b100,fluid);
  const double f2mHz=freq.f(kFirstBin)/base.TimeNorm()*1000.0;
  const double f3mHz=freq.f(kFirstBin+1)/base.TimeNorm()*1000.0;
  const double f1mHz=freq.f1()/base.TimeNorm()*1000.0;
  const double fMaxmHz=freq.f2()/base.TimeNorm()*1000.0;
  const double cornersmHz[4]={freq.f11()/base.TimeNorm()*1000.0,
                              freq.f12()/base.TimeNorm()*1000.0,
                              freq.f21()/base.TimeNorm()*1000.0,
                              freq.f22()/base.TimeNorm()*1000.0};
  const double dampingRatePerSecond=freq.ep()/base.TimeNorm();
  const double dtSeconds=freq.dt()*base.TimeNorm();
  const double windowEndSeconds=(freq.nt()-1)*dtSeconds;
  std::ofstream cfg(out/"run.txt");
  cfg << std::setprecision(17) << "model="<<fs::absolute(modelPath)<<"\nobjective="<<result.objective<<"\n"
      << "source: depth="<<source.Depth()<<" km, lat="<<source.Latitude()<<" deg, lon="<<source.Longitude()<<" deg; tensor order (Mrr,Mrt,Mrp,Mtt,Mtp,Mpp)=("<<source.Mrr()<<','<<source.Mrt()<<','<<source.Mrp()<<','<<source.Mtt()<<','<<source.Mtp()<<','<<source.Mpp()<<") N m; receiver: depth="<<p.receiver_depth()<<" km, lat="<<receiver.first<<" deg, lon="<<receiver.second<<" deg\n"
      << "mode type="<<p.type()<<" (radial+spheroidal+toroidal); attenuation="<<p.attenuation()<<"; gravity formulation="<<p.gravitation()<<"; Cowling cutoff="<<p.cowlingFrequencyMhz()<<" mHz; output=displacement\n"
      << "degrees="<<p.lmin()<<".."<<p.lmax()<<"; NQ="<<mesh.NN()<<"; elements="<<mesh.NE()<<"; total one-sided GLL samples="<<mesh.NN()*mesh.NE()<<"; maximum element width="<<kMaxstep<<"R\n"
      << "configured frequency passband mHz=["<<f1mHz<<','<<fMaxmHz<<"]; configured taper-corner metadata mHz="<<cornersmHz[0]<<','<<cornersmHz[1]<<','<<cornersmHz[2]<<','<<cornersmHz[3]<<"; no frequency taper/filter applied; bins=["<<result.firstSolvedBin<<','<<result.endSolvedBin<<"); bin frequencies mHz="<<f2mHz<<','<<f3mHz<<"\n"
      << "nt="<<freq.nt()<<"; dt_s="<<dtSeconds<<"; timeNorm_s="<<base.TimeNorm()<<"; Hann=0.5*(1-cos(2*pi*n/(nt-1))); support_s=[0,"<<windowEndSeconds<<"]\n"
      << "damping epsilon_s^-1="<<dampingRatePerSecond<<" s^-1; solved complex frequency z=2*pi*f-i*epsilon\n"
      << "weights: q[vertical,"<<kFirstBin<<"]="<<kQBin2<<"; q[vertical,"<<kFirstBin+1<<"]="<<kQBin3<<"; SI displacement; full-gravity; fixed source/receiver/mesh\n"
      << "density bump: +100 kg/m^3*(1-x^2), solid element; fluid bulk amplitude study: +1e10 and +1e11 Pa*(1-x^2); x maps element endpoints to [-1,1]\n"
      << "direction locations: mantle element nearest 0.7R; fluid element nearest 0.3R; actual radii and selected gate row are written to directions.csv\n"
      << "fluid amplitude selection: compare 10 and 100 GPa using the same steps; 100 GPa is the declared gate direction to improve finite-difference signal against the objective baseline\n"
      << "FD derivative=[J(+a)-J(-a)]/(2a), a dimensionless multiplier, steps 0.02,0.01,0.005; prediction=integral kernel*deltaMaterial dr\n"
      << "normalization: L_m="<<sem.lengthScaleMetres()<<" m; P_Pa="<<sem.pressureScalePascals()<<" Pa; G_unit_SI="<<GUnitSI<<" m^3 kg^-1 s^-2; G_hat="<<GHat<<" (dimensionless); G_physical="<<physicalG<<" m^3 kg^-1 s^-2; densityScale_kgm3="<<base.DensityNorm()<<"\n"
      << "objective units: q-contracted processed displacement in m*s under the documented spectral inner product; derivative is per unit bump amplitude\n"
      << "measurement map: y=dtSI*RFFT[Hann*L_m*df_hat*exp(epsilon_hat*t_hat)*C2R(raw displacement spectrum)]; exponential damping compensation is applied once after the inverse real FFT and before the window\n"
      << "density representation: one-sided element-local GLL samples with polynomial interpolation; fixed mesh/interfaces, source, receiver, and absolute moduli; full gravity recomputed\n"
      << "scope: discrete derivative on this fixed mesh and selected [2,4) bins; no mesh or bandwidth convergence claim\n";
  std::ofstream json(out/"run.json");
  json << std::setprecision(17) << "{\n  \"model\": " << std::quoted(fs::absolute(modelPath).string())
       << ",\n  \"objective\": " << result.objective
       << ",\n  \"degree_min\": " << p.lmin() << ", \"degree_max\": " << p.lmax()
       << ",\n  \"nq\": " << mesh.NN() << ", \"element_count\": " << mesh.NE()
       << ",\n  \"one_sided_gll_sample_count\": " << mesh.NN()*mesh.NE()
       << ",\n  \"max_element_width_fraction_of_radius\": " << kMaxstep
       << ",\n  \"first_bin\": " << result.firstSolvedBin << ", \"end_bin_exclusive\": " << result.endSolvedBin
       << ",\n  \"frequency_bin_2_mhz\": " << f2mHz
       << ",\n  \"frequency_bin_3_mhz\": " << f3mHz << ",\n  \"nt\": " << freq.nt()
       << ",\n  \"frequency_min_mhz\": " << f1mHz << ", \"frequency_max_mhz\": " << fMaxmHz
       << ",\n  \"filter_corner_frequencies_mhz\": [" << cornersmHz[0] << ',' << cornersmHz[1] << ',' << cornersmHz[2] << ',' << cornersmHz[3] << ']'
       << ",\n  \"filter_corner_values_are_config_metadata_only\": true, \"frequency_taper_or_filter_applied\": false"
       << ",\n  \"mode_type\": " << p.type() << ", \"attenuation_enabled\": " << (p.attenuation() ? "true" : "false")
       << ",\n  \"gravity_formulation\": " << p.gravitation() << ", \"cowling_cutoff_mhz\": " << p.cowlingFrequencyMhz()
       << ",\n  \"dt_seconds\": " << dtSeconds << ",\n  \"time_norm_seconds\": " << base.TimeNorm()
       << ",\n  \"time_window_start_seconds\": 0, \"time_window_end_seconds\": " << windowEndSeconds
       << ",\n  \"damping_epsilon_per_second\": " << dampingRatePerSecond
       << ",\n  \"complex_frequency_convention\": \"z=2*pi*f-i*epsilon\""
       << ",\n  \"length_scale_metres\": " << sem.lengthScaleMetres()
       << ",\n  \"stress_scale_pascals\": " << sem.pressureScalePascals()
       << ",\n  \"density_scale_kg_m3\": " << base.DensityNorm()
       << ",\n  \"gravitational_constant_unit_m3_kg_s2\": " << GUnitSI
       << ",\n  \"gravitational_constant_nondimensional\": " << GHat
       << ",\n  \"physical_gravitational_constant_m3_kg_s2\": " << physicalG
       << ",\n  \"q_vertical_bin_2_real\": " << kQBin2.real() << ", \"q_vertical_bin_2_imag\": " << kQBin2.imag()
       << ",\n  \"q_vertical_bin_3_real\": " << kQBin3.real() << ", \"q_vertical_bin_3_imag\": " << kQBin3.imag()
       << ",\n  \"time_window\": \"Hann: 0.5*(1-cos(2*pi*n/(nt-1)))\", \"time_window_applied_once\": true"
       << ",\n  \"selected_solid_element\": " << solid << ", \"selected_fluid_element\": " << fluid
       << ",\n  \"solid_element_target_radius_fraction\": 0.7, \"fluid_element_target_radius_fraction\": 0.3"
       << ",\n  \"source_latitude_deg\": " << source.Latitude() << ", \"source_longitude_deg\": " << source.Longitude() << ", \"source_depth_km\": " << source.Depth()
       << ",\n  \"source_tensor_order\": [\"Mrr\",\"Mrt\",\"Mrp\",\"Mtt\",\"Mtp\",\"Mpp\"]"
       << ",\n  \"source_tensor_Nm\": [" << source.Mrr() << ',' << source.Mrt() << ',' << source.Mrp() << ',' << source.Mtt() << ',' << source.Mtp() << ',' << source.Mpp() << ']'
       << ", \"receiver_latitude_deg\": " << receiver.first << ", \"receiver_longitude_deg\": " << receiver.second << ", \"receiver_depth_km\": " << p.receiver_depth()
       << ",\n  \"density_bump_kg_m3\": 100.0, \"fluid_bulk_amplitude_study_pa\": [10000000000.0, 100000000000.0]"
       << ",\n  \"selected_fluid_bulk_gate_amplitude_pa\": 100000000000.0"
       << ",\n  \"fd_steps\": [0.02, 0.01, 0.005], \"relative_error_gate\": 1e-6"
       << ",\n  \"units\": \"objective: m s; density kernel: objective/((kg m^-3) m); bulk kernel: objective/(Pa m)\"\n}\n";
  fs::remove(paramPath);
  std::cout << "Wrote K2 kernel samples and directional checks to " << out << '\n';
  return 0;
} catch(const std::exception &e) { std::cerr << "sensitivity_kernels: "<<e.what()<<'\n'; return 2; }
