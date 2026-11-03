// Isolated Phase B adapter. The reviewed driver supplies physical SEM setup;
// its radial solver and original main are not called by this executable.
#define main phase_a_driver_main
#include "physical_driver.cpp"
#undef main

namespace phase_b {
using Record=dspec_global_gsh::PoleRecord;
using Reader=dspec_global_gsh::PoleResponseSpoolReader;
using Records=dspec_global_gsh::FrequencySideRecords;
using Frame=SPARSESPEC::windowed3d_detail::PoleFrame3D;
using Strain=SPARSESPEC::SpinStrainField;
struct Angle {double theta=0,phi=0,along=0;};
constexpr std::array<int,6> spins{0,0,1,-1,2,-2};

// Same geometry and spin phases as Windowed3DFieldEngine::tangentRotation.
std::array<double,4> localAngle(const Frame &f,const Angle &p) {
  const std::array<double,3> x{std::sin(p.theta)*std::cos(p.phi),std::sin(p.theta)*std::sin(p.phi),std::cos(p.theta)};
  std::array<double,3> local{};
  for(int j=0;j<3;++j)for(int i=0;i<3;++i)local[j]+=f.axes[j][i]*x[i];
  const double h=std::hypot(local[0],local[1]);
  const double t=std::atan2(h,local[2]),q=h==0?0:std::atan2(local[1],local[0]);
  const std::array<double,3> et{std::cos(t)*std::cos(q),std::cos(t)*std::sin(q),-std::sin(t)};
  std::array<double,3> g{};
  for(int i=0;i<3;++i)for(int j=0;j<3;++j)g[i]+=f.axes[j][i]*et[j];
  const std::array<double,3> gt{std::cos(p.theta)*std::cos(p.phi),std::cos(p.theta)*std::sin(p.phi),-std::sin(p.theta)},gp{-std::sin(p.phi),std::cos(p.phi),0};
  double a=0,b=0;for(int i=0;i<3;++i){a+=gt[i]*g[i];b+=gp[i]*g[i];}
  return {t,q,a,b};
}
Strain rotate(const std::array<Complex,6> &v,double a,double b) {
  const Complex c(a,-b),c2=c*c;
  return {v[0],v[1],c*v[2],std::conj(c)*v[3],c2*v[4],std::conj(c2)*v[5]};
}
struct Basis {
  std::vector<std::array<double,4>> geometry;
  std::map<std::pair<int,int>,Eigen::MatrixXd> matrices;
  std::vector<std::pair<int,int>> labels;
  Frame frame;
};
Basis makeBasis(const Record &record,int degree,const std::vector<Angle> &angles) {
  Basis out;out.labels=record.labels;out.frame=record.frame;
  for(auto a:angles)out.geometry.push_back(localAngle(record.frame,a));
  std::set<int> orders;for(auto [m,l]:record.labels){(void)l;orders.insert(m);}
  for(int spin=-2;spin<=2;++spin)for(int m:orders) {
    const int first=std::max({1,std::abs(spin),std::abs(m)});
    out.matrices[{spin,m}].resize(angles.size(),degree-first+1);
  }
  // Generate narrow-order Wigner values in bounded angle chunks. This is not
  // the O(L^3) full-order Stored table. Retain only the section's real bases.
  for(std::size_t begin=0;begin<angles.size();begin+=32) {
    const std::size_t count=std::min<std::size_t>(32,angles.size()-begin);
    std::vector<double> theta;for(std::size_t a=0;a<count;++a)theta.push_back(out.geometry[begin+a][0]);
    SPARSESPEC::Windowed3DFieldEngine::CachedAngularWigner w(degree,std::min(2,degree),std::min(2,degree),theta);
    for(auto &[key,matrix]:out.matrices) {
      const auto [spin,m]=key;const int first=std::max({1,std::abs(spin),std::abs(m)});
      if(std::abs(spin)>degree)continue;
      for(std::size_t a=0;a<count;++a)for(int l=first;l<=degree;++l)
        matrix(begin+a,l-first)=w[spin,static_cast<int>(a)][l,m];
    }
  }
  return out;
}
std::vector<std::array<Complex,6>> synthesize(const Basis &basis,
    const std::vector<Records> &records,bool adjoint,int degree) {
  const std::size_t sides=records.size(),angles=basis.geometry.size();
  std::vector<std::array<Complex,6>> values(angles*sides);
  for(std::size_t component=0;component<6;++component) {
    for(const auto &[key,matrix]:basis.matrices) {
      const auto [spin,m]=key;if(spin!=spins[component]||matrix.cols()==0)continue;
      const int first=std::max({1,std::abs(spin),std::abs(m)});
      Eigen::MatrixXcd coefficients=Eigen::MatrixXcd::Zero(degree-first+1,sides);
      bool nonzero=false;
      for(std::size_t r=0;r<sides;++r) {
        const auto &rec=adjoint?records[r].unitAdjoint:records[r].forward;
        if(rec.labels!=basis.labels||rec.frame.axes!=basis.frame.axes)
          throw std::runtime_error("section record layout/frame differs from the cached basis");
        for(std::size_t q=0;q<rec.labels.size();++q) {
          const auto [order,l]=rec.labels[q];if(order!=m||l<first)continue;
          // Fluid divergence collapse is handled by the caller, after both
          // spin-zero components have been synthesized with the same basis.
          coefficients(l-first,r)=rec.rawStrain[q][component];
          nonzero=nonzero||rec.rawStrain[q][component]!=Complex{};
        }
      }
      if(!nonzero)continue;
      const auto sums=dsfast::synthesize(matrix,coefficients,dsfast::DotBackend::Eigen);
      for(std::size_t a=0;a<angles;++a) {
        const Complex phase=std::exp(Complex(0,m*basis.geometry[a][1]));
        for(std::size_t r=0;r<sides;++r)values[a*sides+r][component]+=phase*sums(a,r);
      }
    }
  }
  return values;
}
std::vector<Angle> readAngles(const fs::path &path) {
  std::ifstream in(path);std::string line;std::getline(in,line);
  if(line!="theta_rad,phi_rad,along_deg")throw std::runtime_error("invalid section angle header");
  std::vector<Angle> out;
  while(std::getline(in,line)){std::replace(line.begin(),line.end(),',',' ');std::istringstream s(line);Angle a;
    if(!(s>>a.theta>>a.phi>>a.along)||!std::isfinite(a.theta)||!std::isfinite(a.phi)||!std::isfinite(a.along)||a.theta<0||a.theta>kPi)
      throw std::runtime_error("invalid section angle");out.push_back(a);}
  if(out.empty()||out.size()>6145)throw std::runtime_error("section angle count must be 1..6145");return out;
}
void section(const fs::path &spool,const Full1D::SEM &sem,const fs::path &out,
    const fs::path &anglePath,std::size_t workers) {
  const auto started=Clock::now();Reader metadata(spool.string());const int degree=metadata.lmax();
  const auto angles=readAngles(anglePath);const auto &sides=metadata.sides();
  const std::size_t bytes=angles.size()*sides.size()*6*sizeof(double);
  if(bytes>128ull*1024*1024||workers==0||workers>24||sides.size()!=179)
    throw std::runtime_error("section replay exceeds bounded output/worker/radial scope");
  auto first=metadata.read(0,0);
  std::cout<<"Building section bases "<<angles.size()<<" angles degree "<<degree<<std::endl;
  const auto forward=makeBasis(first.forward,degree,angles),adjoint=makeBasis(first.unitAdjoint,degree,angles);
  std::uint64_t basisBytes=0;for(const auto *b:{&forward,&adjoint})for(const auto &[key,v]:b->matrices){(void)key;basisBytes+=v.size()*sizeof(double);}
  const double setupSeconds=std::chrono::duration<double>(Clock::now()-started).count();
  std::cout<<"Section bases ready bytes="<<basisBytes<<" wall_s="<<setupSeconds<<std::endl;
  const fs::path binary=out/"section.bin";
  if(fs::exists(binary))throw std::runtime_error("section output exists");
  {std::ofstream f(binary,std::ios::binary);f.seekp(bytes-1);f.put('\0');if(!f)throw std::runtime_error("section allocation failed");}
  std::atomic<std::size_t> completed{0};std::mutex logMutex;
  constexpr std::size_t radialTile=16;
  const std::size_t tiles=(sides.size()+radialTile-1)/radialTile;
  dspec_global_gsh::runOwnedTiles(tiles,workers,[&](std::size_t tile,std::size_t worker){
    (void)worker;Reader reader(spool.string());
    const std::size_t begin=tile*radialTile,count=std::min(radialTile,sides.size()-begin);
    std::vector<std::array<double,6>> accumulator(angles.size()*count);
    for(std::size_t bin=0;bin<reader.frequencyCount();++bin) {
      const auto audit=reader.frequencyAudit(bin);
      if(audit.basisReconstructionRelative>1e-12)throw std::runtime_error("receiver basis reconstruction gate failed");
      if(audit.adjointMixing==Complex{})continue;
      std::vector<Records> records;records.reserve(count);
      for(std::size_t r=0;r<count;++r)records.push_back(reader.read(bin,begin+r));
      const auto f=synthesize(forward,records,false,degree),a=synthesize(adjoint,records,true,degree);
      for(std::size_t p=0;p<angles.size();++p)for(std::size_t r=0;r<count;++r) {
        const std::size_t index=p*count+r;auto fv=f[index],av=a[index];const auto &l0=records[r].l0;
        auto fs=rotate(fv,forward.geometry[p][2],forward.geometry[p][3]);
        auto as=rotate(av,adjoint.geometry[p][2],adjoint.geometry[p][3]);
        fs.rr+=l0.forwardRr;fs.horizontalTrace+=l0.forwardHorizontalTrace;
        as.rr+=l0.adjointRr;as.horizontalTrace+=l0.adjointHorizontalTrace;
        const auto mix=audit.adjointMixing;
        as={as.rr*mix,as.horizontalTrace*mix,as.radialPlus*mix,as.radialMinus*mix,as.horizontalPlus*mix,as.horizontalMinus*mix};
        SPARSESPEC::KernelField3DSample sample;sample.radialElement=sides[begin+r].element;sample.fluid=sides[begin+r].fluid;
        sample.forwardRawStrain=fs;sample.adjointDualStrain=as;
        const auto k=SPARSESPEC::Windowed3DKernelEngine::pointKernel(sem,sample);
        const std::array<double,6> values{k.A,k.C,k.F,k.L,k.N,k.kappa};
        for(std::size_t c=0;c<6;++c){if(!std::isfinite(values[c]))throw std::runtime_error("nonfinite section kernel");accumulator[index][c]+=values[c];}
      }
    }
    std::fstream f(binary,std::ios::binary|std::ios::in|std::ios::out);
    std::vector<std::array<double,6>> row(angles.size());
    for(std::size_t r=0;r<count;++r){for(std::size_t p=0;p<angles.size();++p)row[p]=accumulator[p*count+r];
      f.seekp((begin+r)*angles.size()*sizeof(row.front()));f.write(reinterpret_cast<const char*>(row.data()),row.size()*sizeof(row.front()));}
    if(!f)throw std::runtime_error("section write failed");f.close();
    const auto n=++completed;std::lock_guard<std::mutex> lock(logMutex);std::cout<<"Section radial tile "<<n<<'/'<<tiles<<" completed"<<std::endl;
  });
  fs::copy_file(anglePath,out/"section_angles.csv");
  std::ofstream csv(out/"section_radial_sides.csv");csv<<"native_side_index,element,node,radius_m,width_m,fluid\n"<<std::setprecision(17);
  for(std::size_t r=0;r<sides.size();++r){const auto &s=sides[r];csv<<r<<','<<s.element<<','<<s.node<<','<<s.radiusMetres<<','<<s.widthMetres<<','<<s.fluid<<'\n';}
  std::ofstream m(out/"section.json");m<<std::setprecision(17)<<"{\"degree\":"<<degree<<",\"angles\":"<<angles.size()<<",\"radial_sides\":"<<sides.size()<<",\"frequency_records\":"<<metadata.frequencyCount()<<",\"channels\":[\"A\",\"C\",\"F\",\"L\",\"N\",\"kappa\"],\"layout\":\"radial-side,angle,channel; float64 native endian\",\"source_contacts_included\":false,\"full_global_volume\":false,\"basis_bytes\":"<<basisBytes<<",\"basis_setup_s\":"<<setupSeconds<<",\"elapsed_s\":"<<std::chrono::duration<double>(Clock::now()-started).count()<<",\"workers_requested\":"<<workers<<",\"workers_actual\":"<<std::min(workers,tiles)<<",\"output_bytes\":"<<bytes<<"}\n";
}
} // namespace phase_b

int main(int argc,char **argv) try {
  if(argc!=8)throw std::runtime_error("usage: figure-replay section|slices MODEL SPOOL EMPTY_OUTPUT ANGLES_CSV|SIDE_INDICES THREADS");
  // argc includes seven supplied arguments; reserved final argument is the
  // immutable degree assertion, preventing accidental cross-degree replay.
  const std::string mode=argv[1];const fs::path spool=argv[3],output=argv[4];
  dspec_global_gsh::PoleResponseSpoolReader metadata(spool.string());
  if(std::stoi(argv[7])!=metadata.lmax())throw std::runtime_error("spool degree assertion mismatch");
  if(fs::exists(output))throw std::runtime_error("EMPTY_OUTPUT already exists");fs::create_directories(output);
  Options o;o.model=argv[2];o.output=output;o.lmax=metadata.lmax();o.semDegree=6;o.qex=1;o.radialRefinement=.05;
  Eigen::setNbThreads(1);auto context=prepare(o);
  const auto traces=nativeTraces(*context.sem);if(traces.size()!=metadata.sides().size())throw std::runtime_error("SEM/spool side count mismatch");
  for(std::size_t r=0;r<traces.size();++r){const auto &a=traces[r];const auto &b=metadata.sides()[r];
    if(a.element!=b.element||a.node!=b.node||a.radius!=b.radiusMetres||a.fluid!=b.fluid)throw std::runtime_error("SEM/spool native side mismatch");}
  if(mode=="section")phase_b::section(spool,*context.sem,output,argv[5],std::stoul(argv[6]));
  else if(mode=="slices") {
    dspec_global_gsh::ConsumerOptions options;options.workers=1;options.gshThreads=std::stoul(argv[6]);options.storedWignerValues=false;
    options.maxRadialSides=3;options.maxOutputGiB=4;options.maxWorkingGiB=64;options.transformScratchGiBPerWorker=8;
    std::istringstream s(argv[5]);std::string token;while(std::getline(s,token,','))options.radialSideIndices.push_back(std::stoul(token));
    if(options.radialSideIndices!=std::vector<std::size_t>{147,130,114}&&options.radialSideIndices!=std::vector<std::size_t>{59,100})throw std::runtime_error("slice selection outside approved production/control scope");
    dspec_global_gsh::runPoleResponseSpoolKernel(spool.string(),*context.sem,(output/"slices.bin").string(),options);
  } else throw std::runtime_error("mode must be section or slices");
  return 0;
} catch(const std::exception &e){std::cerr<<"Phase B replay failed: "<<e.what()<<'\n';return 1;}
