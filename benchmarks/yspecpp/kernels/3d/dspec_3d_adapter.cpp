#include <DSpecM1D/All>
#include <DSpecM1D/ModelInput>
#include <DSpecM1D/Windowed3DKernels>
#include <DSpecM1D/src/InputParser.h>
#include <DSpecM1D/src/NormClass.h>
#include <DSpecM1D/src/SourceInfo.h>

#include <array>
#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numbers>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;
using Complex = std::complex<double>;
constexpr int nt=128, bins=65, firstBin=1, endBin=4, nangTheta=8, nangPhi=24;
constexpr double r0=4200000.0, r1=6200000.0, sourceRadius=6321000.0;
constexpr double sourceLocalLower=6311000.0, sourceLocalUpper=6331000.0;
constexpr std::array<double,8> glx{-0.9602898564975363,-0.7966664774136267,-0.5255324099163290,-0.1834346424956498,0.1834346424956498,0.5255324099163290,0.7966664774136267,0.9602898564975363};
constexpr std::array<double,8> glw{0.1012285362903763,0.2223810344533745,0.3137066458778873,0.3626837833783620,0.3626837833783620,0.3137066458778873,0.2223810344533745,0.1012285362903763};
constexpr std::array<const char*,5> modNames{"A","C","F","L","N"};
constexpr std::array<double,5> baselineModuliPa{352960000000.0,352960000000.0,129602500000.0,111678750000.0,111678750000.0};
std::string paramsText(const fs::path& model,const fs::path& out){
  return "\""+(out/"unused").string()+"\"\n\""+model.string()+"\"\n4\n0\n2\n0\n0\n1e-5\n0\n3\n0.1\n0.8\n30\n20\n0\n0\n0\n0\n50\n10\n20\n1e19\n2e18\n-3e18\n-5e18\n4e18\n-5e18\n3\n1\n45 90\n";
}
Eigen::VectorXd hann(){ Eigen::VectorXd w=Eigen::VectorXd::Zero(nt); for(int i=0;i<=90;++i) w[i]=0.5-0.5*std::cos(2*std::numbers::pi*i/90.0); return w; }
Eigen::MatrixXcd qweights(int component){ Eigen::MatrixXcd q=Eigen::MatrixXcd::Zero(3,bins); q(component,1)=Complex(0.7,0.2); q(component,2)=Complex(-0.25,0.45); return q; }
double chi(double theta,double phi){
  const double x=std::cos(theta),s=std::sin(theta);
  const double n21=std::sqrt(5.0/(24.0*std::numbers::pi));
  const double n32=std::sqrt(7.0/(480.0*std::numbers::pi));
  return -3.0*x*s*n21*std::cos(phi)+7.5*x*(1.0-x*x)*n32*std::cos(2.0*phi);
}
double radialWeight(double r){ const double tol=8.0*std::numeric_limits<double>::epsilon()*std::max(std::abs(r0),std::abs(r1)); if(r<=r0+tol||r>=r1-tol)return 0.0; return std::pow(std::sin(std::numbers::pi*(r-r0)/(r1-r0)),2); }
double sourceLocalWeight(double r){
  const double tol=8.0*std::numeric_limits<double>::epsilon()*
      std::max(std::abs(sourceLocalLower),std::abs(sourceLocalUpper));
  if(r<=sourceLocalLower+tol||r>=sourceLocalUpper-tol)return 0.0;
  return std::pow(std::cos(std::numbers::pi*(r-sourceRadius)/
                           (sourceLocalUpper-sourceLocalLower)),2);
}
Eigen::SparseMatrix<Complex> oracleH0(const Full1D::SEM&sem,SPARSESPEC::SensitivitySolveFamily family,int l,Complex z){
  Eigen::SparseMatrix<double> h,p;
  if(family==SPARSESPEC::SensitivitySolveFamily::Radial){h=sem.hR();p=sem.pR();}
  else if(family==SPARSESPEC::SensitivitySolveFamily::Spheroidal){h=sem.hS(l);p=sem.pS(l);}
  else {h=sem.hTk(l);p=sem.pTk(l);}
  Eigen::SparseMatrix<Complex> out=h.cast<Complex>()-z*z*p.cast<Complex>();out.makeCompressed();return out;
}
void writeMatrixMarket(const fs::path&path,const Eigen::SparseMatrix<Complex>&matrix){
  std::ofstream f(path);f<<"%%MatrixMarket matrix coordinate complex general\n"<<matrix.rows()<<' '<<matrix.cols()<<' '<<matrix.nonZeros()<<'\n'<<std::setprecision(17);
  for(int col=0;col<matrix.outerSize();++col)for(Eigen::SparseMatrix<Complex>::InnerIterator it(matrix,col);it;++it)
    f<<it.row()+1<<' '<<it.col()+1<<' '<<it.value().real()<<' '<<it.value().imag()<<'\n';
}
void exportSourceElementOperators(Full1D::SEM&sem,SourceInfo::EarthquakeCMT&source,
                                 const SpectraSolver::FreqFull&freq,const fs::path&dir){
  const int element=sem.sourceElement(source);
  for(int k=firstBin;k<endBin;++k){
    const Complex z=freq.w(k)+Complex(0.0,-freq.ep());
    const auto writeLocal=[&](const std::string&family,int degree,const Eigen::MatrixXcd&local){
      Eigen::SparseMatrix<Complex> sparse=local.sparseView();sparse.makeCompressed();
      writeMatrixMarket(dir/("source_local_H0_k"+std::to_string(k)+"_"+family+"_l"+std::to_string(degree)+".mtx"),sparse);
    };
    writeLocal("R",0,Full1D::source_treatment_detail::radialElementOperator(sem,element,z));
    for(int l=1;l<=sem.maxDegree();++l){
      writeLocal("S",l,Full1D::source_treatment_detail::spheroidalElementOperator(sem,element,l,z,{},false));
      writeLocal("T",l,Full1D::source_treatment_detail::toroidalElementOperator(sem,element,l,z));
    }
  }
}
void exportOracleBlocks(Full1D::SEM&sem,SourceInfo::EarthquakeCMT&source,InputParameters&parameters,
                       const SpectraSolver::FreqFull&freq,double gravitySI,const fs::path&out){
  const fs::path dir=out/"oracle";fs::create_directories(dir);
  std::ofstream modes(dir/"mesh_dofs.csv");modes<<"family,degree,element,node,field,index\n"<<std::setprecision(17);
  for(int e=0;e<static_cast<int>(sem.mesh().NE());++e)for(int q=0;q<sem.mesh().NN();++q){
    modes<<"R,0,"<<e<<','<<q<<",U,"<<sem.ltgR(0,e,q)<<'\n';
    modes<<"R,0,"<<e<<','<<q<<",phi,"<<sem.ltgR(1,e,q)<<'\n';
    for(int l=1;l<=parameters.lmax();++l){
      modes<<"S,"<<l<<','<<e<<','<<q<<",U,"<<sem.ltgS(0,e,q)<<'\n';
      modes<<"S,"<<l<<','<<e<<','<<q<<",V,"<<sem.ltgS(1,e,q)<<'\n';
      modes<<"S,"<<l<<','<<e<<','<<q<<",phi,"<<sem.ltgS(2,e,q)<<'\n';
      modes<<"T,"<<l<<','<<e<<','<<q<<",W,"<<sem.ltgT(e,q)<<'\n';
    }
  }
  std::ofstream src(dir/"source.csv");src<<"bin,family,degree,order,dof,rhs_real,rhs_imag,lift_real,lift_imag\n"<<std::setprecision(17);
  std::ofstream sourceBasis(dir/"source_basis.csv");
  sourceBasis<<"bin,family,degree,order,source_column,dof,coefficient_real,coefficient_imag,force_real,force_imag,lift_real,lift_imag\n"<<std::setprecision(17);
  std::ofstream receiver(dir/"receiver.csv");receiver<<"bin,family,degree,order,dof,component,coefficient_real,coefficient_imag\n"<<std::setprecision(17);
  std::ofstream solvedReceiver(dir/"solved_receiver.csv");
  solvedReceiver<<"bin,family,degree,order,component,coefficient_real,coefficient_imag,residual_relative\n"<<std::setprecision(17);
  const int nfreq=freq.nt()/2+1;Eigen::VectorXcd zeroAdj=Eigen::VectorXcd::Zero(3);
  for(int k=firstBin;k<endBin;++k){
    const Complex z=freq.w(k)+Complex(0.0,-freq.ep());
    auto workspace=SPARSESPEC::Windowed3DFieldEngine::prepareFrequencyWorkspace(sem,source,parameters,z,zeroAdj);
    for(const auto&block:workspace.modalBlocks()){
      const int l=block.degree;Eigen::MatrixXcd force,coeff;
      std::string family;
      Eigen::MatrixXcd receiverBase;
      if(block.family==SPARSESPEC::SensitivitySolveFamily::Radial){
        family="R";force=sem.sourceTreatedRadial(source,z).force;coeff=Eigen::MatrixXcd::Ones(1,1);
        receiverBase=Eigen::MatrixXcd::Zero(3,sem.hR().rows());
        receiverBase.row(0)=sem.rvZR(parameters,0).transpose();
      } else if(block.family==SPARSESPEC::SensitivitySolveFamily::Spheroidal){
        family="S";force=sem.sourceTreatedSpheroidal(source,l,z).force;coeff=block.sourceCoefficients;
        const auto elems=sem.receiverElements(parameters,false);
        const auto local=sem.rvBaseFull(parameters,l,false);
        const int low=sem.ltgS(0,elems.front(),0);
        receiverBase=Eigen::MatrixXcd::Zero(local.rows(),sem.hS(l).rows());
        receiverBase.block(0,low,local.rows(),local.cols())=local;
      } else {
        family="T";force=sem.sourceTreatedToroidal(source,l,z).force;coeff=block.sourceCoefficients;
        const auto elems=sem.receiverElements(parameters,true);
        const auto local=sem.rvBaseFullT(parameters,l);
        const int low=sem.ltgT(elems.front(),0);
        receiverBase=Eigen::MatrixXcd::Zero(local.rows(),sem.hTk(l).rows());
        receiverBase.block(0,low,local.rows(),local.cols())=local;
      }
      if(force.cols()!=coeff.cols()||coeff.rows()!=(family=="R"?1:2*l+1))throw std::runtime_error("oracle source basis dimensions differ from modal coefficient table");
      const auto H=oracleH0(sem,block.family,l,z);
      writeMatrixMarket(dir/("H0_k"+std::to_string(k)+"_"+family+"_l"+std::to_string(l)+".mtx"),H);
      for(int im=0;im<coeff.rows();++im){
        const Eigen::VectorXcd c=coeff.row(im).transpose();
        const Eigen::VectorXcd rhs=force*c,lift=block.sourceLift*c;
        const Eigen::VectorXcd solved=block.forward*c;
        const double residual=(H*solved-rhs).norm()/std::max(rhs.norm(),1.0e-300);
        const int m=family=="R"?0:im-l;
        for(Eigen::Index col=0;col<force.cols();++col){
          const Complex coefficient=coeff.rows()==1?Complex(1.0,0.0):coeff(im,col);
          for(Eigen::Index d=0;d<force.rows();++d)
            if(force(d,col)!=Complex{}||block.sourceLift(d,col)!=Complex{})
              sourceBasis<<k<<','<<family<<','<<l<<','<<m<<','<<col<<','<<d<<','
                         <<coefficient.real()<<','<<coefficient.imag()<<','
                         <<force(d,col).real()<<','<<force(d,col).imag()<<','
                         <<block.sourceLift(d,col).real()<<','<<block.sourceLift(d,col).imag()<<'\n';
        }
        for(Eigen::Index d=0;d<rhs.size();++d)src<<k<<','<<family<<','<<l<<','<<m<<','<<d<<','<<rhs[d].real()<<','<<rhs[d].imag()<<','<<lift[d].real()<<','<<lift[d].imag()<<'\n';
        const Eigen::VectorXcd displacement=receiverBase*solved;
        for(int component=0;component<3;++component){
          Complex value=displacement[component];
          if(family!="R")value*=block.modeReceiver(component,im);
          solvedReceiver<<k<<','<<family<<','<<l<<','<<m<<','<<component<<','
                        <<value.real()<<','<<value.imag()<<','<<residual<<'\n';
        }
        for(int component=0;component<3;++component)for(Eigen::Index d=0;d<receiverBase.cols();++d){
          Complex value=receiverBase(component,d);
          if(family!="R")value*=block.modeReceiver(component,im);
          if(value!=Complex{})receiver<<k<<','<<family<<','<<l<<','<<m<<','<<d<<','<<component<<','<<value.real()<<','<<value.imag()<<'\n';
        }
      }
    }
  }
  exportSourceElementOperators(sem,source,freq,dir);
  std::ofstream meta(dir/"meta.json");meta<<std::setprecision(17)<<"{\"frequency_bins\":[1,2,3],\"length_scale_m\":"<<sem.lengthScaleMetres()<<",\"pressure_scale_pa\":"<<sem.pressureScalePascals()<<",\"time_scale_s\":"<<sem.timeScaleSeconds()<<",\"physical_G\":"<<sem.gravityCoefficient()*gravitySI<<",\"radius_m\":"<<sem.mesh().PR()*sem.lengthScaleMetres()<<",\"source_radius_m\":"<<sem.sourceRadius()*sem.lengthScaleMetres()<<",\"source_element\":"<<sem.sourceElement(source)<<",\"source_local_element\":"<<sem.sourceElement(source)<<",\"source_local_nq\":"<<sem.mesh().NN()<<",\"source_element_bounds_m\":["<<sem.mesh().ELR(sem.sourceElement(source))*sem.lengthScaleMetres()<<','<<sem.mesh().EUR(sem.sourceElement(source))*sem.lengthScaleMetres()<<"],\"lmax\":"<<parameters.lmax()<<",\"nt\":"<<freq.nt()<<",\"dt_model_s\":"<<freq.dt()<<",\"df_model\":"<<freq.df()<<",\"epsilon_model\":"<<freq.ep()<<"}\n";
}
int main(int argc,char**argv) try {
  if(argc!=5) throw std::runtime_error("usage: phase_d_adapter MODEL OUTPUT MAXSTEP NQ");
  const fs::path modelPath=fs::absolute(argv[1]),out=fs::absolute(argv[2]);
  if(!fs::exists(out)||!fs::is_directory(out)||fs::directory_iterator(out)!=fs::directory_iterator{}) throw std::runtime_error("adapter output must be a newly created empty directory");
  const double maxstep=std::stod(argv[3]); const int nq=std::stoi(argv[4]);
  prem_norm<double> norm; using Base=EarthModels::ModelInput<double>; Base base(modelPath.string(),norm);
  auto pf=out/"params.txt"; {std::ofstream f(pf); f<<paramsText(modelPath,out);}
  InputParametersNew p(pf.string(),nq,1,maxstep,1.0,0.1,0.0,1);
  auto &parameters=p.inputParameters(); auto &freq=p.freqFull(); SourceInfo::EarthquakeCMT source(parameters);
  Full1D::SEM sem(base,maxstep,nq,3,source); const auto window=hann();
  if(freq.nt()!=nt||freq.i1()>firstBin||freq.i2()<endBin) throw std::runtime_error("resolved DSpec grid does not cover required bins 1..3");
  std::vector<SPARSESPEC::WeightedField3DPoint> quadrature,sourceLocalQuadrature;
  std::vector<SPARSESPEC::Field3DPoint> points,sourceLocalPoints;
  std::vector<std::array<double,2>> angles;
  std::vector<std::array<double,4>> radialNodes,meshNodes;
  std::vector<int> sourceLocalElements;
  int sourceLocalLowerNodes=0;
  int sourceLocalUpperNodes=0;
  for(int it=0;it<nangTheta;++it) for(int ip=0;ip<nangPhi;++ip){
    const double theta=std::acos(glx[it]),phi=2.0*std::numbers::pi*ip/nangPhi; angles.push_back({theta,phi});
  }
  for(int e=0;e<static_cast<int>(sem.mesh().NE());++e) for(int q=0;q<sem.mesh().NN();++q){
    const double radius=sem.mesh().NodeRadius(e,q)*sem.lengthScaleMetres(); const double wr=radialWeight(radius);
    const double dr=0.5*sem.mesh().EW(e)*const_cast<EarthMesh::RadialMesh&>(sem.mesh()).GLL().W(q)*sem.lengthScaleMetres();
    meshNodes.push_back({static_cast<double>(e),static_cast<double>(q),radius,dr});
    const double elementLower=sem.mesh().ELR(e)*sem.lengthScaleMetres();
    const double elementUpper=sem.mesh().EUR(e)*sem.lengthScaleMetres();
    const bool inClosedSupport=elementLower>=r0-8.0*std::numeric_limits<double>::epsilon()*r0 &&
                               elementUpper<=r1+8.0*std::numeric_limits<double>::epsilon()*r1;
    if(!inClosedSupport) continue;
    radialNodes.push_back({static_cast<double>(e),static_cast<double>(q),radius,dr});
    for(int it=0;it<nangTheta;++it) for(int ip=0;ip<nangPhi;++ip){
      const double theta=std::acos(glx[it]),phi=2.0*std::numbers::pi*ip/nangPhi;
      SPARSESPEC::Field3DPoint point{radius,theta,phi,e};
      quadrature.push_back({point,dr,glw[it]*(2.0*std::numbers::pi/nangPhi)}); points.push_back(point);
    }
  }
  const double sourceSupportTol=8.0*std::numeric_limits<double>::epsilon()*
      std::max(std::abs(sourceLocalLower),std::abs(sourceLocalUpper));
  for(int e=0;e<static_cast<int>(sem.mesh().NE());++e){
    const double elementLower=sem.mesh().ELR(e)*sem.lengthScaleMetres();
    const double elementUpper=sem.mesh().EUR(e)*sem.lengthScaleMetres();
    if(elementLower<sourceLocalLower-sourceSupportTol||
       elementUpper>sourceLocalUpper+sourceSupportTol)continue;
    sourceLocalElements.push_back(e);
    for(int q=0;q<sem.mesh().NN();++q){
      const double radius=sem.mesh().NodeRadius(e,q)*sem.lengthScaleMetres();
      if(radius<sourceRadius)++sourceLocalLowerNodes;
      else if(radius>sourceRadius)++sourceLocalUpperNodes;
      const double dr=0.5*sem.mesh().EW(e)*
          const_cast<EarthMesh::RadialMesh&>(sem.mesh()).GLL().W(q)*sem.lengthScaleMetres();
      for(int it=0;it<nangTheta;++it)for(int ip=0;ip<nangPhi;++ip){
        const double theta=std::acos(glx[it]),phi=2.0*std::numbers::pi*ip/nangPhi;
        SPARSESPEC::Field3DPoint point{radius,theta,phi,e};
        sourceLocalQuadrature.push_back({point,dr,glw[it]*(2.0*std::numbers::pi/nangPhi)});
        sourceLocalPoints.push_back(point);
      }
    }
  }
  if(sourceLocalElements.empty()||sourceLocalLowerNodes==0||sourceLocalUpperNodes==0)
    throw std::runtime_error("source-local closed support lacks elements or interior nodes on both sides of source");
  if(quadrature.empty()) throw std::runtime_error("regular profile has no resolved quadrature points");
  std::ofstream spectra(out/"spectra.csv"),processed(out/"processed.csv");
  std::ofstream profiles(out/"profiles.csv"),integrals(out/"integrals.csv");
  std::ofstream contacts(out/"contacts.csv"),radial(out/"radial_nodes.csv");
  std::ofstream mesh(out/"mesh_nodes.csv"),angular(out/"angles.csv");
  std::ofstream sourceLocalProfiles(out/"source_local_profiles.csv");
  std::ofstream sourceLocalIntegrals(out/"source_local_integrals.csv");
  radial<<"element,node,radius_m,dr_m\n"<<std::setprecision(17); for(const auto&v:radialNodes) radial<<v[0]<<','<<v[1]<<','<<v[2]<<','<<v[3]<<'\n';
  mesh<<"element,node,radius_m,dr_m,regular_weight,source_local_weight\n"<<std::setprecision(17); for(const auto&v:meshNodes){const double r=v[2],wr=radialWeight(r),x=(r-sourceRadius)/10000.0,ws=std::abs(x)<=1.0?std::pow(std::cos(std::numbers::pi*x/2.0),2):0.0;mesh<<v[0]<<','<<v[1]<<','<<r<<','<<v[3]<<','<<wr<<','<<ws<<'\n';}
  angular<<"theta_rad,phi_rad,dOmega,chi\n"<<std::setprecision(17); for(int it=0;it<nangTheta;++it) for(int ip=0;ip<nangPhi;++ip){double th=std::acos(glx[it]),ph=2.0*std::numbers::pi*ip/nangPhi; angular<<th<<','<<ph<<','<<glw[it]*2.0*std::numbers::pi/nangPhi<<','<<chi(th,ph)<<'\n';}
  spectra<<"objective,component,bin,frequency_hz,raw_real_m_s,raw_imag_m_s\n"<<std::setprecision(17);
  processed<<"objective,component,bin,real_m_s,imag_m_s\n"<<std::setprecision(17);
  profiles<<"objective,modulus,radius_m,theta_rad,phi_rad,chi,w,dV_m3,deltaX_Pa,kernel_per_Pa_m3\n"<<std::setprecision(17);
  integrals<<"objective,modulus,S_m_s,J0_m_s,Vsupport_m3,Kscale_per_Pa_m3,I_regular_m_s,Kweighted_per_Pa_m3\n"<<std::setprecision(17);
  contacts<<"objective,modulus,theta_rad,phi_rad,chi,contact_per_Pa,deltaX_Pa,dOmega,contribution_m_s\n"<<std::setprecision(17);
  sourceLocalProfiles
      <<"objective,modulus,element,radius_m,theta_rad,phi_rad,chi,w,dV_m3,"
        "deltaX_Pa,kernel_per_Pa_m3,contribution_m_s_per_h\n"
      <<std::setprecision(17);
  sourceLocalIntegrals
      <<"objective,modulus,S_m_s,Vsupport_m3,Kscale_per_Pa_m3,"
        "I_volume_m_s_per_h,I_contact_m_s_per_h,Kpred_total_m_s_per_h\n"
      <<std::setprecision(17);
  // Keep the regular samples first so their existing profile indices stay fixed.
  auto evaluationPoints=points;
  evaluationPoints.insert(evaluationPoints.end(),sourceLocalPoints.begin(),sourceLocalPoints.end());
  const double sourceLocalSupportVolume=(4.0*std::numbers::pi/3.0)*
      (sourceLocalUpper*sourceLocalUpper*sourceLocalUpper-
       sourceLocalLower*sourceLocalLower*sourceLocalLower);
  for(int obs=0;obs<2;++obs){
    auto q=qweights(obs==0?0:1);
    const auto result=SPARSESPEC::Windowed3DKernelEngine::evaluate(
        sem,source,parameters,freq,window,q,firstBin,endBin,evaluationPoints,angles);
    const std::string label=obs==0?"J_Z":"J_N";
    double qnorm=q.norm(),processedNorm=result.processedSISpectrum.norm(); const double S=qnorm*processedNorm;
    if(!(std::isfinite(S)&&S>0.0)) throw std::runtime_error("baseline objective scale S is zero or invalid");
    for(int c=0;c<3;++c) for(int k=0;k<bins;++k){
      const auto u=result.rawDisplacementSpectrum(c,k)*sem.lengthScaleMetres()*sem.timeScaleSeconds();
      spectra<<label<<','<<(c==0?"Z":c==1?"N":"E")<<','<<k<<','<<freq.f(k)/freq.timeNorm()<<','<<u.real()<<','<<u.imag()<<'\n';
      const auto v=result.processedSISpectrum(c,k); processed<<label<<','<<(c==0?"Z":c==1?"N":"E")<<','<<k<<','<<v.real()<<','<<v.imag()<<'\n';
    }
    const double dOmega=2.0*std::numbers::pi/nangPhi;
    double Vsupport=0.0; std::array<double,5> I{},W{},K2{};
    for(std::size_t i=0;i<quadrature.size();++i){
      const auto&qpt=quadrature[i]; const double r=qpt.point.radiusMetres,wr=radialWeight(r),ch=chi(qpt.point.thetaRadians,qpt.point.phiRadians);
      const double dV=r*r*qpt.radialWeightMetres*qpt.angularWeightSteradians; Vsupport+=dV;
      const auto&k=result.pointKernels[i]; const std::array<double,5> kval{k.A,k.C,k.F,k.L,k.N};
      for(int m=0;m<5;++m){ const double delta=baselineModuliPa[m]*wr*ch; I[m]+=kval[m]*delta*dV; W[m]+=delta*delta*dV; K2[m]+=kval[m]*kval[m]*delta*delta*dV;
        profiles<<label<<','<<modNames[m]<<','<<r<<','<<qpt.point.thetaRadians<<','<<qpt.point.phiRadians<<','<<ch<<','<<wr<<','<<dV<<','<<delta<<','<<kval[m]<<'\n'; }
    }
    for(int m=0;m<5;++m){ const double scale=S/(baselineModuliPa[m]*Vsupport); const double kw=std::sqrt(K2[m]/W[m]); integrals<<label<<','<<modNames[m]<<','<<S<<','<<result.objective<<','<<Vsupport<<','<<scale<<','<<I[m]<<','<<kw<<'\n'; }
    std::array<double,5> sourceLocalVolume{},sourceLocalContact{};
    double sourceLocalQuadratureVolume=0.0;
    for(std::size_t i=0;i<sourceLocalQuadrature.size();++i){
      const auto&qpt=sourceLocalQuadrature[i];
      const double r=qpt.point.radiusMetres,theta=qpt.point.thetaRadians,
                   phi=qpt.point.phiRadians,w=sourceLocalWeight(r),ch=chi(theta,phi);
      const double dV=r*r*qpt.radialWeightMetres*qpt.angularWeightSteradians;
      sourceLocalQuadratureVolume+=dV;
      const auto&kernel=result.pointKernels[points.size()+i];
      const std::array<double,5> kval{kernel.A,kernel.C,kernel.F,kernel.L,kernel.N};
      for(int m:{1,2,3}){
        const double delta=baselineModuliPa[m]*w*ch;
        const double contribution=kval[m]*delta*dV;
        sourceLocalVolume[m]+=contribution;
        sourceLocalProfiles<<label<<','<<modNames[m]<<','<<qpt.point.radialElement<<','<<r<<','
          <<theta<<','<<phi<<','<<ch<<','<<w<<','<<dV<<','<<delta<<','<<kval[m]<<','
          <<contribution<<'\n';
      }
    }
    if(std::abs(sourceLocalQuadratureVolume-sourceLocalSupportVolume)>
       1e-12*sourceLocalSupportVolume)
      throw std::runtime_error("source-local closed-support volume does not match analytic shell volume");
    for(std::size_t i=0;i<angles.size();++i){
      const auto&c=result.sourceContact[i];
      const std::array<double,5>kv{c.A,c.C,c.F,c.L,c.N};
      const auto [theta,phi]=angles[i];
      const double ch=chi(theta,phi),dOm=glw[static_cast<int>(i/nangPhi)]*dOmega;
      for(int m=0;m<5;++m){
        const double delta=baselineModuliPa[m]*ch;
        const double contribution=kv[m]*delta*dOm;
        contacts<<label<<','<<modNames[m]<<','<<theta<<','<<phi<<','<<ch<<','
          <<kv[m]<<','<<delta<<','<<dOm<<','<<contribution<<'\n';
        if(m==1||m==2||m==3)sourceLocalContact[m]+=contribution;
      }
    }
    for(int m:{1,2,3}){
      const double scale=S/(baselineModuliPa[m]*sourceLocalSupportVolume);
      sourceLocalIntegrals<<label<<','<<modNames[m]<<','<<S<<','<<sourceLocalSupportVolume<<','
        <<scale<<','<<sourceLocalVolume[m]<<','<<sourceLocalContact[m]<<','
        <<sourceLocalVolume[m]+sourceLocalContact[m]<<'\n';
    }
  }
  exportOracleBlocks(sem,source,parameters,freq,base.GravitationalConstant(),out);
  const double analyticVsupport=(4.0*std::numbers::pi/3.0)*(r1*r1*r1-r0*r0*r0);
  const double quadratureVsupport=std::accumulate(radialNodes.begin(),radialNodes.end(),0.0,[&](double sum,const auto&v){return sum+v[2]*v[2]*v[3]*4.0*std::numbers::pi;});
  if(std::abs(quadratureVsupport-analyticVsupport)>1e-12*analyticVsupport)
    throw std::runtime_error("element-local regular-support volume does not match analytic shell volume");
  std::ofstream meta(out/"meta.json"); meta<<std::setprecision(17)<<"{\"solver\":\"DSpecM1D\",\"nt\":"<<freq.nt()<<",\"dt_s\":"<<freq.dt()*freq.timeNorm()<<",\"df_hz\":"<<freq.df()/freq.timeNorm()<<",\"epsilon_s_inv\":"<<freq.ep()/freq.timeNorm()<<",\"physical_G\":"<<sem.gravityCoefficient()*base.GravitationalConstant()<<",\"radius_m\":"<<sem.mesh().PR()*sem.lengthScaleMetres()<<",\"source_radius_m\":"<<sem.sourceRadius()*sem.lengthScaleMetres()<<",\"source_element\":"<<sem.sourceElement(source)<<",\"source_element_bounds_m\":["<<sem.mesh().ELR(sem.sourceElement(source))*sem.lengthScaleMetres()<<','<<sem.mesh().EUR(sem.sourceElement(source))*sem.lengthScaleMetres()<<"],\"receiver_geocentric\":[45,90,3],\"source_geocentric\":[10,20,50],\"moment_tensor_Nm\":[1e19,2e18,-3e18,-5e18,4e18,-5e18],\"gravity\":\"full\",\"attenuation\":false,\"degrees\":[0,3],\"solver_bins\":[1,2,3],\"window_samples\":91,\"window\":\"Hann\",\"maxstep_fraction\":"<<maxstep<<",\"nq\":"<<nq<<",\"angular_rule\":\"GL8 x periodic24\",\"regular_support_m\":[4200000,6200000],\"regular_support_volume_m3\":"<<quadratureVsupport<<",\"analytic_regular_support_volume_m3\":"<<analyticVsupport<<",\"source_local_support_m\":[6311000,6331000]}\n";
  return 0;
} catch(const std::exception&e){ std::cerr<<"Phase D DSpec adapter: "<<e.what()<<'\n'; return 1; }
