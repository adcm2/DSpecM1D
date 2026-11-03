#include "fast/v5/fast_engine.hpp"
#include "fast/v5/engine_headers/SensitivityKernels.h"
#include <DSpecM1D/All>
#include "fast/v5/engine_headers/Windowed3DKernels.h"
#include "test_utils.h"
#include <iomanip>
#include <iostream>
using Complex=std::complex<double>;
void strain(const SPARSESPEC::SpinStrainField &s) {
  for(auto v:{s.rr,s.horizontalTrace,s.radialPlus,s.radialMinus,s.horizontalPlus,s.horizontalMinus})
    std::cout<<','<<v.real()<<','<<v.imag();
}
int main(int argc,char **argv) {
  if(argc!=3)throw std::invalid_argument("contact output and backend required");
  SPARSESPEC::FixedMeshAdjointSolver::setDiagnosticBackend(std::stoi(argv[2]));
  std::ofstream contacts(argv[1]);contacts<<std::setprecision(17);std::cout<<std::setprecision(17);
  prem_norm<double> norm;
  EarthModels::ModelInput model(DSpecMTest::modelPath().string(),norm);
  int caseIndex=0;
  for(auto location:std::array<std::array<double,2>,3>{{{15.,110.},{90.,37.},{-89.99999,-23.}}}) {
    SourceInfo::EarthquakeCMT source(location[0],location[1],30.,1.,2.,3.,4.,5.,6.);
    Full1D::SEM sem(model,.25,4,3,source);const double length=sem.lengthScaleMetres();
    std::vector<SPARSESPEC::Field3DPoint> points;
    for(int e=0;e<sem.mesh().NE();++e) {
      double r=.5*(sem.mesh().ELR(e)+sem.mesh().EUR(e))*length;
      for(auto a:std::array<std::array<double,2>,3>{{{.63,1.17},{0.,0.},{EIGEN_PI,2.1}}})points.push_back({r,a[0],a[1],e});
      if(e>0){double r0=sem.mesh().ELR(e)*length;points.push_back({r0,.9,1.7,e-1});points.push_back({r0,.9,1.7,e});}
    }
    const int sourceElement=sem.sourceElement(source);const double rs=sem.sourceRadius()*length;
    const double thetaS=(90.-location[0])*EIGEN_PI/180.,phiS=location[1]*EIGEN_PI/180.;
    std::vector<std::array<double,2>> angles{{thetaS,phiS},{EIGEN_PI-thetaS,phiS+EIGEN_PI},{0.,.3},{EIGEN_PI,1.2},{EIGEN_PI/3,65.*EIGEN_PI/180.},{.73,1.17},{std::nextafter(.73,1.),1.17}};
    for(auto a:angles){points.push_back({rs,a[0],a[1],sourceElement});points.push_back({rs,a[0],a[1],sourceElement-1});}
    for(int i=0;i<70;++i)points.push_back({rs,.2+2.6*i/69.,.17+.11*i,sourceElement});
    points.push_back({0.,.7,1.1,0});points.push_back(points.front());std::reverse(points.begin(),points.end());
    for(int receiverSet=0;receiverSet<2;++receiverSet){
      DSpecMTest::TempDir temp;DSpecMTest::ParameterOptions options;options.lmax=3;options.receiverDepth=3.;
      options.receivers=receiverSet?std::vector<std::pair<double,double>>{{30.,65.},{90.,17.},{-90.,123.}}:std::vector<std::pair<double,double>>{{30.,65.}};
      options.numReceivers=options.receivers.size();
      auto path=DSpecMTest::writeFile(temp.path()/"parameters.txt",DSpecMTest::makeParameterText(DSpecMTest::modelPath().string(),options));InputParameters parameters(path.string());
      for(int weights=0;weights<4;++weights,++caseIndex){
        Eigen::VectorXcd q=Eigen::VectorXcd::Zero(3*options.numReceivers);
        for(int r=0;r<options.numReceivers;++r){
          double by=1.+.2*r;
          if(weights==1||weights==3)q[3*r]=by*Complex{.7,-.2};
          if(weights==2||weights==3){q[3*r+1]=by*Complex{-.3,.5};q[3*r+2]=by*Complex{.2,.4};}
        }
        int frequency=0;
        for(Complex z:{Complex{.31,-.02},Complex{.47,-.03}}){
          auto workspace=SPARSESPEC::Windowed3DFieldEngine::prepareFrequencyWorkspace(sem,source,parameters,z,q);
          const auto fields=workspace.sampleWeightedChunk(points);
          std::vector<SPARSESPEC::Field3DPoint> noncenter;
          std::vector<std::size_t> original;
          for(std::size_t i=0;i<points.size();++i)if(points[i].radiusMetres>0){noncenter.push_back(points[i]);original.push_back(i);}
          std::map<std::pair<int,double>,SPARSESPEC::Field3DPoint> tracePoints;
          std::map<std::pair<double,double>,SPARSESPEC::Field3DPoint> anglePoints;
          for(const auto &p:noncenter){tracePoints.try_emplace(std::pair{p.radialElement,p.radiusMetres},p);anglePoints.try_emplace(std::pair{p.thetaRadians,p.phiRadians},p);}
          noncenter.clear();
          for(const auto &[tk,tp]:tracePoints)for(const auto &[ak,ap]:anglePoints){(void)tk;(void)ak;noncenter.push_back({tp.radiusMetres,ap.thetaRadians,ap.phiRadians,tp.radialElement});}
          std::reverse(noncenter.begin(),noncenter.end()); // Exercise a permuted complete tensor.
          const auto tensorFull=workspace.sampleWeightedChunk(noncenter);
          SPARSESPEC::Windowed3DFieldEngine::SamplingPlanCache packedPlan;
          if(std::getenv("FORCE_PLAN_FULL"))packedPlan.bytes=SPARSESPEC::Windowed3DFieldEngine::SamplingPlanCache::byteBudget;
          const bool realTensor=!(std::getenv("TENSOR_COMPLEX"));
          const auto packed=workspace.sampleForKernelsTensor(noncenter,noncenter.size(),&packedPlan,3,{},realTensor);
          std::array<double,3> numerator{},denominator{};
          for(std::size_t j=0;j<packed.size();++j){
            const auto &a=packed[j];const auto &b=tensorFull[j];
            for(int c=0;c<3;++c){numerator[0]+=0.0;denominator[0]+=0.0;}
            const auto av=a.forwardRawStrain,bv=b.forwardRawStrain,aa=a.adjointDualStrain,ba=b.adjointDualStrain;
            if(a.fluid){numerator[1]+=std::norm(av.rr+av.horizontalTrace-bv.rr-bv.horizontalTrace);denominator[1]+=std::norm(bv.rr+bv.horizontalTrace);
              numerator[2]+=std::norm(aa.rr+aa.horizontalTrace-ba.rr-ba.horizontalTrace);denominator[2]+=std::norm(ba.rr+ba.horizontalTrace);}
            else{const std::array<Complex,6> af{av.rr,av.horizontalTrace,av.radialPlus,av.radialMinus,av.horizontalPlus,av.horizontalMinus},bf{bv.rr,bv.horizontalTrace,bv.radialPlus,bv.radialMinus,bv.horizontalPlus,bv.horizontalMinus}, ax{aa.rr,aa.horizontalTrace,aa.radialPlus,aa.radialMinus,aa.horizontalPlus,aa.horizontalMinus},bx{ba.rr,ba.horizontalTrace,ba.radialPlus,ba.radialMinus,ba.horizontalPlus,ba.horizontalMinus};
              for(int c=0;c<6;++c){numerator[1]+=std::norm(af[c]-bf[c]);denominator[1]+=std::norm(bf[c]);numerator[2]+=std::norm(ax[c]-bx[c]);denominator[2]+=std::norm(bx[c]);}}
          }
          std::cerr<<std::setprecision(17)<<"PACKERR {\"case\":"<<caseIndex<<",\"frequency\":"<<frequency<<",\"relative_norms\":[";
          for(int c=0;c<3;++c)std::cerr<<(c?",":"")<<(denominator[c]==0?(numerator[c]==0?0:std::numeric_limits<double>::infinity()):std::sqrt(numerator[c]/denominator[c]));
          std::cerr<<"]}\n";

          for(std::size_t i=0;i<fields.size();++i){const auto &s=fields[i];std::cout<<caseIndex<<','<<frequency<<','<<i<<','<<s.radialElement<<','<<s.fluid;
            for(const auto &a:{s.forwardRawDisplacement,s.observableForwardDisplacementSpectrumMetreSeconds,s.adjointDualDisplacement})for(auto v:a)std::cout<<','<<v.real()<<','<<v.imag();
            strain(s.forwardRawStrain);strain(s.adjointDualStrain);strain(s.rScaledForwardStrainMetres);strain(s.rScaledAdjointDualStrainMetres);std::cout<<'\n';
          }
          const auto density=SPARSESPEC::Windowed3DKernelEngine::sourceShellContact(workspace,source,angles);
          for(std::size_t i=0;i<density.size();++i){const auto &v=density[i];contacts<<caseIndex<<','<<frequency<<','<<i<<','<<v.A<<','<<v.C<<','<<v.F<<','<<v.L<<','<<v.N<<'\n';}
          ++frequency;
        }
      }
    }
  }
}
