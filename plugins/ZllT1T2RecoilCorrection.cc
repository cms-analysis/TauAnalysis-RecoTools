#include "TauAnalysis/RecoTools/plugins/ZllT1T2RecoilCorrection.h"

#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/Common/interface/RefToPtr.h"

#include "TauAnalysis/CandidateTools/interface/candidateAuxFunctions.h"

#include <TMath.h>

const std::string instNameMEtObjects = "met";
const std::string instNameDiTauToMEtAssociations = "diTauToMEtAssociations";

const double sigmaU1min = 5.0;
const double sigmaU2min = 5.0;

template <typename T1, typename T2>
ZllT1T2RecoilCorrection<T1,T2>::ZllT1T2RecoilCorrection(const edm::ParameterSet& cfg)
  : corrParameterData_(0),
    corrParameterMC_(0)
{
  src_ = cfg.getParameter<edm::InputTag>("src");

  edm::ParameterSet cfgCorrParameter = cfg.getParameter<edm::ParameterSet>("parameter"); 
  corrParameterData_ = new corrParameterType(cfgCorrParameter.getParameter<edm::ParameterSet>("data"));
  corrParameterMC_ = new corrParameterType(cfgCorrParameter.getParameter<edm::ParameterSet>("mc"));

  //std::cout << "corr. Parameter(Data):" << std::endl;
  //corrParameterData_->print();
  //std::cout << "corr. Parameter(MC):" << std::endl;
  //corrParameterMC_->print();

  produces<pat::METCollection>(instNameMEtObjects);
  produces<diTauToMEtAssociation>(instNameDiTauToMEtAssociations);
}

template <typename T1, typename T2>
ZllT1T2RecoilCorrection<T1,T2>::~ZllT1T2RecoilCorrection()
{
  delete corrParameterData_;
  delete corrParameterMC_;
}

template <typename T1, typename T2>
std::pair<double, double> ZllT1T2RecoilCorrection<T1,T2>::compSigma(const corrParameterType& corrParameter, double qT)
{
  double sigmaU1 = corrParameter.sigma1_*(1. + corrParameter.b1_*qT + corrParameter.c1_*qT*qT);
  double sigmaU2 = corrParameter.sigma2_*(1. + corrParameter.b2_*qT + corrParameter.c2_*qT*qT);
  if ( sigmaU1 < sigmaU1min ) sigmaU1 = sigmaU1min;
  if ( sigmaU2 < sigmaU2min ) sigmaU2 = sigmaU2min;
  return std::pair<double, double> (sigmaU1, sigmaU2);
}

template <typename T1, typename T2>
void ZllT1T2RecoilCorrection<T1,T2>::produce(edm::Event& evt, const edm::EventSetup& es)
{
  //std::cout << "<ZllT1T2RecoilCorrection::produce>:" << std::endl;

  std::auto_ptr<pat::METCollection> correctedMETs(new pat::METCollection);

  edm::Handle<CompositePtrCandidateCollection> diTauCollection;
  evt.getByLabel(src_, diTauCollection);
  //std::cout << " num. diTau objects = " << diTauCollection->size() << std::endl;

  // association between diTau candidates and corrected MEt values
  std::auto_ptr<diTauToMEtAssociation> correctedMEtAssociations(new diTauToMEtAssociation(CompositePtrCandidateRefProd(diTauCollection)));

  edm::RefProd<pat::METCollection> correctedMEtRefProd = evt.getRefBeforePut<pat::METCollection>("met");

  size_t numDiTauCandidates = diTauCollection->size();
  for ( size_t iDiTauCandidate = 0; iDiTauCandidate < numDiTauCandidates; ++iDiTauCandidate ) {
    const CompositePtrCandidate& diTauCandidate = diTauCollection->at(iDiTauCandidate);

    //std::cout << "diTauCandidate:" << std::endl;
    //std::cout << " leg1: pt = " << diTauCandidate.p4Leg1gen().pt() << "," 
    //	        << " eta = " << diTauCandidate.p4Leg1gen().eta() << ", phi = " << diTauCandidate.p4Leg1gen().phi() << std::endl;
    //std::cout << " leg2: pt = " << diTauCandidate.p4Leg2gen().pt() << "," 
    //	        << " eta = " << diTauCandidate.p4Leg2gen().eta() << ", phi = " << diTauCandidate.p4Leg2gen().phi() << std::endl;
    //std::cout << "mass(SVfit) **before** Z-recoil correction = " 
    //	        << diTauCandidate.svFitSolution("psKine_MEt_ptBalance")->mass() << std::endl;

    pat::MET correctedMEt(*diTauCandidate.met());

    //std::cout << " MEt: px = " << diTauCandidate.met()->px() << ", py = " << diTauCandidate.met()->py() << std::endl;

    double qX = diTauCandidate.p4gen().px();
    double qY = diTauCandidate.p4gen().py();
    double qT = diTauCandidate.p4gen().pt();

    //std::cout << " qT = " << qT << " (qX = " << qX << ", qY = " << qY << ")" << std::endl;

    if ( qT > 0. ) {

      int errorFlag = 0;
      std::pair<double, double> uTrec = compMEtProjU(diTauCandidate.p4gen(), 
						     diTauCandidate.met()->px(), diTauCandidate.met()->py(), errorFlag);
      std::pair<double, double> uTgen = compMEtProjU(diTauCandidate.p4gen(), 
						     0., 0., errorFlag); // MEt resolution parameters determined 
                                                                         // from Z --> mu+ mu- events;
                                                                         // assume MEt to be zero in **those** events
      double u1rec = uTrec.first;
      double u1gen = uTgen.first;
      double u2rec = uTrec.second;
      double u2gen = uTgen.second;
      
      //std::cout << " gen: u1 = " << u1gen << ", u2 = " << u2gen << std::endl;
      //std::cout << " rec: u1 = " << u1rec << ", u2 = " << u2rec << std::endl;

      std::pair<double, double> uTsigmaData = compSigma(*corrParameterData_, qT);
      double u1sigmaData = uTsigmaData.first;
      double u2sigmaData = uTsigmaData.second;

      //std::cout << " sigma(Data): u1 = " << u1sigmaData << ", u2 = " << u2sigmaData << std::endl;

      std::pair<double, double> uTsigmaMC = compSigma(*corrParameterMC_, qT);
      double u1sigmaMC = uTsigmaMC.first;
      double u2sigmaMC = uTsigmaMC.second;
      
      //std::cout << " sigma(MC): u1 = " << u1sigmaMC << ", u2 = " << u2sigmaMC << std::endl;

      double u1sigmaCorr = ( u1sigmaMC != 0. ) ? (u1sigmaData/u1sigmaMC) : 1.;
      double u2sigmaCorr = ( u2sigmaMC != 0. ) ? (u2sigmaData/u2sigmaMC) : 1.;
            
      double u1rec_corrected = (u1rec - u1gen)*u1sigmaCorr + u1gen;         
      double u2rec_corrected = (u2rec - u2gen)*u2sigmaCorr + u2gen;
      
      double dDiff = corrParameterData_->d_ - corrParameterMC_->d_;
      double kDiff = corrParameterData_->k_ - corrParameterMC_->k_;
      u1rec_corrected += dDiff + kDiff*qT;
      
      //std::cout << "--> rec(corrected): u1 = " << u1rec_corrected << ", u2 = " << u2rec_corrected << std::endl;

      double uXrec_corrected = (u1rec_corrected*qX + u2rec_corrected*qY)/qT;
      double uYrec_corrected = (u1rec_corrected*qY - u2rec_corrected*qX)/qT;

      double correctedMETpx = -(uXrec_corrected + qX);
      double correctedMETpy = -(uYrec_corrected + qY);
      double correctedMETpt = TMath::Sqrt(correctedMETpx*correctedMETpx + correctedMETpy*correctedMETpy);

      //std::cout << "--> corrected MEt: px = " << correctedMETpx << ", py = " << correctedMETpy << std::endl;

      correctedMEt.setP4(math::XYZTLorentzVector(correctedMETpx, correctedMETpy, 0., correctedMETpt));
    }

    correctedMETs->push_back(correctedMEt);

    pat::METRef correctedMEtRef(correctedMEtRefProd, iDiTauCandidate);
    int correctedMEt_index = correctedMEtRef.index();

    correctedMEtAssociations->setValue(iDiTauCandidate, correctedMEt_index);
  }

  //std::cout << "--> num. MEt objects = " << correctedMETs->size() << std::endl;
  
  evt.put(correctedMETs, instNameMEtObjects);
  evt.put(correctedMEtAssociations, instNameDiTauToMEtAssociations);  
}

#include "FWCore/Framework/interface/MakerMacros.h"

typedef ZllT1T2RecoilCorrection<pat::Electron, pat::Tau> ZllRecoilCorrectionElecTauPair;
typedef ZllT1T2RecoilCorrection<pat::Muon, pat::Tau> ZllRecoilCorrectionMuTauPair;
typedef ZllT1T2RecoilCorrection<pat::Tau, pat::Tau> ZllRecoilCorrectionDiTauPair;
typedef ZllT1T2RecoilCorrection<pat::Electron, pat::Muon> ZllRecoilCorrectionElecMuPair;
typedef ZllT1T2RecoilCorrection<pat::Electron, pat::Electron> ZllRecoilCorrectionDiElecPair;
typedef ZllT1T2RecoilCorrection<pat::Muon, pat::Muon> ZllRecoilCorrectionDiMuPair;

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(ZllRecoilCorrectionElecTauPair);
DEFINE_FWK_MODULE(ZllRecoilCorrectionMuTauPair);
DEFINE_FWK_MODULE(ZllRecoilCorrectionDiTauPair);
DEFINE_FWK_MODULE(ZllRecoilCorrectionElecMuPair);
DEFINE_FWK_MODULE(ZllRecoilCorrectionDiElecPair);
DEFINE_FWK_MODULE(ZllRecoilCorrectionDiMuPair);

