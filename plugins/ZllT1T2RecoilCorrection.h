#ifndef TauAnalysis_RecoTools_ZllT1T2RecoilCorrection_h
#define TauAnalysis_RecoTools_ZllT1T2RecoilCorrection_h

/** \class ZllT1T2RecoilCorrection
 *
 * Correct for difference in hadronic recoil between Monte Carlo simulation and Data
 * in Z --> l+ l- events
 * (cf. CMS AN-10-264 for a description of the method)
 * 
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.4 $
 *
 * $Id: ZllT1T2RecoilCorrection.h,v 1.4 2010/12/12 08:53:09 veelken Exp $
 *
 */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "TauAnalysis/CandidateTools/interface/CompositePtrCandidateT1T2MEtAlgorithm.h"

#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEt.h"
#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEtFwd.h"

#include "DataFormats/PatCandidates/interface/MET.h"

#include "DataFormats/Common/interface/AssociationVector.h"
#include "DataFormats/Common/interface/RefProd.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/Ptr.h"

template <typename T1, typename T2>
class ZllT1T2RecoilCorrection : public edm::EDProducer 
{
  typedef CompositePtrCandidateT1T2MEt<T1,T2> CompositePtrCandidate;
  typedef std::vector<CompositePtrCandidate> CompositePtrCandidateCollection;
  typedef edm::RefProd<CompositePtrCandidateCollection> CompositePtrCandidateRefProd;
  typedef std::vector<int> vint;
  typedef edm::AssociationVector<CompositePtrCandidateRefProd, vint> diTauToMEtAssociation;

 public:
  explicit ZllT1T2RecoilCorrection(const edm::ParameterSet&);
  ~ZllT1T2RecoilCorrection();

  void produce(edm::Event&, const edm::EventSetup&);

 private:

//--- configuration parameters

  std::string moduleLabel_;

  // collection of diTau objects representing Z --> l+ l- candidates
  edm::InputTag src_;

  // collection of gen. particles for determining Z/W/H Pt
  // used in parametrization of recoil correction
  edm::InputTag srcGenParticles_;
  vint genParticlePdgIds_;

  edm::InputTag srcGenMET_;

  // MEt resolution parameters for MC/Data
  struct corrParameterType
  {
    corrParameterType(const edm::ParameterSet& cfg)
      : d_(cfg.getParameter<double>("d")),
        k_(cfg.getParameter<double>("k")),
	sigma1_(cfg.getParameter<double>("sigma1")),
        b1_(cfg.getParameter<double>("b1")),
	c1_(cfg.getParameter<double>("c1")),
        sigma2_(cfg.getParameter<double>("sigma2")),
	b2_(cfg.getParameter<double>("b2")),
        c2_(cfg.getParameter<double>("c2"))
    {}
    ~corrParameterType() {}
    void print()
    {
      std::cout << " d = " << d_ << std::endl;
      std::cout << " k = " << d_ << std::endl;
      std::cout << " sigma1 = " << sigma1_ << std::endl;
      std::cout << " b1 = " << b1_ << std::endl;
      std::cout << " c1 = " << c1_ << std::endl;
      std::cout << " sigma2 = " << sigma2_ << std::endl;
      std::cout << " b2 = " << b2_ << std::endl;
      std::cout << " c2 = " << c2_ << std::endl;
    }
    double d_;
    double k_;
    double sigma1_;
    double b1_;
    double c1_;
    double sigma2_;
    double b2_;
    double c2_;
  };

  std::pair<double, double> compSigma(const corrParameterType&, double);

  corrParameterType* corrParameterData_;
  corrParameterType* corrParameterMC_;

  corrParameterType* corrUncertaintyData_;
  corrParameterType* corrUncertaintyMC_;

  double shiftByUncertainty_;
};

#endif


