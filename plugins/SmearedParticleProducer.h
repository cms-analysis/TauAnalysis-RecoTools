#ifndef TauAnalysis_SmearedParticleProducer
#define TauAnalysis_SmearedParticleProducer
/* Smeared Particle Producer 
Module that changes the energy scale PT, ETA 
for a general PAT Particle 
Author : Michail Bachtis 
University of Wisconsin
*/
// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"
#include <TFile.h>
#include <TH1F.h>
/*Helper classes to return the generated Lorentz Vector in case of leptons and jets*/

template <typename T>
class GenParticleRetriever {
 public:
  GenParticleRetriever() {}
  ~GenParticleRetriever() {}
  GenParticleRetriever(const T& particle) {
    if(particle.genParticle()!=0)
      p4_ = particle.genParticle()->p4();
    else
      p4_ = particle.p4();
  }

  math::XYZTLorentzVector genP4() { return p4_; }
  
 private:
  math::XYZTLorentzVector p4_;

};

template <typename T>
class GenJetRetriever {
 public:
  GenJetRetriever() {}
  ~GenJetRetriever() {}
  GenJetRetriever(const T& particle) {
    if(particle.genJet()!=0)
      p4_ = particle.genJet()->p4();
    else
      p4_ = particle.p4();
  }
  math::XYZTLorentzVector genP4() { return p4_; }
 private:
  math::XYZTLorentzVector p4_;
};


/*end helper classes*/


template <typename T,typename G>
class SmearedParticleProducer : public edm::EDProducer {
   public:
  explicit SmearedParticleProducer(const edm::ParameterSet& iConfig):
    src_(iConfig.getParameter<edm::InputTag>("src")),  
    fileName_(iConfig.getParameter<edm::FileInPath>("fileName")),  
    smearMCParticle_(iConfig.getParameter<bool>("smearMCParticle")),  
    smearFromPtHistogram_(iConfig.getParameter<bool>("smearFromPtHistogram")),  
    smearFromEtaHistogram_(iConfig.getParameter<bool>("smearFromEtaHistogram")),  
    ptHisto_(iConfig.getParameter<std::string>("ptHistogram")),			  
    etaHisto_(iConfig.getParameter<std::string>("etaHistogram")),			  
    energyScale_(iConfig.getParameter<double>("energyScale")),
    deltaEta_(iConfig.getParameter<double>("deltaEta")),
    deltaPhi_(iConfig.getParameter<double>("deltaPhi")),
    deltaPt_(iConfig.getParameter<double>("deltaPt")),
    gSigmaPt_(iConfig.getParameter<double>("gaussianSmearingSigmaPt")),
    gSigmaEta_(iConfig.getParameter<double>("gaussianSmearingSigmaEta")), 
    gSigmaPhi_(iConfig.getParameter<double>("gaussianSmearingSigmaPhi")) ,
    gSigmaEScale_(iConfig.getParameter<double>("gaussianSmearingSigmaEScale")) 
    {
      //Check if the file with the histograms is there
      if(smearFromPtHistogram_||smearFromEtaHistogram_) {
	if(fileName_.isLocal()) {
	  f = new TFile(fileName_.fullPath().c_str());
	  if(smearFromPtHistogram_) {
	    TH1F *tmp = (TH1F*)f->Get(ptHisto_.c_str());
	    ptHisto = integratedHisto(tmp);
	  }
	  if(smearFromEtaHistogram_) {
	    TH1F *tmp = (TH1F*)f->Get(etaHisto_.c_str());
	    etaHisto = integratedHisto(tmp);
	  }
	}
	else {
	  throw cms::Exception("FileNotFound") << "File Not found\n";
	}
	
      }

      //Initalize random number service
	if ( ! rng.isAvailable()) {
	  throw cms::Exception("Configuration")
	    << "DistortedTauProducer requires the RandomNumberGeneratorService\n"
	    "which is not present in the configuration file.  You must add the service\n"
	    "in the configuration file or remove the modules that require it.";
	}

      CLHEP::HepRandomEngine& randomEngine = rng->getEngine();
      flatDistribution = new CLHEP::RandFlat(randomEngine, 0., 1.);
      gaussPt = new CLHEP::RandGauss(randomEngine, 0.,gSigmaPt_ );
      gaussEta = new CLHEP::RandGauss(randomEngine, 0.,gSigmaEta_ );
      gaussPhi = new CLHEP::RandGauss(randomEngine, 0.,gSigmaPhi_ );
      gaussEScale = new CLHEP::RandGauss(randomEngine, 1.,gSigmaEScale_ );

      produces<std::vector<T> >();
    }

  ~SmearedParticleProducer() {
    if(smearFromPtHistogram_||smearFromEtaHistogram_) 
	f->Close();

    if(ptHisto !=0) 
      delete ptHisto;
    if(etaHisto !=0) 
      delete etaHisto;
    if(flatDistribution !=0)
      delete flatDistribution;
    if(gaussPt !=0)
      delete gaussPt;
    if(gaussEta !=0)
      delete gaussEta;
    if(gaussPhi !=0)
      delete gaussPhi;
    if(gaussEScale !=0)
      delete gaussEScale;
    if(etaHisto !=0)
      delete gaussEta;
  }
    
   protected:
     void smear(T& object)
       {


	 math::XYZTLorentzVector vToSmear;
	 
	 //Decide which LV to smear
	 if(smearMCParticle_){
	   //check if it matched to a gen jet or a gen particle
	   G mcRetriever(object);
	   vToSmear = mcRetriever.genP4(); //If it is not matched it will return the reco p4
	 }
	 else {
	   vToSmear = object.p4();
	 }
      
	 //apply energy scale!
	 //following Christian's suggestion and applying energy scale in the whole Lorentz Vector
	 math::XYZTLorentzVector cartVecES = energyScale_*vToSmear;
	 object.setP4(cartVecES);
	 
	 //apply pt and eta and phi DISPLACEMENTS
	 math::PtEtaPhiMLorentzVector vec(object.pt()+deltaPt_,object.eta()+deltaEta_,object.phi()+deltaPhi_,object.mass());
	 math::XYZTLorentzVector cartVec(vec.px(),vec.py(),vec.pz(),vec.energy());
	 object.setP4(cartVec);
	 
	 //Apply gaussian smearing
	 double deltaPt = gaussPt->fire();
	 double deltaEta = gaussEta->fire();
	 double deltaPhi = gaussPhi->fire();
	 double deltaEScale = gaussEScale->fire();

	 math::PtEtaPhiMLorentzVector vecGauss(object.pt()+deltaPt,object.eta()+deltaEta,object.phi()+deltaPhi,object.mass());
	 math::XYZTLorentzVector cartVecGauss(vecGauss.px(),vecGauss.py(),vecGauss.pz(),vecGauss.energy());
	 object.setP4(deltaEScale*cartVecGauss);
	
	 //OK Now histogram Smearing
	 double hDeltaPt = 0.0;
	 double hDeltaEta = 0.0;
	 
	 if(smearFromPtHistogram_) {
	   double randomNumber = flatDistribution->fire();
	   hDeltaPt = sampleFromHisto(ptHisto,randomNumber);
	 }
	 if(smearFromEtaHistogram_) {
	   double randomNumber2 = flatDistribution->fire();
	  hDeltaEta = sampleFromHisto(etaHisto,randomNumber2);
	 }
	 math::PtEtaPhiMLorentzVector vecH(object.pt()+hDeltaPt,object.eta()+hDeltaEta,object.phi(),object.mass());
	 math::XYZTLorentzVector cartVecH(vecH.px(),vecH.py(),vecH.pz(),vecH.energy());
	 object.setP4(cartVecH);
       }
     
     virtual void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
       using namespace edm;
       
       std::auto_ptr<std::vector<T> > out(new std::vector<T> );
       Handle<std::vector<T> > srcH;
       
       if(iEvent.getByLabel(src_,srcH)) 
	 for(unsigned int i=0;i<srcH->size();++i) {
	   T  object = srcH->at(i);
	   smear(object);
	   out->push_back(object);
	 }
       
       iEvent.put(out);
       
     } 


     TH1F* integratedHisto(TH1F* histo) {
       histo->Scale(1/histo->Integral());
       TH1F * intHisto = new TH1F("intHisto","",histo->GetNbinsX(),histo->GetXaxis()->GetXmin(),histo->GetXaxis()->GetXmax());
       float integral = 0.0;
       for(int i=1;i<=histo->GetNbinsX();++i) {
	 integral+=histo->GetBinContent(i);
	 intHisto->SetBinContent(i,integral);
       }
       return intHisto;
     }


     double sampleFromHisto(TH1F * histo, double x) {
       double min = histo->GetMinimum();
       double max = histo->GetMaximum();
       int  bins  = histo->GetNbinsX();
       int bin=0;
       for(int i=1;i<=bins;++i) {
	 
	 double low;
	 double high;
	 
	 if(i==1) {
	   low=min;
	   high =  histo->GetBinContent(2);
	 }
	 else if (i==bins) {
	   low= histo->GetBinContent(bins-1);
	   high= max;
	 }
	 else {
	   low  =  histo->GetBinContent(i-1);
	   high =  histo->GetBinContent(i+1);
	 }
	 
	 if(x>low&&x<=high)
	   bin=i;
	 
       }
    
    return histo->GetBinCenter(bin);
    
}

  
      // ----------member data ---------------------------
      edm::InputTag src_;           //input Collection
      edm::FileInPath fileName_;    //Name of the file with the histograms
      bool smearMCParticle_;        //Smear MC particle instead of reconstructed particle
      bool smearFromPtHistogram_;   //Smear from an ET  resolution histogram
      bool smearFromEtaHistogram_;  //Smear from an Eta resolution histogram 
      std::string ptHisto_;         //Histogram name in File
      std::string etaHisto_;        //Histogram name in File

      //Flat smearing
      double energyScale_;          //energy Scale 
      double deltaEta_;             //delta Eta
      double deltaPt_;              //delta Pt
      double deltaPhi_;              //delta phi

      //Random number generator services
      edm::Service<edm::RandomNumberGenerator> rng;
      CLHEP::RandFlat* flatDistribution; 
      CLHEP::RandGauss* gaussPt;  
      CLHEP::RandGauss* gaussEta; 
      CLHEP::RandGauss* gaussPhi; 
      CLHEP::RandGauss* gaussEScale; 

      //Gaussian Smearing
      double gSigmaPt_;
      double gSigmaEta_;
      double gSigmaPhi_;
      double gSigmaEScale_;


      //Private data
      TFile *f;
      TH1F *ptHisto;
      TH1F *etaHisto;
};


#endif
