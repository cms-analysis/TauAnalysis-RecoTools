#include "TauAnalysis/RecoTools/plugins/DistortedTauProducer.h"

DistortedTauProducer::~DistortedTauProducer() 
{
	
}


DistortedTauProducer::DistortedTauProducer(const edm::ParameterSet& iConfig):
      src_(iConfig.getParameter<edm::InputTag>("src")),  
      fileName_(iConfig.getParameter<edm::FileInPath>("fileName")),  
      smearMCParticle_(iConfig.getParameter<bool>("smearMCParticle")),  
      smearConstituents_(iConfig.getParameter<bool>("smearConstituents")),  
      smearFromPtHistogram_(iConfig.getParameter<bool>("smearFromPtHistogram")),  
      smearFromEtaHistogram_(iConfig.getParameter<bool>("smearFromEtaHistogram")),  
      ptHisto_(iConfig.getParameter<std::string>("ptHistogram")),			  
      etaHisto_(iConfig.getParameter<std::string>("etaHistogram")),			  
      energyScale_(iConfig.getParameter<double>("energyScale")),
      hadronEnergyScale_(iConfig.getParameter<double>("hadronEnergyScale")),
      gammaEnergyScale_(iConfig.getParameter<double>("gammaEnergyScale")),
      deltaEta_(iConfig.getParameter<double>("deltaEta")),
      deltaPt_(iConfig.getParameter<double>("deltaPt")),
      gSigmaPt_(iConfig.getParameter<double>("gaussianSmearingSigmaPt")),
      gSigmaEta_(iConfig.getParameter<double>("gaussianSmearingSigmaEta")) 
{

      //Check if the file is there
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

      produces<pat::TauCollection >();
}

void 
DistortedTauProducer::beginJob(const edm::EventSetup& es)
{}


void 
DistortedTauProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) 
{
  //get the collection
  using namespace edm;
  using namespace reco;
  
  std::auto_ptr<pat::TauCollection > out(new pat::TauCollection);
  Handle<pat::TauCollection> srcH;
  
  if(iEvent.getByLabel(src_,srcH) &&srcH->size()>0) 
    for(unsigned int i=0;i<srcH->size();++i) {
      pat::Tau  object = srcH->at(i);
      
      math::XYZTLorentzVector vToSmear;

      
      //Decide which LV to smear
      if(smearMCParticle_&&object.genJet()!=0)
	vToSmear = object.genJet()->p4();
      else vToSmear = object.p4();

      
      //apply energy scale!
      math::PtEtaPhiMLorentzVector vecES(vToSmear.pt()*energyScale_,vToSmear.eta(),vToSmear.phi(),vToSmear.M());
      math::XYZTLorentzVector cartVecES(vecES.px(),vecES.py(),vecES.pz(),vecES.energy());
      object.setP4(cartVecES);

      
      //If smearing reconstructed particle apply the eScales of hadrons and gammas seperately
      if(smearConstituents_) {
	
	math::XYZTLorentzVector hadronLV;
	PFCandidateRefVector hadrons = object.signalPFChargedHadrCands();
	if(hadrons.size()>0)
	  for(unsigned int i=0;i<hadrons.size();++i)
	    hadronLV+=hadrons.at(i)->p4();
	      //apply hadron energy scale

	math::PtEtaPhiMLorentzVector vechES(hadronLV.pt()*hadronEnergyScale_,hadronLV.eta(),hadronLV.phi(),hadronLV.M());
	math::XYZTLorentzVector cartVechES(vechES.px(),vechES.py(),vechES.pz(),vechES.energy());
	hadronLV=cartVechES;
	
	math::XYZTLorentzVector gammaLV;
	PFCandidateRefVector gammas = object.signalPFGammaCands();
	if(gammas.size()>0)
	  for(unsigned int i=0;i<gammas.size();++i)
	    gammaLV+=gammas.at(i)->p4();


	math::PtEtaPhiMLorentzVector vecgES(gammaLV.pt()*gammaEnergyScale_,gammaLV.eta(),gammaLV.phi(),gammaLV.M());
	math::XYZTLorentzVector cartVecgES(vecgES.px(),vecgES.py(),vecgES.pz(),vecgES.energy());
	gammaLV=cartVecgES;

	object.setP4(gammaLV+hadronLV);
      }

      
      //apply pt and eta
      //In eta since we need it for acceptabnce reasons use a sign convension
      
      double deta=0;
      if(object.eta()>0)
	deta = deltaEta_;
      else
	deta = -deltaEta_;
      math::PtEtaPhiMLorentzVector vec(object.pt()+deltaPt_,object.eta()+deta,object.phi(),object.mass());
      math::XYZTLorentzVector cartVec(vec.px(),vec.py(),vec.pz(),vec.energy());
      object.setP4(cartVec);
      
      //Apply gaussian smearing
      if ( ! rng.isAvailable()) {
	throw cms::Exception("Configuration")
	  << "DistortedTauProducer requires the RandomNumberGeneratorService\n"
		  "which is not present in the configuration file.  You must add the service\n"
	  "in the configuration file or remove the modules that require it.";
      }
      
      CLHEP::HepRandomEngine& engine = rng->getEngine();
      
      CLHEP::RandFlat* flatDistribution = new CLHEP::RandFlat(engine, 0., 1.);
      CLHEP::RandGauss* gaussPt = new CLHEP::RandGauss(engine, 0.,gSigmaPt_ );
      CLHEP::RandGauss* gaussEta = new CLHEP::RandGauss(engine, 0.,gSigmaEta_ );
      
      double deltaPt = gaussPt->fire();
      double deltaEta = gaussEta->fire();
      math::PtEtaPhiMLorentzVector vecGauss(object.pt()+deltaPt,object.eta()+deltaEta,object.phi(),object.mass());
      math::XYZTLorentzVector cartVecGauss(vecGauss.px(),vecGauss.py(),vecGauss.pz(),vecGauss.energy());
      object.setP4(cartVecGauss);

      
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
      
      delete flatDistribution;
      delete gaussPt;
      delete gaussEta;

      printf("Final Object Pt %f\n",object.pt());
      out->push_back(object);
    }
  

  iEvent.put(out);
}


void
DistortedTauProducer::endJob()
{}

	
TH1F* 
DistortedTauProducer::integratedHisto(TH1F* histo) {
  histo->Scale(1/histo->Integral());
  TH1F * intHisto = new TH1F("intHisto","",histo->GetNbinsX(),histo->GetXaxis()->GetXmin(),histo->GetXaxis()->GetXmax());
  float integral = 0.0;
  for(int i=1;i<=histo->GetNbinsX();++i) {
    integral+=histo->GetBinContent(i);
    intHisto->SetBinContent(i,integral);
  }
  return intHisto;
}

double 
DistortedTauProducer::sampleFromHisto(TH1F * histo, double x)
{
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

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_ANOTHER_FWK_MODULE(DistortedTauProducer);
