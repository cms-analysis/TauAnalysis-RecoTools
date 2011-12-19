#include "TauAnalysis/RecoTools/interface/ZllRecoilCorrectionHistManager.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "TauAnalysis/CandidateTools/interface/candidateAuxFunctions.h"

#include <TMath.h>

ZllRecoilCorrectionHistManager::ZllRecoilCorrectionHistManager(const edm::ParameterSet& cfg)
{}

ZllRecoilCorrectionHistManager::~ZllRecoilCorrectionHistManager()
{
  for ( std::vector<histogramsUvsQtNumObjType*>::iterator it = histogramsUvsQtNumVtxBinned_.begin();
	it != histogramsUvsQtNumVtxBinned_.end(); ++it ) {
    delete (*it);
  }

  for ( std::vector<histogramsUvsQtNumObjType*>::iterator it = histogramsUvsQtNumJetsBinned_.begin();
	it != histogramsUvsQtNumJetsBinned_.end(); ++it ) {
    delete (*it);
  }
}

void ZllRecoilCorrectionHistManager::bookHistograms(TFileDirectory& dir)
{
  histogramLepPlusPt_          = book1D(dir, "lPlusPt",            "P_{T}^{l+}",                              40,          0. ,         100.);
  histogramLepPlusEta_         = book1D(dir, "lPlusEta",           "#eta_{l+}",                               50,         -2.5,         +2.5);
  histogramLepPlusPhi_         = book1D(dir, "lPlusPhi",           "#phi_{l+}",                               36, -TMath::Pi(), +TMath::Pi());

  histogramLepMinusPt_         = book1D(dir, "lMinusPt",           "P_{T}^{l-}",                              40,          0. ,         100.);
  histogramLepMinusEta_        = book1D(dir, "lMinusEta",          "#eta_{l-}",                               50,         -2.5,         +2.5);
  histogramLepMinusPhi_        = book1D(dir, "lMinusPhi",          "#phi_{l-}",                               36, -TMath::Pi(), +TMath::Pi());

  histogramZllCandPt_          = book1D(dir, "ZllCandPt",          "P_{T}^{Z}",                               40,          0. ,         100.);
  histogramZllCandEta_         = book1D(dir, "ZllCandEta",         "#eta_{Z}",                                50,         -2.5,         +2.5);
  histogramZllCandPhi_         = book1D(dir, "ZllCandPhi",         "#phi_{Z}",                                36, -TMath::Pi(), +TMath::Pi());
  histogramZllCandMass_        = book1D(dir, "ZllCandMass",        "M(l+ l-)",                                60,         60. ,         120.);
  
  histogramNumJetsRawPtGt10_   = book1D(dir, "numJetsRawPtGt10",   "Num. Jets (P_{T}^{raw} > 10 GeV)",        50,         -0.5,         49.5);
  histogramNumJetsCorrPtGt10_  = book1D(dir, "numJetsCorrPtGt10",  "Num. Jets (P_{T}^{corr} > 10 GeV)",       20,         -0.5,         19.5);
  histogramNumJetsRawPtGt12_   = book1D(dir, "numJetsRawPtGt12",   "Num. Jets (P_{T}^{raw} > 12 GeV)",        50,         -0.5,         49.5);
  histogramNumJetsCorrPtGt12_  = book1D(dir, "numJetsCorrPtGt12",  "Num. Jets (P_{T}^{corr} > 12 GeV)",       20,         -0.5,         19.5);
  histogramNumJetsRawPtGt15_   = book1D(dir, "numJetsRawPtGt15",   "Num. Jets (P_{T}^{raw} > 15 GeV)",        50,         -0.5,         49.5);
  histogramNumJetsCorrPtGt15_  = book1D(dir, "numJetsCorrPtGt15",  "Num. Jets (P_{T}^{corr} > 15 GeV)",       20,         -0.5,         19.5);
  histogramNumJetsRawPtGt17_   = book1D(dir, "numJetsRawPtGt17",   "Num. Jets (P_{T}^{raw} > 17 GeV)",        50,         -0.5,         49.5);
  histogramNumJetsCorrPtGt17_  = book1D(dir, "numJetsCorrPtGt17",  "Num. Jets (P_{T}^{corr} > 17 GeV)",       20,         -0.5,         19.5);
  histogramNumJetsRawPtGt20_   = book1D(dir, "numJetsRawPtGt20",   "Num. Jets (P_{T}^{raw} > 20 GeV)",        50,         -0.5,         49.5);
  histogramNumJetsCorrPtGt20_  = book1D(dir, "numJetsCorrPtGt20",  "Num. Jets (P_{T}^{corr} > 20 GeV)",       20,         -0.5,         19.5);

  histogramJetPtAbsEtaLt11_    = book1D(dir, "jetPtAbsEtaLt11",    "P_{T}^{jet} (|#eta_{jet}| < 1.1)",       100,          0. ,         250.);
  histogramJetResAbsEtaLt11_   = book1D(dir, "jetResAbsEtaLt11",   "Jet res. (|#eta_{jet}| < 1.1)",           40,         -1. ,         +1.);
  histogramJetPtAbsEta11to17_  = book1D(dir, "jetPtAbsEta11to17",  "P_{T}^{jet} (1.1 < |#eta_{jet}| < 1.7)", 100,          0. ,         250.);
  histogramJetResAbsEta11to17_ = book1D(dir, "jetResAbsEta11to17", "Jet res. (1.1 < |#eta_{jet}| < 1.7)",     40,         -1. ,         +1.);
  histogramJetPtAbsEta17to23_  = book1D(dir, "jetPtAbsEta17to23",  "P_{T}^{jet} (1.7 < |#eta_{jet}| < 2.3)", 100,          0. ,         250.);
  histogramJetResAbsEta17to23_ = book1D(dir, "jetResAbsEta17to23", "Jet res. (1.7 < |#eta_{jet}| < 2.3)",     40,         -1. ,         +1.);
  histogramJetPtAbsEtaGt23_    = book1D(dir, "jetPtAbsEtaGt23",    "P_{T}^{jet} (|#eta_{jet}| > 2.3)",       100,          0. ,         250.);
  histogramJetResAbsEtaGt23_   = book1D(dir, "jetResAbsEtaGt23",   "Jet res. (|#eta_{jet}| > 2.3)",           40,         -1. ,         +1.);

  histogramMEtS_               = book1D(dir, "metS",               "E_{T}^{miss}",                            30,          0.0,         60.0);
  histogramMEtL_               = book1D(dir, "metL",               "E_{T}^{miss}",                            75,          0.0,        150.0);
  histogramMEtProjParlZ_       = book1D(dir, "metProjParlZ",       "E_{T}^{miss} Proj. parallel Z",           75,        -75.0,        +75.0);
  histogramMEtProjPerpZ_       = book1D(dir, "metProjPerpZ",       "E_{T}^{miss} Proj. perp. Z",              50,        -50.0,        +50.0);

  const int qTnumBins = 34;
  double qTbinning[qTnumBins + 1] = { 
    0., 2.5, 5., 7.5, 10., 12.5, 15., 17.5, 20., 22.5, 25., 27.5, 30., 35., 40., 45., 50., 
    60., 70., 80., 90., 100., 110., 120., 130., 140., 150., 160., 170., 180., 200., 220., 240., 260., 300.
  };
  histogramUparlDivQtVsQt_ = book2D(dir, "uParlDivQtVsQt", "u_{#parallel}/q_{T} vs q_{T}",                           
				    qTnumBins, qTbinning, 400,  -5.0,   +5.0);
  histogramUparlVsQt_      = book2D(dir, "uParlVsQt",      "u_{#parallel} vs q_{T}",                           
				    qTnumBins, qTbinning, 230, -500.0,  +75.0);
  histogramUperpDivQtVsQt_ = book2D(dir, "uPerpDivQtVsQt", "u_{#perp}/q_{T} vs q_{T}",                           
				    qTnumBins, qTbinning, 400,  -5.0,   +5.0);
  histogramUperpVsQt_      = book2D(dir, "uPerpVsQt",      "u_{#perp} vs q_{T}",                           
				    qTnumBins, qTbinning,  60, -75.0,  +75.0);
  histogramQt_             = book1D(dir, "qT",             "q_{T}",  
				    600, 0., 300.);
  
  histogramsUvsQtNumVtxBinned_.push_back(
    new histogramsUvsQtNumObjType(this, dir, qTnumBins, qTbinning, "NumVertices", -1,  2));
  histogramsUvsQtNumVtxBinned_.push_back(
    new histogramsUvsQtNumObjType(this, dir, qTnumBins, qTbinning, "NumVertices",  3,  5));
  histogramsUvsQtNumVtxBinned_.push_back(
    new histogramsUvsQtNumObjType(this, dir, qTnumBins, qTbinning, "NumVertices",  6,  8));
  histogramsUvsQtNumVtxBinned_.push_back(
    new histogramsUvsQtNumObjType(this, dir, qTnumBins, qTbinning, "NumVertices",  9, 11));
  histogramsUvsQtNumVtxBinned_.push_back(
    new histogramsUvsQtNumObjType(this, dir, qTnumBins, qTbinning, "NumVertices", 12, -1));

  for ( int iNumVtx = 1; iNumVtx <= 25; ++iNumVtx ) {
    histogramsUvsQtNumVtxBinned_.push_back(
      new histogramsUvsQtNumObjType(this, dir, qTnumBins, qTbinning, "NumVertices", iNumVtx, iNumVtx));
  }

  histogramsUvsQtNumJetsBinned_.push_back(
    new histogramsUvsQtNumObjType(this, dir, qTnumBins, qTbinning, "NumJets",  0,  0));
  histogramsUvsQtNumJetsBinned_.push_back(
    new histogramsUvsQtNumObjType(this, dir, qTnumBins, qTbinning, "NumJets",  1,  1));
  histogramsUvsQtNumJetsBinned_.push_back(
    new histogramsUvsQtNumObjType(this, dir, qTnumBins, qTbinning, "NumJets",  2,  2));
  histogramsUvsQtNumJetsBinned_.push_back(
    new histogramsUvsQtNumObjType(this, dir, qTnumBins, qTbinning, "NumJets",  3, -1));
  
  histogramVtxMultiplicity_  = book1D(dir, "numVertices",        "Num. Vertices",                             25,         -0.5,         24.5);
  histogramVtxZ_             = book1D(dir, "vertexZ",            "z_{vtx}",                                  100,        -25.,         +25.);
  histogramVtxRho_           = book1D(dir, "vertexRho",          "#rho_{vtx}",                               100,         -2.5,         +2.5);
  histogramRhoNeutral_       = book1D(dir, "rhoNeutral",         "#rho_{neutral}",                            50,          0. ,          12.);
}

void ZllRecoilCorrectionHistManager::fillHistograms(
       const reco::CompositeCandidate& ZllCand, const std::vector<pat::Muon>& muons, 
       const std::vector<pat::Jet>& jets, const pat::MET& met, 
       const reco::VertexCollection& vertices, double rhoNeutral, double evtWeight)
{
  assert(ZllCand.numberOfDaughters() == 2);

  reco::Candidate::LorentzVector p4LepPlus;
  if      ( ZllCand.daughter(0)->charge() > +0.5 ) p4LepPlus = ZllCand.daughter(0)->p4();
  else if ( ZllCand.daughter(1)->charge() > +0.5 ) p4LepPlus = ZllCand.daughter(1)->p4();
  else assert(0);

  histogramLepPlusPt_->Fill(p4LepPlus.pt(), evtWeight);
  histogramLepPlusEta_->Fill(p4LepPlus.eta(), evtWeight); 
  histogramLepPlusPhi_->Fill(p4LepPlus.phi(), evtWeight); 

  reco::Candidate::LorentzVector p4LepMinus;
  if      ( ZllCand.daughter(0)->charge() < -0.5 ) p4LepMinus = ZllCand.daughter(0)->p4();
  else if ( ZllCand.daughter(1)->charge() < -0.5 ) p4LepMinus = ZllCand.daughter(1)->p4();
  else assert(0);

  histogramLepMinusPt_->Fill(p4LepMinus.pt(), evtWeight);
  histogramLepMinusEta_->Fill(p4LepMinus.eta(), evtWeight); 
  histogramLepMinusPhi_->Fill(p4LepMinus.phi(), evtWeight); 

  histogramZllCandPt_->Fill(ZllCand.pt(), evtWeight);
  histogramZllCandEta_->Fill(ZllCand.eta(), evtWeight);
  histogramZllCandPhi_->Fill(ZllCand.phi(), evtWeight);
  histogramZllCandMass_->Fill(ZllCand.mass(), evtWeight);
  
  int vtxMultiplicity = vertices.size();
  
  int numJetsRawPtGt10 = 0;
  int numJetsCorrPtGt10 = 0;
  int numJetsRawPtGt12 = 0;
  int numJetsCorrPtGt12 = 0;
  int numJetsRawPtGt15 = 0;
  int numJetsCorrPtGt15 = 0;
  int numJetsRawPtGt17 = 0;
  int numJetsCorrPtGt17 = 0;
  int numJetsRawPtGt20 = 0;
  int numJetsCorrPtGt20 = 0;
  
  for ( std::vector<pat::Jet>::const_iterator jet = jets.begin();
	jet != jets.end(); ++jet ) {
    bool isMuonOverlap = false;
    for ( std::vector<pat::Muon>::const_iterator muon = muons.begin();
	  muon != muons.end(); ++muon ) {
      if ( deltaR(jet->p4(), muon->p4()) < 0.5 ) isMuonOverlap = true;
    }
    if ( isMuonOverlap ) continue;

    reco::Candidate::LorentzVector rawJetP4 = jet->correctedP4("Uncorrected");
    if ( rawJetP4.pt() > 10. ) ++numJetsRawPtGt10;
    if ( rawJetP4.pt() > 12. ) ++numJetsRawPtGt12;
    if ( rawJetP4.pt() > 15. ) ++numJetsRawPtGt15;
    if ( rawJetP4.pt() > 17. ) ++numJetsRawPtGt17;
    if ( rawJetP4.pt() > 20. ) ++numJetsRawPtGt20;
    
    reco::Candidate::LorentzVector corrJetP4 = jet->p4();
    // CV: require pseudo-rapidity of "raw" and corrected jet to match,
    //     in order to avoid problem with "unphysical" corrected jets of negative energy
    //     for which the energy becomes positive and the sign of eta "flips"
    //    (cf. https://hypernews.cern.ch/HyperNews/CMS/get/JetMET/1311/1.html)
    if ( TMath::Abs(rawJetP4.eta() - corrJetP4.eta()) < 0.1 ) {
      
      if      ( TMath::Abs(jet->eta()) < 1.1 ) histogramJetPtAbsEtaLt11_->Fill(jet->pt(), evtWeight);
      else if ( TMath::Abs(jet->eta()) < 1.7 ) histogramJetPtAbsEta11to17_->Fill(jet->pt(), evtWeight);
      else if ( TMath::Abs(jet->eta()) < 2.3 ) histogramJetPtAbsEta17to23_->Fill(jet->pt(), evtWeight);
      else                                     histogramJetPtAbsEtaGt23_->Fill(jet->pt(), evtWeight);
      
      if ( jet->genJet() && jet->genJet()->pt() > 0. ) {
	double jetRes = (jet->pt() - jet->genJet()->pt())/jet->genJet()->pt();
	if      ( TMath::Abs(jet->eta()) < 1.1 ) histogramJetResAbsEtaLt11_->Fill(jetRes, evtWeight);
	else if ( TMath::Abs(jet->eta()) < 1.7 ) histogramJetResAbsEta11to17_->Fill(jetRes, evtWeight);
	else if ( TMath::Abs(jet->eta()) < 2.3 ) histogramJetResAbsEta17to23_->Fill(jetRes, evtWeight);
	else                                     histogramJetResAbsEtaGt23_->Fill(jetRes, evtWeight);
      }
      
      if ( corrJetP4.pt() > 10. ) ++numJetsCorrPtGt10;
      if ( corrJetP4.pt() > 12. ) ++numJetsCorrPtGt12;
      if ( corrJetP4.pt() > 15. ) ++numJetsCorrPtGt15;
      if ( corrJetP4.pt() > 17. ) ++numJetsCorrPtGt17;
      if ( corrJetP4.pt() > 20. ) ++numJetsCorrPtGt20;
    }
  }

  //std::cout << " numJetsRawPtGt10 = " << numJetsRawPtGt10 << std::endl;
  //std::cout << " numJetsCorrPtGt10 = " << numJetsCorrPtGt10 << std::endl;

  histogramNumJetsRawPtGt10_->Fill(numJetsRawPtGt10, evtWeight);
  histogramNumJetsRawPtGt12_->Fill(numJetsRawPtGt12, evtWeight);
  histogramNumJetsRawPtGt15_->Fill(numJetsRawPtGt15, evtWeight);
  histogramNumJetsRawPtGt17_->Fill(numJetsRawPtGt17, evtWeight);
  histogramNumJetsRawPtGt20_->Fill(numJetsRawPtGt20, evtWeight);
  histogramNumJetsCorrPtGt10_->Fill(numJetsCorrPtGt10, evtWeight);
  histogramNumJetsCorrPtGt12_->Fill(numJetsCorrPtGt12, evtWeight);
  histogramNumJetsCorrPtGt15_->Fill(numJetsCorrPtGt15, evtWeight);
  histogramNumJetsCorrPtGt17_->Fill(numJetsCorrPtGt17, evtWeight);
  histogramNumJetsCorrPtGt20_->Fill(numJetsCorrPtGt20, evtWeight);

  histogramMEtS_->Fill(met.pt(), evtWeight);
  histogramMEtL_->Fill(met.pt(), evtWeight);

  if ( ZllCand.pt() > 0. && met.pt() < 60. ) { // CV: fit Z-recoil correction parameters for events with MEt < 60 GeV only,
                                               //     as di-boson and TTbar backgrounds dominate in high MEt tail
    double qT = ZllCand.pt();
    double qX = ZllCand.px();
    double qY = ZllCand.py();
    
    double metPx = met.px();
    double metPy = met.py();

    double metProjParlZ = (qX*metPx + qX*metPy)/qT;
    double metProjPerpZ = (qX*metPy - qY*metPx)/qT;
    
    histogramMEtProjParlZ_->Fill(metProjParlZ, evtWeight);
    histogramMEtProjPerpZ_->Fill(metProjPerpZ, evtWeight);

    int errorFlag = 0;
    std::pair<double, double> uT = compMEtProjU(ZllCand.p4(), met.px(), met.py(), errorFlag);
    if ( !errorFlag ) {
      double uParl = uT.first;
      double uPerp = uT.second;
      if ( qT > 0. ) histogramUparlDivQtVsQt_->Fill(qT, uParl/qT, evtWeight);
      histogramUparlVsQt_->Fill(qT, uParl, evtWeight);
      if ( qT > 0. ) histogramUperpDivQtVsQt_->Fill(qT, uPerp/qT, evtWeight);
      histogramUperpVsQt_->Fill(qT, uPerp, evtWeight);
      histogramQt_->Fill(qT, evtWeight);
     
      for ( std::vector<histogramsUvsQtNumObjType*>::iterator it = histogramsUvsQtNumVtxBinned_.begin();
	    it != histogramsUvsQtNumVtxBinned_.end(); ++it ) {
	if ( ((*it)->numObjMin_ == -1 || vtxMultiplicity >= (*it)->numObjMin_) &&
	     ((*it)->numObjMax_ == -1 || vtxMultiplicity <= (*it)->numObjMax_) ) {
	  if ( qT > 0. ) (*it)->histogramUparlDivQtVsQt_->Fill(qT, uParl/qT, evtWeight);
	  (*it)->histogramUparlVsQt_->Fill(qT, uParl, evtWeight);
	  if ( qT > 0. ) (*it)->histogramUperpDivQtVsQt_->Fill(qT, uPerp/qT, evtWeight);
	  (*it)->histogramUperpVsQt_->Fill(qT, uPerp, evtWeight);
	  (*it)->histogramQt_->Fill(qT, evtWeight);
	}
      }
      
      for ( std::vector<histogramsUvsQtNumObjType*>::iterator it = histogramsUvsQtNumJetsBinned_.begin();
	    it != histogramsUvsQtNumJetsBinned_.end(); ++it ) {
	if ( ((*it)->numObjMin_ == -1 || numJetsCorrPtGt10 >= (*it)->numObjMin_) &&
	     ((*it)->numObjMax_ == -1 || numJetsCorrPtGt10 <= (*it)->numObjMax_) ) {
	  if ( qT > 0. ) (*it)->histogramUparlDivQtVsQt_->Fill(qT, uParl/qT, evtWeight);
	  (*it)->histogramUparlVsQt_->Fill(qT, uParl, evtWeight);
	  if ( qT > 0. ) (*it)->histogramUperpDivQtVsQt_->Fill(qT, uPerp/qT, evtWeight);
	  (*it)->histogramUperpVsQt_->Fill(qT, uPerp, evtWeight);
	  (*it)->histogramQt_->Fill(qT, evtWeight);
	}
      }
    }
  }

  histogramVtxMultiplicity_->Fill(vtxMultiplicity, evtWeight);
  for ( reco::VertexCollection::const_iterator vertex = vertices.begin();
	vertex != vertices.end(); ++vertex ) {
    histogramVtxZ_->Fill(vertex->z(), evtWeight);
    histogramVtxRho_->Fill(vertex->position().Rho(), evtWeight);
  }
  histogramRhoNeutral_->Fill(rhoNeutral, evtWeight);
}

void ZllRecoilCorrectionHistManager::scaleHistograms(double factor)
{
  for ( std::vector<TH1*>::iterator histogram = histograms_.begin();
	histogram != histograms_.end(); ++histogram ) {
    if ( !(*histogram)->GetSumw2N() ) (*histogram)->Sumw2(); // CV: compute "proper" errors before scaling histogram
    (*histogram)->Scale(factor);
  }
}

TH1* ZllRecoilCorrectionHistManager::book1D(
       TFileDirectory& dir, const std::string& distribution, const std::string& title, int numBins, double min, double max)
{
  TH1* retVal = dir.make<TH1D>(distribution.data(), title.data(), numBins, min, max);
  histograms_.push_back(retVal);
  return retVal;
}

TH2* ZllRecoilCorrectionHistManager::book2D(
       TFileDirectory& dir, const std::string& distribution, const std::string& title, 
       int numBinsX, double xMin, double xMax, int numBinsY, double yMin, double yMax)				 
{
  TH2* retVal = dir.make<TH2D>(distribution.data(), title.data(), numBinsX, xMin, xMax, numBinsY, yMin, yMax);
  histograms_.push_back(retVal);
  return retVal;
}

TH2* ZllRecoilCorrectionHistManager::book2D(
       TFileDirectory& dir, const std::string& distribution, const std::string& title, 
       int numBinsX, double* xBinning, int numBinsY, double yMin, double yMax)				 
{
  TH2* retVal = dir.make<TH2D>(distribution.data(), title.data(), numBinsX, xBinning, numBinsY, yMin, yMax);
  histograms_.push_back(retVal);
  return retVal;
}
