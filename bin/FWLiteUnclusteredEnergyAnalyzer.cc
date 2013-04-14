
/** \executable FWLiteUnclusteredEnergyAnalyzer
 *
 * Determine data/MC residual corrections to "unclustered energy
 *
 * \author Christian Veelken, LLR
 *
 * \version $Revision: 1.13 $
 *
 * $Id: FWLiteUnclusteredEnergyAnalyzer.cc,v 1.13 2012/08/28 15:01:36 veelken Exp $
 *
 */

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"

#include "FWCore/Utilities/interface/Exception.h"

#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "DataFormats/FWLite/interface/InputSource.h"
#include "DataFormats/FWLite/interface/OutputFiles.h"

#include "DataFormats/Candidate/interface/Candidate.h"

#include "TauAnalysis/CandidateTools/interface/candidateAuxFunctions.h"

#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TBenchmark.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TArrayD.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TSystem.h>
#include <TROOT.h>

#include <vector>
#include <string>

typedef std::vector<std::string> vstring;
typedef std::vector<edm::ParameterSet> vParameterSet;

struct binningEntryType
{
  binningEntryType(TTree* tree, const std::string& branchName, TFileDirectory& dir, const edm::ParameterSet& cfg, bool subtract_qT, double shiftBy)
    : binLabel_(cfg.getParameter<std::string>("binLabel")),
      binCenter_(cfg.getParameter<double>("binCenter")),
      binEdgeLow_(cfg.getParameter<double>("binEdgeLow")),
      binEdgeHigh_(cfg.getParameter<double>("binEdgeHigh")),
      subtract_qT_(subtract_qT),
      shiftBy_(shiftBy)
  {
    //std::cout << "<binningEntryType::binningEntryType>:" << std::endl;
    //std::cout << " binLabel = " << binLabel_ << std::endl;
    tree->SetBranchAddress(Form("%s_%sPx", branchName.data(), binLabel_.data()), &sumPx_);
    tree->SetBranchAddress(Form("%s_%sPy", branchName.data(), binLabel_.data()), &sumPy_);
    tree->SetBranchAddress(Form("%s_%sSumEt", branchName.data(), binLabel_.data()), &sumEt_);
    std::string histogramName_metParl_vs_qT = Form("metParl_vs_qT_%s", binLabel_.data());
    histogram_metParl_vs_qT_ = dir.make<TH2D>(histogramName_metParl_vs_qT.data(), histogramName_metParl_vs_qT.data(), 300, 0., 300., 500, -400., +100.);
    std::string histogramName_metParl_vs_uParl = Form("metParl_vs_uParl_%s", binLabel_.data());
    histogram_metParl_vs_uParl_ = dir.make<TH2D>(histogramName_metParl_vs_uParl.data(), histogramName_metParl_vs_uParl.data(), 500, -400., +100., 500, -400., +100.);
  }
  ~binningEntryType() 
  {
    //delete histogram_metParl_vs_qT_;
    //delete histogram_metParl_vs_uParl_;
  }
  void update(double qX, double qY, const reco::Candidate::LorentzVector& muPlusP4, const reco::Candidate::LorentzVector& muMinusP4)
  {
    sumPx_corrected_ = sumPx_;
    sumPy_corrected_ = sumPy_;
    sumEt_corrected_ = sumEt_;
    if ( subtract_qT_ ) {
      if ( muPlusP4.eta() > binEdgeLow_ && muPlusP4.eta() < binEdgeHigh_ ) {
	sumPx_corrected_ -= muPlusP4.px();
	sumPy_corrected_ -= muPlusP4.py();
	sumEt_corrected_ -= muPlusP4.pt();
      }
      if ( muMinusP4.eta() > binEdgeLow_ && muMinusP4.eta() < binEdgeHigh_ ) {
	sumPx_corrected_ -= muMinusP4.px();
	sumPy_corrected_ -= muMinusP4.py();
	sumEt_corrected_ -= muMinusP4.pt();
      }      
    }
    sumPx_corrected_ *= (1. + shiftBy_);
    sumPy_corrected_ *= (1. + shiftBy_);
    sumEt_corrected_ *= (1. + shiftBy_);
    errorFlag_ = 0;
    std::pair<double, double> uProj = compMEtProjU(qX, qY, -sumPx_corrected_, -sumPy_corrected_, errorFlag_, false);
    if ( errorFlag_ ) return;
    uParl_ = uProj.first;
    uPerp_ = uProj.second;
  }
  void fillHistograms(double qX, double qY, double qT, double metX, double metY, double evtWeight) 
  {
    int errorFlag = 0;
    std::pair<double, double> metProj = compMEtProjU(qX, qY, metX, metY, errorFlag, false);
    if ( errorFlag ) return;
    metParl_ = metProj.first;
    metPerp_ = metProj.second;
    histogram_metParl_vs_qT_->Fill(qT, metParl_, evtWeight);
    histogram_metParl_vs_uParl_->Fill(uParl_, metParl_, evtWeight);
  }
  double uParl() const { return uParl_; }
  double uPerp() const { return uPerp_; }
  double sumPx() const { return sumPx_corrected_; }
  double sumPy() const { return sumPy_corrected_; }
  double sumEt() const { return sumEt_corrected_; }
  double metParl() const { return metParl_; }
  double metPerp() const { return metPerp_; }
  std::string binLabel_;
  double binCenter_;
  double binEdgeLow_;
  double binEdgeHigh_;
  bool subtract_qT_; // CV: true for PF, false for Calo 
  double shiftBy_;
  int errorFlag_;
  Float_t sumPx_;
  Float_t sumPy_;
  Float_t sumEt_;
  double sumPx_corrected_;
  double sumPy_corrected_;
  double sumEt_corrected_;
  double uParl_;
  double uPerp_;
  double metParl_;
  double metPerp_;
  TH2* histogram_metParl_vs_qT_;
  TH2* histogram_metParl_vs_uParl_;
};

struct weightEntryType
{
  weightEntryType(TTree* tree, const std::string& branchName)
    : branchName_(branchName)
  {
    tree->SetBranchAddress(branchName_.data(), &weight_);
  }
  ~weightEntryType() 
  {}
  std::string branchName_;
  Float_t weight_;
};

TArrayD convertToArray(const std::vector<double>& binning)
{
  TArrayD binning_array(binning.size());
  for ( size_t i = 0; i < binning.size(); ++i ) {
    binning_array[i] = binning[i];
  }
  return binning_array;
}

int main(int argc, char* argv[]) 
{
//--- parse command-line arguments
  if ( argc < 2 ) {
    std::cout << "Usage: " << argv[0] << " [parameters.py]" << std::endl;
    return 0;
  }

  std::cout << "<FWLiteUnclusteredEnergyAnalyzer>:" << std::endl;

//--- load framework libraries
  gSystem->Load("libFWCoreFWLite");
  AutoLibraryLoader::enable();

//--- keep track of time it takes the macro to execute
  TBenchmark clock;
  clock.Start("FWLiteUnclusteredEnergyAnalyzer");

//--- read python configuration parameters
  if ( !edm::readPSetsFrom(argv[1])->existsAs<edm::ParameterSet>("process") ) 
    throw cms::Exception("FWLiteUnclusteredEnergyAnalyzer") 
      << "No ParameterSet 'process' found in configuration file = " << argv[1] << " !!\n";

  edm::ParameterSet cfg = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("process");

  edm::ParameterSet cfgCalibUnclusteredEnergy = cfg.getParameter<edm::ParameterSet>("UnclusteredEnergyAnalyzer");

  std::string type = cfgCalibUnclusteredEnergy.getParameter<std::string>("type");

  std::string treeName = cfgCalibUnclusteredEnergy.getParameter<std::string>("treeName");

  fwlite::InputSource inputFiles(cfg); 
  int maxEvents = inputFiles.maxEvents();
  std::cout << "maxEvents = " << maxEvents << std::endl;

  std::vector<TTree*> treesToDelete;

  TChain* tree = new TChain(treeName.c_str());  
  for ( vstring::const_iterator inputFileName = inputFiles.files().begin();
	inputFileName != inputFiles.files().end(); ++inputFileName ) {
    tree->AddFile(inputFileName->c_str());
  }
  treesToDelete.push_back(tree);
  std::cout << "numInputFiles = " << inputFiles.files().size() << std::endl;

  int numEvents = tree->GetEntries();  
  std::cout << "numEvents = " << numEvents << std::endl;
  if ( !(numEvents >= 10000) ) 
    throw cms::Exception("FWLiteUnclusteredEnergyAnalyzer") 
      << " Failed to find Tree containing sufficient event statistics !!\n";

  fwlite::OutputFiles outputFile(cfg);
  fwlite::TFileService fs = fwlite::TFileService(outputFile.file().data());
  std::string directory = cfgCalibUnclusteredEnergy.getParameter<std::string>("directory");
  TFileDirectory dir = ( directory != "" ) ? fs.mkdir(directory) : fs;

  std::string branchName_recZ = cfgCalibUnclusteredEnergy.getParameter<std::string>("branchName_recZ");
  std::string branchName_recMuPlus = cfgCalibUnclusteredEnergy.getParameter<std::string>("branchName_recMuPlus");
  std::string branchName_recMuMinus = cfgCalibUnclusteredEnergy.getParameter<std::string>("branchName_recMuMinus");
  std::string branchName_met = cfgCalibUnclusteredEnergy.getParameter<std::string>("branchName_met");
  std::string branchName_unclEnSum = cfgCalibUnclusteredEnergy.getParameter<std::string>("branchName_unclEnSum");
  std::string branchName_numVertices = cfgCalibUnclusteredEnergy.getParameter<std::string>("branchName_numVertices");
  std::string branchName_zVtx = cfgCalibUnclusteredEnergy.getParameter<std::string>("branchName_zVtx");
  vstring branchNames_weights = cfgCalibUnclusteredEnergy.getParameter<vstring>("branchNames_weights");
  
  bool subtract_qT = cfgCalibUnclusteredEnergy.getParameter<bool>("subtract_qT");
  double shiftBy = cfgCalibUnclusteredEnergy.getParameter<double>("shiftBy");

  vParameterSet cfgBinning = cfgCalibUnclusteredEnergy.getParameter<vParameterSet>("binning");

  std::vector<double> etaBinning;
  for ( vParameterSet::const_iterator cfgBinningEntry = cfgBinning.begin();
	cfgBinningEntry != cfgBinning.end(); ++cfgBinningEntry ) {
    double binEdgeLow = cfgBinningEntry->getParameter<double>("binEdgeLow");
    double binEdgeHigh = cfgBinningEntry->getParameter<double>("binEdgeHigh");
    if ( etaBinning.size() == 0 ) etaBinning.push_back(binEdgeLow);
    else assert(TMath::Abs(etaBinning.back() - binEdgeLow) < 1.e-3);
    etaBinning.push_back(binEdgeHigh);
  }
  int etaNumBins = etaBinning.size() - 1;
  assert(etaNumBins >= 1);
  TArrayD etaBinning_array = convertToArray(etaBinning);

  std::vector<binningEntryType*> binningEntries;
  for ( vParameterSet::const_iterator cfgBinningEntry = cfgBinning.begin();
	cfgBinningEntry != cfgBinning.end(); ++cfgBinningEntry ) {
    binningEntryType* binningEntry = new binningEntryType(tree, branchName_unclEnSum, dir, *cfgBinningEntry, subtract_qT, shiftBy);
    binningEntries.push_back(binningEntry);
  }
  std::vector<weightEntryType*> weightEntries;
  for ( vstring::const_iterator branchName_weight = branchNames_weights.begin();
	branchName_weight != branchNames_weights.end(); ++branchName_weight ) {
    weightEntryType* weightEntry = new weightEntryType(tree, *branchName_weight);
    weightEntries.push_back(weightEntry);
  }

  bool applyZvtxReweight = cfgCalibUnclusteredEnergy.getParameter<bool>("applyZvtxReweight");
  TFile* zVtxReweightInputFile = 0;
  TH1* zVtxReweightHistogram = 0;
  if ( applyZvtxReweight ) {
    edm::FileInPath zVtxReweightInputFileName = cfgCalibUnclusteredEnergy.getParameter<edm::FileInPath>("zVtxReweightInputFileName");
    if ( !zVtxReweightInputFileName.isLocal() ) 
      throw cms::Exception("FWLiteUnclusteredEnergyAnalyzer") 
	<< " Failed to find File = " << zVtxReweightInputFileName << " !!\n";
    std::string zVtxReweightHistogramName = cfgCalibUnclusteredEnergy.getParameter<std::string>("zVtxReweightHistogramName");
    zVtxReweightInputFile = new TFile(zVtxReweightInputFileName.fullPath().data());
    zVtxReweightHistogram = dynamic_cast<TH1*>(zVtxReweightInputFile->Get(zVtxReweightHistogramName.data()));
    if ( !zVtxReweightHistogram )
      throw cms::Exception("FWLiteUnclusteredEnergyAnalyzer") 
	<< " Failed to load zVtxReweightHistogram = " << zVtxReweightHistogramName.data() << " from file = " << zVtxReweightInputFileName.fullPath().data() << " !!\n";
  }
  
  const int qTnumBins = 34;
  double qTbinning[qTnumBins + 1] = { 
    0., 2.5, 5., 7.5, 10., 12.5, 15., 17.5, 20., 22.5, 25., 27.5, 30., 35., 40., 45., 50., 
    60., 70., 80., 90., 100., 110., 120., 130., 140., 150., 160., 170., 180., 200., 220., 240., 260., 300.
  };

  TH2* histogram_uParl_vs_qT_absZvtxLt5 = dir.make<TH2D>("uParl_vs_qT_absZvtxLt5",   "u_{#parallel} vs q_{T}", qTnumBins, qTbinning, 230, -500.0, +75.0);
  TH2* histogram_uParlDivQt_vs_qT_absZvtxLt5 = dir.make<TH2D>("uParlDivQt_vs_qT_absZvtxLt5",   "u_{#parallel}/q_{T} vs q_{T}", qTnumBins, qTbinning, 400, -5.0, +5.0);
  TH1* histogramQt_absZvtxLt5 = dir.make<TH1D>("qT_absZvtxLt5", "q_{T}", 600, 0., 300.);
  TH2* histogram_uParl_vs_qT_absZvtx5to15 = dir.make<TH2D>("uParl_vs_qT_absZvtx5to15", "u_{#parallel} vs q_{T}", qTnumBins, qTbinning, 230, -500.0, +75.0);
  TH2* histogram_uParlDivQt_vs_qT_absZvtx5to15 = dir.make<TH2D>("uParlDivQt_vs_qT_absZvtx5to15", "u_{#parallel}/q_{T} vs q_{T}", qTnumBins, qTbinning, 400, -5.0, +5.0);
  TH1* histogramQt_absZvtx5to15 = dir.make<TH1D>("qT_absZvtx5to15", "q_{T}", 600, 0., 300.);
  TH2* histogram_uParl_vs_qT_absZvtxGt15 = dir.make<TH2D>("uParl_vs_qT_absZvtxGt15",  "u_{#parallel} vs q_{T}", qTnumBins, qTbinning, 230, -500.0, +75.0);
  TH2* histogram_uParlDivQt_vs_qT_absZvtxGt15 = dir.make<TH2D>("uParlDivQt_vs_qT_absZvtxGt15",  "u_{#parallel}/q_{T} vs q_{T}", qTnumBins, qTbinning, 400, -5.0, +5.0);
  TH1* histogramQt_absZvtxGt15 = dir.make<TH1D>("qT_absZvtxGt15", "q_{T}", 600, 0., 300.);
  TH2* histogram_uParl_vs_qT = dir.make<TH2D>("uParl_vs_qT",  "u_{#parallel} vs q_{T}", qTnumBins, qTbinning, 230, -500.0, +75.0);
  TH2* histogram_uParlDivQt_vs_qT = dir.make<TH2D>("uParlDivQt_vs_qT",  "u_{#parallel}/q_{T} vs q_{T}", qTnumBins, qTbinning, 400, -5.0, +5.0);
  TH1* histogramQt = dir.make<TH1D>("qT", "q_{T}", 600, 0., 300.);

  std::vector<double> uParlBinning;
  double uParlMin = -500.0;
  double uParlMax = +75.0;
  int uParlNumBins = 230;
  double binEdgeLow = -500.0;
  double binWidth = (uParlMax - uParlMin)/uParlNumBins;
  double binEdgeHigh = binEdgeLow + binWidth;
  while ( binEdgeHigh <= (uParlMax + 1.e-3) ) {
    if ( uParlBinning.size() == 0 ) uParlBinning.push_back(binEdgeLow);
    uParlBinning.push_back(binEdgeHigh);
    binEdgeLow = binEdgeHigh;
    binEdgeHigh = binEdgeLow + binWidth;
  }
  assert((uParlBinning.size() - 1) == uParlNumBins);
  TArrayD uParlBinning_array = convertToArray(uParlBinning);

  TH3* histogram_uParl_vs_eta_vs_qT = 
    new TH3D("uParl_vs_eta_vs_qT", "u_{#parallel} vs #eta vs q_{T}", 
	     qTnumBins, qTbinning, etaNumBins, etaBinning_array.GetArray(), uParlNumBins, uParlBinning_array.GetArray());
  TH3* histogram_uParlDivQt_vs_eta_vs_qT = 
    new TH3D("uParlDivQt_vs_eta_vs_qT", "u_{#parallel}/q_{T} vs #eta vs q_{T}", 
	     qTnumBins, qTbinning, etaNumBins, etaBinning_array.GetArray(), uParlNumBins, uParlBinning_array.GetArray());

  TH1* histogramZvtx = dir.make<TH1D>("zVtx", "zVtx", 500, -25., +25.);

  Float_t qX, qY, qT;
  tree->SetBranchAddress(Form("%sPx", branchName_recZ.data()), &qX);
  tree->SetBranchAddress(Form("%sPy", branchName_recZ.data()), &qY);
  tree->SetBranchAddress(Form("%sPt", branchName_recZ.data()), &qT);
  Float_t muPlusEn, muPlusPx, muPlusPy, muPlusPz;
  tree->SetBranchAddress(Form("%sEn", branchName_recMuPlus.data()), &muPlusEn);
  tree->SetBranchAddress(Form("%sPx", branchName_recMuPlus.data()), &muPlusPx);
  tree->SetBranchAddress(Form("%sPy", branchName_recMuPlus.data()), &muPlusPy);
  tree->SetBranchAddress(Form("%sPz", branchName_recMuPlus.data()), &muPlusPz);
  Float_t muMinusEn, muMinusPx, muMinusPy, muMinusPz;
  tree->SetBranchAddress(Form("%sEn", branchName_recMuMinus.data()), &muMinusEn);
  tree->SetBranchAddress(Form("%sPx", branchName_recMuMinus.data()), &muMinusPx);
  tree->SetBranchAddress(Form("%sPy", branchName_recMuMinus.data()), &muMinusPy);
  tree->SetBranchAddress(Form("%sPz", branchName_recMuMinus.data()), &muMinusPz);
  Float_t metPx, metPy;
  tree->SetBranchAddress(Form("%sPx", branchName_met.data()), &metPx);
  tree->SetBranchAddress(Form("%sPy", branchName_met.data()), &metPy);
  Int_t numVertices;
  tree->SetBranchAddress(branchName_numVertices.data(), &numVertices);
  Float_t zVtx;
  tree->SetBranchAddress(branchName_zVtx.data(), &zVtx);

  int numEvents_selected = 0;
  for ( int iEvent = 0; iEvent < numEvents && (maxEvents == -1 || iEvent < maxEvents); ++iEvent ) {
    if ( iEvent > 0 && (iEvent % 100000) == 0 ) {
      std::cout << "processing Event " << iEvent << std::endl;
    }
    tree->GetEntry(iEvent);
    //std::cout << " qT = " << qT << " (qX = " << qX << ", qY = " << qY << ")" << std::endl;
    reco::Candidate::LorentzVector muPlusP4(muPlusPx, muPlusPy, muPlusPz, muPlusEn);
    reco::Candidate::LorentzVector muMinusP4(muMinusPx, muMinusPy, muMinusPz, muMinusEn);
    for ( std::vector<binningEntryType*>::iterator binningEntry = binningEntries.begin();
	  binningEntry != binningEntries.end(); ++binningEntry ) {
      (*binningEntry)->update(qX, qY, muPlusP4, muMinusP4);
    }
    double metPx_summed = 0.;
    double metPy_summed = 0.;
    double sumEt_total = 0.;    
    binningEntryType* binningEntry_maxActivity = 0;
    double sumEt_maxActivity = 0.;
    for ( std::vector<binningEntryType*>::iterator binningEntry = binningEntries.begin();
	  binningEntry != binningEntries.end(); ++binningEntry ) {
      metPx_summed -= (*binningEntry)->sumPx();
      metPy_summed -= (*binningEntry)->sumPy();
      sumEt_total += (*binningEntry)->sumEt();
      if ( (*binningEntry)->sumEt() > sumEt_maxActivity ) {
	binningEntry_maxActivity = (*binningEntry);
	sumEt_maxActivity = (*binningEntry)->sumEt();
      }
    }
    double evtWeight = 1.0;
    for ( std::vector<weightEntryType*>::const_iterator weightEntry = weightEntries.begin();
	  weightEntry != weightEntries.end(); ++weightEntry ) {
      evtWeight *= (*weightEntry)->weight_;
    }
    if ( applyZvtxReweight ) {
      int bin = zVtxReweightHistogram->FindBin(zVtx);
      if ( bin <=  1                                 ) bin = 1;
      if ( bin >= zVtxReweightHistogram->GetNbinsX() ) bin = zVtxReweightHistogram->GetNbinsX();
      evtWeight *= zVtxReweightHistogram->GetBinContent(bin);
    }
    //std::cout << "evtWeight = " << evtWeight << std::endl;
    if ( evtWeight < 1.e-3 || evtWeight > 1.e+3 ) continue;
    if ( qT > 10. &&
	 sumEt_maxActivity > (0.25*sumEt_total) && 
	 sumEt_maxActivity < (0.75*sumEt_total) ) {      
      binningEntry_maxActivity->fillHistograms(qX, qY, qT, metPx_summed, metPy_summed, evtWeight);
      //std::cout << "qT = " << qT << ": uParl = " << binningEntry_maxActivity->uParl() << ", metParl = " << binningEntry_maxActivity->metParl() << std::endl;
    }
    if ( numVertices >= 1 ) {
      int errorFlag = 0;
      std::pair<double, double> uProj = compMEtProjU(qX, qY, metPx, metPy, errorFlag, subtract_qT);
      if ( !errorFlag ) {
	double uParl = uProj.first;
	double uPerp = uProj.second;
	double absZvtx = TMath::Abs(zVtx);
	if ( absZvtx < 5. ) {
	  histogram_uParl_vs_qT_absZvtxLt5->Fill(qT, uParl, evtWeight);
	  histogramQt_absZvtxLt5->Fill(qT, evtWeight);
	} else if ( absZvtx < 15. ) {
	  histogram_uParl_vs_qT_absZvtx5to15->Fill(qT, uParl, evtWeight);
	  histogramQt_absZvtx5to15->Fill(qT, evtWeight);
	} else {
	  histogram_uParl_vs_qT_absZvtxGt15->Fill(qT, uParl, evtWeight);
	  histogramQt_absZvtxGt15->Fill(qT, evtWeight);
	}	
	histogram_uParl_vs_qT->Fill(qT, uParl, evtWeight);
	histogramQt->Fill(qT, evtWeight);      
      }
      for ( std::vector<binningEntryType*>::iterator binningEntry = binningEntries.begin();
	    binningEntry != binningEntries.end(); ++binningEntry ) {
	//std::cout << "eta = " << (*binningEntry)->binCenter_ << ": qT = " << qT << ", uParl = " << (*binningEntry)->uParl() << std::endl;
	histogram_uParl_vs_eta_vs_qT->Fill(qT, (*binningEntry)->binCenter_, (*binningEntry)->uParl(), evtWeight);
	histogram_uParlDivQt_vs_eta_vs_qT->Fill(qT, (*binningEntry)->binCenter_, (*binningEntry)->uParl()/qT, evtWeight);
      }
      histogramZvtx->Fill(zVtx, evtWeight);
    }
    ++numEvents_selected;
  }

  std::cout << "numEvents_selected = " << numEvents_selected << std::endl;

  for ( std::vector<TTree*>::iterator it = treesToDelete.begin();
	it != treesToDelete.end(); ++it ) {
    delete (*it);
  }

  for ( std::vector<binningEntryType*>::iterator it = binningEntries.begin();
	it != binningEntries.end(); ++it ) {
    delete (*it);
  }

  delete zVtxReweightInputFile;
  //delete zVtxReweightHistogram;

  clock.Show("FWLiteUnclusteredEnergyAnalyzer");

  return 0;
}
