/** \executable FWLiteUnclusteredEnergyAnalyzer
 *
 * Determine data/MC residual corrections to "unclustered energy
 *
 * \author Christian Veelken, LLR
 *
 * \version $Revision: 1.1 $
 *
 * $Id: FWLiteUnclusteredEnergyAnalyzer.cc,v 1.1 2013/04/14 11:56:32 veelken Exp $
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

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"

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
typedef std::vector<double> vdouble;
typedef std::vector<edm::ParameterSet> vParameterSet;

struct binningEntryType
{
  binningEntryType(TTree* tree, const std::string& branchName_jetSum, const std::string& branchName_unclEnSum, const std::string& branchName_jetCorr, const std::string& branchName_jetCorr_offset, 
		   TFileDirectory& dir, const edm::ParameterSet& cfg, bool subtract_qT, double shiftBy)
    : branchName_jetSum_(branchName_jetSum),
      branchName_unclEnSum_(branchName_unclEnSum),    
      branchName_jetCorr_(branchName_jetCorr),
      branchName_jetCorr_offset_(branchName_jetCorr_offset),
      binLabel_(cfg.getParameter<std::string>("binLabel")),
      binCenter_(cfg.getParameter<double>("binCenter")),
      binEdgeLow_(cfg.getParameter<double>("binEdgeLow")),
      binEdgeHigh_(cfg.getParameter<double>("binEdgeHigh")),
      subtract_qT_(subtract_qT),
      shiftBy_(shiftBy)
  {
    //std::cout << "<binningEntryType::initialize>:" << std::endl;
    //std::cout << " binLabel = " << binLabel_ << std::endl;
    if ( branchName_jetSum_ != "" ) {
      //std::cout << "branchName_jetSum = " << branchName_jetSum << std::endl;
      std::string branchName_sumJetPx = Form("%s_%sPx", branchName_jetSum_.data(), binLabel_.data());
      tree->SetBranchAddress(branchName_sumJetPx.data(), &sumJetPx_);
      std::string branchName_sumJetPy = Form("%s_%sPy", branchName_jetSum_.data(), binLabel_.data());
      tree->SetBranchAddress(branchName_sumJetPy.data(), &sumJetPy_);
      std::string branchName_sumJetEt = Form("%s_%sSumEt", branchName_jetSum_.data(), binLabel_.data());
      tree->SetBranchAddress(branchName_sumJetEt.data(), &sumJetEt_);
    } else {
      sumJetPx_ = 0.;
      sumJetPy_ = 0.;
      sumJetEt_ = 0.;
    }
    if ( branchName_unclEnSum_ != "" ) {
      //std::cout << "branchName_unclEnSum = " << branchName_unclEnSum << std::endl;
      std::string branchName_sumUnclEnPx = Form("%s_%sPx", branchName_unclEnSum_.data(), binLabel_.data());
      tree->SetBranchAddress(branchName_sumUnclEnPx.data(), &sumUnclEnPx_);
      std::string branchName_sumUnclEnPy = Form("%s_%sPy", branchName_unclEnSum_.data(), binLabel_.data());
      tree->SetBranchAddress(branchName_sumUnclEnPy.data(), &sumUnclEnPy_);
      std::string branchName_sumUnclEnEt = Form("%s_%sSumEt", branchName_unclEnSum_.data(), binLabel_.data());
      tree->SetBranchAddress(branchName_sumUnclEnEt.data(), &sumUnclEnEt_);
    } else {
      throw cms::Exception("binningEntryType") 
	<< "Invalid Configuration Parameter 'branchName_unclEnSum' = " << branchName_unclEnSum << " !!\n";
    }
    if ( branchName_jetCorr != "" && branchName_jetCorr_offset != "" ) {
      //std::cout << "branchName_jetCorr = " << branchName_jetCorr << std::endl;
      std::string branchName_jetCorrPx = Form("%s_%sPx", branchName_jetCorr.data(), binLabel_.data());
      tree->SetBranchAddress(branchName_jetCorrPx.data(), &jetCorrPx_);
      std::string branchName_jetCorrPy = Form("%s_%sPy", branchName_jetCorr.data(), binLabel_.data());
      tree->SetBranchAddress(branchName_jetCorrPy.data(), &jetCorrPy_);
      std::string branchName_jetCorrSumEt = Form("%s_%sSumEt", branchName_jetCorr.data(), binLabel_.data());
      tree->SetBranchAddress(branchName_jetCorrSumEt.data(), &jetCorrSumEt_);
      //std::cout << "branchName_jetCorr_offset = " << branchName_jetCorr_offset << std::endl;
      std::string branchName_jetCorrPx_offset = Form("%s_%sPx", branchName_jetCorr_offset.data(), binLabel_.data());
      tree->SetBranchAddress(branchName_jetCorrPx_offset.data(), &jetCorrPx_offset_);
      std::string branchName_jetCorrPy_offset = Form("%s_%sPy", branchName_jetCorr_offset.data(), binLabel_.data());
      tree->SetBranchAddress(branchName_jetCorrPy_offset.data(), &jetCorrPy_offset_);
      std::string branchName_jetCorrSumEt_offset = Form("%s_%sSumEt", branchName_jetCorr_offset.data(), binLabel_.data());
      tree->SetBranchAddress(branchName_jetCorrSumEt_offset.data(), &jetCorrSumEt_offset_);
    } else {
      jetCorrPx_           = 0.;
      jetCorrPy_           = 0.;
      jetCorrSumEt_        = 0.;
      jetCorrPx_offset_    = 0.;
      jetCorrPy_offset_    = 0.;
      jetCorrSumEt_offset_ = 0.;
    }
  }
  ~binningEntryType() {}
  void setVerbosity(int verbosity)
  {
    verbosity_ = verbosity;
  }
  void update(double qX, double qY, const reco::Candidate::LorentzVector& muPlusP4, const reco::Candidate::LorentzVector& muMinusP4, double residualCorrection)
  {
    if ( verbosity_ ) {
      std::cout << "<binningEntryType::update>:" << std::endl;
      std::cout << " eta = " << binEdgeLow_ << ".." << binEdgeHigh_ << " (center = " << binCenter_ << ")" << std::endl;
      if ( subtract_qT_ ) {
	std::cout << "subtracting muPlus: Pt = " << muPlusP4.pt() << ", eta = " << muPlusP4.eta() << ", phi = " << muPlusP4.phi() 
		  << " (Px = " << muPlusP4.px() << ", Py = " << muPlusP4.py() << ")" << std::endl;
	std::cout << "subtracting muMinus: Pt = " << muMinusP4.pt() << ", eta = " << muMinusP4.eta() << ", phi = " << muMinusP4.phi() 
		  << " (Px = " << muMinusP4.px() << ", Py = " << muMinusP4.py() << ")" << std::endl;
      }
      std::cout << "qX = " << qX << ", qY = " << qY << " (qT = " << TMath::Sqrt(qX*qX + qY*qY) << ")" << std::endl;
      std::cout << "sumJet: Px = " << sumJetPx_ << ", Py = " << sumJetPy_ << ", sumEt = " << sumJetEt_ << std::endl;
      std::cout << "jetCorr: Px = " << jetCorrPx_ << ", Py = " << jetCorrPy_ << ", sumEt = " << jetCorrSumEt_ << std::endl;
      std::cout << "jetCorr(offset): Px = " << jetCorrPx_offset_ << ", Py = " << jetCorrPy_offset_ << ", sumEt = " << jetCorrSumEt_offset_ << std::endl;
      std::cout << "sumUnclEn: Px = " << sumUnclEnPx_ << ", Py = " << sumUnclEnPy_ << ", sumEt = " << sumUnclEnEt_ << std::endl;
    }

    sumJetPx_exclMuons_ = sumJetPx_;
    sumJetPy_exclMuons_ = sumJetPy_;
    sumJetEt_exclMuons_ = sumJetEt_;

    sumUnclEnPx_exclMuons_ = sumUnclEnPx_;
    sumUnclEnPy_exclMuons_ = sumUnclEnPy_;
    sumUnclEnEt_exclMuons_ = sumUnclEnEt_;

    if ( subtract_qT_ ) {
      if ( muPlusP4.eta() >= binEdgeLow_ && muPlusP4.eta() < binEdgeHigh_ ) {
	if ( branchName_jetSum_ != "" ) {
	  sumJetPx_exclMuons_ -= muPlusP4.px();
	  sumJetPy_exclMuons_ -= muPlusP4.py();
	  sumJetEt_exclMuons_ -= muPlusP4.pt();
	} else {
	  sumUnclEnPx_exclMuons_ -= muPlusP4.px();
	  sumUnclEnPy_exclMuons_ -= muPlusP4.py();
	  sumUnclEnEt_exclMuons_ -= muPlusP4.pt();
	}
      }
      if ( muMinusP4.eta() >= binEdgeLow_ && muMinusP4.eta() < binEdgeHigh_ ) {
	if ( branchName_jetSum_ != "" ) {
	  sumJetPx_exclMuons_ -= muMinusP4.px();
	  sumJetPy_exclMuons_ -= muMinusP4.py();
	  sumJetEt_exclMuons_ -= muMinusP4.pt();
	} else {
	  sumUnclEnPx_exclMuons_ -= muMinusP4.px();
	  sumUnclEnPy_exclMuons_ -= muMinusP4.py();
	  sumUnclEnEt_exclMuons_ -= muMinusP4.pt();
	}
      }      
    }

    sumHadPx_uncorrected_ = sumJetPx_exclMuons_ + (jetCorrPx_ - jetCorrPx_offset_) + sumUnclEnPx_exclMuons_;
    sumHadPy_uncorrected_ = sumJetPy_exclMuons_ + (jetCorrPy_ - jetCorrPy_offset_) + sumUnclEnPy_exclMuons_;
    sumHadEt_uncorrected_ = sumJetEt_exclMuons_ + (jetCorrSumEt_ - jetCorrSumEt_offset_) + sumUnclEnEt_exclMuons_;
    if ( verbosity_ ) {
      std::cout << " sumHad (before corrections): Px = " << sumHadPx_uncorrected_ << ", Py = " << sumHadPy_uncorrected_ << ", Et = " << sumHadEt_uncorrected_ << std::endl;
    }

    sumHadPx_corrected_ = sumHadPx_uncorrected_;
    sumHadPy_corrected_ = sumHadPy_uncorrected_;
    sumHadEt_corrected_ = sumHadEt_uncorrected_;
    sumHadPx_corrected_ += ((residualCorrection + shiftBy_ - 1.)*sumUnclEnPx_exclMuons_);
    sumHadPy_corrected_ += ((residualCorrection + shiftBy_ - 1.)*sumUnclEnPy_exclMuons_);
    sumHadEt_corrected_ += ((residualCorrection + shiftBy_ - 1.)*sumUnclEnEt_exclMuons_);
    if ( verbosity_ ) {
      std::cout << " sumHad (after corrections): Px = " << sumHadPx_corrected_ << ", Py = " << sumHadPy_corrected_ << ", Et = " << sumHadEt_corrected_ << std::endl;
    }

    errorFlag_ = 0;
    std::pair<double, double> uProj = compMEtProjU(qX, qY, -sumHadPx_corrected_, -sumHadPy_corrected_, errorFlag_, false);
    if ( !errorFlag_ ) {
      uParl_ = uProj.first;
      uPerp_ = uProj.second;
    }
    if ( verbosity_ ) {
      std::cout << "--> uParl = " << uParl_ << ", uPerp = " << uPerp_ << " (errorFlag = " << errorFlag_ << ")" << std::endl;
    }
    
    int dummy = 0;
    std::pair<double, double> uProj_unclustered = compMEtProjU(qX, qY, -sumUnclEnPx_exclMuons_, -sumUnclEnPy_exclMuons_, dummy, false);
    uParl_unclustered_ = uProj_unclustered.first;
    uPerp_unclustered_ = uProj_unclustered.second;
    uPx_unclCorr_ = (sumHadPx_corrected_ - sumHadPx_uncorrected_);
    uPy_unclCorr_ = (sumHadPy_corrected_ - sumHadPy_uncorrected_);
    std::pair<double, double> uProj_unclCorr = compMEtProjU(qX, qY, -uPx_unclCorr_, -uPy_unclCorr_, dummy, false);
    uParl_unclCorr_ = uProj_unclCorr.first;
    uPerp_unclCorr_ = uProj_unclCorr.second;
    std::pair<double, double> uProj_jets = compMEtProjU(qX, qY, -sumJetPx_exclMuons_, -sumJetPy_exclMuons_, dummy, false);
    uParl_jets_ = uProj_jets.first;
    uPerp_jets_ = uProj_jets.second;
    uPx_jetCorr_ = (jetCorrPx_ - jetCorrPx_offset_);
    uPy_jetCorr_ = (jetCorrPy_ - jetCorrPy_offset_);
    std::pair<double, double> uProj_jetCorr = compMEtProjU(qX, qY, -uPx_jetCorr_, -uPy_jetCorr_, dummy, false);
    uParl_jetCorr_ = uProj_jetCorr.first;
    uPerp_jetCorr_ = uProj_jetCorr.second;    
  }
  double uParl() const { return uParl_; }
  double uParl_unclustered() const { return uParl_unclustered_; } 
  double uParl_unclCorr() const { return uParl_unclCorr_; }
  double uParl_jets() const { return uParl_jets_; } 
  double uParl_jetCorr() const { return uParl_jetCorr_; } 
  double uPerp() const { return uPerp_; }
  double uPerp_unclustered() const { return uPerp_unclustered_; } 
  double uPerp_unclCorr() const { return uPerp_unclCorr_; } 
  double uPerp_jets() const { return uPerp_jets_; } 
  double uPerp_jetCorr() const { return uPerp_jetCorr_; } 
  double sumHadPx() const { return sumHadPx_corrected_; }
  double sumHadPy() const { return sumHadPy_corrected_; }
  double sumHadEt() const { return sumHadEt_corrected_; }
  std::string branchName_jetSum_;
  std::string branchName_unclEnSum_;
  std::string branchName_jetCorr_;
  std::string branchName_jetCorr_offset_;
  std::string binLabel_;
  double binCenter_;
  double binEdgeLow_;
  double binEdgeHigh_;
  bool subtract_qT_; // CV: true for PF, false for Calo 
  double shiftBy_;
  int verbosity_;
  int errorFlag_;
  Float_t sumJetPx_;
  Float_t sumJetPy_;
  Float_t sumJetEt_;
  Float_t sumJetPx_exclMuons_;
  Float_t sumJetPy_exclMuons_;
  Float_t sumJetEt_exclMuons_;
  Float_t sumUnclEnPx_;
  Float_t sumUnclEnPy_;
  Float_t sumUnclEnEt_;
  Float_t sumUnclEnPx_exclMuons_;
  Float_t sumUnclEnPy_exclMuons_;
  Float_t sumUnclEnEt_exclMuons_;
  Float_t jetCorrPx_;
  Float_t jetCorrPy_;
  Float_t jetCorrSumEt_;
  Float_t jetCorrPx_offset_;
  Float_t jetCorrPy_offset_;
  Float_t jetCorrSumEt_offset_;
  double sumHadPx_uncorrected_;
  double sumHadPy_uncorrected_;
  double sumHadEt_uncorrected_;
  double sumHadPx_corrected_;
  double sumHadPy_corrected_;
  double sumHadEt_corrected_;
  double uParl_;
  double uParl_unclustered_;
  double uParl_unclCorr_;
  double uParl_jets_;
  double uParl_jetCorr_;
  double uPerp_;
  double uPerp_unclustered_;
  double uPerp_unclCorr_;
  double uPerp_jets_;
  double uPerp_jetCorr_;  
  double uPx_unclCorr_;
  double uPx_jetCorr_;
  double uPy_unclCorr_;
  double uPy_jetCorr_;
};

struct histogramEntryType
{
  histogramEntryType(double numPileUpMin, double numPileUpMax, double qTmin, double qTmax)
    : numPileUpMin_(numPileUpMin),
      numPileUpMax_(numPileUpMax),
      qTmin_(qTmin),
      qTmax_(qTmax)
  {
    TString subDirectory_tstring = Form("numPileUp%1.1fto%1.1f_qT%1.1fto%1.1f", numPileUpMin_, numPileUpMax_, qTmin_, qTmax_);
    subDirectory_tstring.ReplaceAll(".", "_");
    subDirectory_ = subDirectory_tstring.Data();
  }
  ~histogramEntryType() {}
  void bookHistograms(TFileDirectory& dir, const TArrayD& etaBinning, const TArrayD& uParlBinning)
  {
    TFileDirectory subdir = dir.mkdir(subDirectory_);
    int etaNumBins = etaBinning.GetSize() - 1;
    int uParlNumBins = uParlBinning.GetSize() - 1;
    histogram_uParl_vs_eta_ = subdir.make<TH2D>("uParl_vs_eta", "u_{#parallel} vs #eta", etaNumBins, etaBinning.GetArray(), uParlNumBins, uParlBinning.GetArray());
    histogram_uParlDivQt_vs_eta_ = subdir.make<TH2D>("uParlDivQt_vs_eta", "u_{#parallel}/q_{T} vs #eta", etaNumBins, etaBinning.GetArray(), 401, -5.0, +5.0);
    histogram_uParl_ = subdir.make<TH1D>("uParl", "u_{#parallel}", uParlNumBins, uParlBinning.GetArray());
    histogram_uParlDivQt_ = subdir.make<TH1D>("uParlDivQt", "u_{#parallel}/q_{T}", 401, -5.0, +5.0);
    histogram_numPileUp_ = subdir.make<TH1D>("numPileUp", "N_{PU} (mean)", 600, 0., 60.);
    histogram_qT_ = subdir.make<TH1D>("qT", "q_{T}", 600, 0., 300.);
  }
  void fillHistograms(double numPileUp, double qT, const std::vector<binningEntryType*>& etaBinningEntries, double evtWeight)
  {
    if ( numPileUp > numPileUpMin_ && numPileUp < numPileUpMax_ && 
	 qT        > qTmin_        && qT        < qTmax_        ) {
      double uParl_total = 0.;
      for ( std::vector<binningEntryType*>::const_iterator etaBinningEntry = etaBinningEntries.begin();
	    etaBinningEntry != etaBinningEntries.end(); ++etaBinningEntry ) {
	histogram_uParl_vs_eta_->Fill((*etaBinningEntry)->binCenter_, (*etaBinningEntry)->uParl(), evtWeight);
	if ( qT > 5. ) histogram_uParlDivQt_vs_eta_->Fill((*etaBinningEntry)->binCenter_, (*etaBinningEntry)->uParl()/qT, evtWeight);
	uParl_total += (*etaBinningEntry)->uParl();
      }
      histogram_uParl_->Fill(uParl_total, evtWeight);
      if ( qT > 5. ) histogram_uParlDivQt_->Fill(uParl_total/qT, evtWeight);
      histogram_numPileUp_->Fill(numPileUp, evtWeight);
      histogram_qT_->Fill(qT, evtWeight);
    }
  }
  double numPileUpMin_;
  double numPileUpMax_;
  double qTmin_;
  double qTmax_;
  std::string subDirectory_;
  TH2* histogram_uParl_vs_eta_;
  TH2* histogram_uParlDivQt_vs_eta_;
  TH1* histogram_uParl_;
  TH1* histogram_uParlDivQt_;
  TH1* histogram_numPileUp_;
  TH1* histogram_qT_;
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

FactorizedJetCorrector* getResidualCorrector(const edm::FileInPath& residualCorrFileName)
{
  FactorizedJetCorrector* residualCorrector = 0;
  if ( !residualCorrFileName.isLocal()) 
    throw cms::Exception("FWLiteUnclusteredEnergyAnalyzer") 
      << " Failed to find File = " << residualCorrFileName << " !!\n";
  JetCorrectorParameters residualCorr(residualCorrFileName.fullPath().data());
  std::vector<JetCorrectorParameters> jetCorrections;
  jetCorrections.push_back(residualCorr);
  residualCorrector = new FactorizedJetCorrector(jetCorrections);
  return residualCorrector;
}

double getResidualCorrection(FactorizedJetCorrector* residualCorrector, double eta)
{
  residualCorrector->setJetEta(eta);
  residualCorrector->setJetPt(10.);
  residualCorrector->setJetA(0.25);
  residualCorrector->setRho(10.); 
  double residualCorrection = residualCorrector->getCorrection();
  return residualCorrection;
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

  std::string metType = cfgCalibUnclusteredEnergy.getParameter<std::string>("metType");

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
  //if ( !(numEvents >= 10000) ) 
  //  throw cms::Exception("FWLiteUnclusteredEnergyAnalyzer") 
  //    << " Failed to find Tree containing sufficient event statistics !!\n";

  fwlite::OutputFiles outputFile(cfg);
  fwlite::TFileService fs = fwlite::TFileService(outputFile.file().data());
  std::string directory = cfgCalibUnclusteredEnergy.getParameter<std::string>("directory");
  TFileDirectory dir = ( directory != "" ) ? fs.mkdir(directory) : fs;

  std::string branchName_recZ = cfgCalibUnclusteredEnergy.getParameter<std::string>("branchName_recZ");
  std::string branchName_recMuPlus = cfgCalibUnclusteredEnergy.getParameter<std::string>("branchName_recMuPlus");
  std::string branchName_recMuMinus = cfgCalibUnclusteredEnergy.getParameter<std::string>("branchName_recMuMinus");
  std::string branchName_met = cfgCalibUnclusteredEnergy.getParameter<std::string>("branchName_met");
  std::string branchName_jetSum = cfgCalibUnclusteredEnergy.getParameter<std::string>("branchName_jetSum");
  std::string branchName_unclEnSum = cfgCalibUnclusteredEnergy.getParameter<std::string>("branchName_unclEnSum");
  std::string branchName_jetCorr = cfgCalibUnclusteredEnergy.getParameter<std::string>("branchName_jetCorr");
  std::string branchName_jetCorr_offset = cfgCalibUnclusteredEnergy.getParameter<std::string>("branchName_jetCorr_offset");
  std::string branchName_numVertices = cfgCalibUnclusteredEnergy.getParameter<std::string>("branchName_numVertices");
  std::string branchName_zVtx = cfgCalibUnclusteredEnergy.getParameter<std::string>("branchName_zVtx");
  std::string branchName_numPileUp = cfgCalibUnclusteredEnergy.getParameter<std::string>("branchName_numPileUp");
  vstring branchNames_weights = cfgCalibUnclusteredEnergy.getParameter<vstring>("branchNames_weights");
  
  bool subtract_qT = cfgCalibUnclusteredEnergy.getParameter<bool>("subtract_qT");
  double shiftBy = cfgCalibUnclusteredEnergy.getParameter<double>("shiftBy");

  vParameterSet cfgEtaBinning = cfgCalibUnclusteredEnergy.getParameter<vParameterSet>("etaBinning");
 
  std::vector<double> etaBinning;
  for ( vParameterSet::const_iterator cfgBinningEntry = cfgEtaBinning.begin();
	cfgBinningEntry != cfgEtaBinning.end(); ++cfgBinningEntry ) {
    double binEdgeLow = cfgBinningEntry->getParameter<double>("binEdgeLow");
    double binEdgeHigh = cfgBinningEntry->getParameter<double>("binEdgeHigh");
    if ( etaBinning.size() == 0 ) etaBinning.push_back(binEdgeLow);
    else assert(TMath::Abs(etaBinning.back() - binEdgeLow) < 1.e-3);
    etaBinning.push_back(binEdgeHigh);
  }
  int etaNumBins = etaBinning.size() - 1;
  assert(etaNumBins >= 1);
  TArrayD etaBinning_array = convertToArray(etaBinning);

  std::vector<binningEntryType*> etaBinningEntries;
  for ( vParameterSet::const_iterator cfgBinningEntry = cfgEtaBinning.begin();
	cfgBinningEntry != cfgEtaBinning.end(); ++cfgBinningEntry ) {
    binningEntryType* etaBinningEntry = new binningEntryType(tree, branchName_jetSum, branchName_unclEnSum, branchName_jetCorr, branchName_jetCorr_offset, dir, *cfgBinningEntry, subtract_qT, shiftBy);
    etaBinningEntries.push_back(etaBinningEntry);
  }

  vdouble numPileUpBinning = cfgCalibUnclusteredEnergy.getParameter<vdouble>("numPileUpBinning");

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
  
  bool applyQtReweight = cfgCalibUnclusteredEnergy.getParameter<bool>("applyQtReweight");
  TFile* qTreweightInputFile = 0;
  TH1* qTreweightHistogram = 0;
  if ( applyQtReweight ) {
    edm::FileInPath qTreweightInputFileName = cfgCalibUnclusteredEnergy.getParameter<edm::FileInPath>("qTreweightInputFileName");
    if ( !qTreweightInputFileName.isLocal() ) 
      throw cms::Exception("FWLiteUnclusteredEnergyAnalyzer") 
	<< " Failed to find File = " << qTreweightInputFileName << " !!\n";
    std::string qTreweightHistogramName = cfgCalibUnclusteredEnergy.getParameter<std::string>("qTreweightHistogramName");
    qTreweightInputFile = new TFile(qTreweightInputFileName.fullPath().data());
    qTreweightHistogram = dynamic_cast<TH1*>(qTreweightInputFile->Get(qTreweightHistogramName.data()));
    if ( !qTreweightHistogram )
      throw cms::Exception("FWLiteUnclusteredEnergyAnalyzer") 
	<< " Failed to load qTreweightHistogram = " << qTreweightHistogramName.data() << " from file = " << qTreweightInputFileName.fullPath().data() << " !!\n";
  }

  bool applyResidualCorr = cfgCalibUnclusteredEnergy.getParameter<bool>("applyResidualCorr");
  FactorizedJetCorrector* residualCorrector_data_offset = 0;
  FactorizedJetCorrector* residualCorrector_data_slope  = 0;
  FactorizedJetCorrector* residualCorrector_mc_offset   = 0;
  FactorizedJetCorrector* residualCorrector_mc_slope    = 0;  
  if ( applyResidualCorr ) {
    edm::ParameterSet cfgResidualCorrFileNames = cfgCalibUnclusteredEnergy.getParameter<edm::ParameterSet>("residualCorrFileNames");
    edm::ParameterSet cfgResidualCorr_data = cfgResidualCorrFileNames.getParameter<edm::ParameterSet>("data");
    residualCorrector_data_offset = getResidualCorrector(cfgResidualCorr_data.getParameter<edm::FileInPath>("offset"));
    residualCorrector_data_slope  = getResidualCorrector(cfgResidualCorr_data.getParameter<edm::FileInPath>("slope"));
    edm::ParameterSet cfgResidualCorr_mc = cfgResidualCorrFileNames.getParameter<edm::ParameterSet>("mc");
    residualCorrector_mc_offset = getResidualCorrector(cfgResidualCorr_mc.getParameter<edm::FileInPath>("offset"));
    residualCorrector_mc_slope  = getResidualCorrector(cfgResidualCorr_mc.getParameter<edm::FileInPath>("slope"));
  }
  
  int verbosity = ( cfgCalibUnclusteredEnergy.exists("verbosity") ) ?
    cfgCalibUnclusteredEnergy.getParameter<int>("verbosity") : 0;
  for ( std::vector<binningEntryType*>::iterator etaBinningEntry = etaBinningEntries.begin();
	etaBinningEntry != etaBinningEntries.end(); ++etaBinningEntry ) {
    (*etaBinningEntry)->setVerbosity(verbosity);
  }

  const int qTnumBins = 34;
  double qTbinning[qTnumBins + 1] = { 
    0., 2.5, 5., 7.5, 10., 12.5, 15., 17.5, 20., 22.5, 25., 27.5, 30., 35., 40., 45., 50., 
    60., 70., 80., 90., 100., 110., 120., 130., 140., 150., 160., 170., 180., 200., 220., 240., 260., 300.
  };

  TH2* histogram_uParl_vs_qT                   = dir.make<TH2D>("uParl_vs_qT",                   "u_{#parallel} vs q_{T}",       qTnumBins, qTbinning, 231, -500.0, +75.0);
  TH2* histogram_uParlDivQt_vs_qT              = dir.make<TH2D>("uParlDivQt_vs_qT",              "u_{#parallel}/q_{T} vs q_{T}", qTnumBins, qTbinning, 401,   -5.0,  +5.0);
  TH2* histogram_uParl_unclustered_vs_qT       = dir.make<TH2D>("uParl_unclustered_vs_qT",       "u_{#parallel} vs q_{T}",       qTnumBins, qTbinning, 231, -500.0, +75.0);
  TH2* histogram_uParlDivQt_unclustered_vs_qT  = dir.make<TH2D>("uParlDivQt_unclustered_vs_qT",  "u_{#parallel}/q_{T} vs q_{T}", qTnumBins, qTbinning, 401,   -5.0,  +5.0);
  TH2* histogram_uParl_unclCorr_vs_qT          = dir.make<TH2D>("uParl_unclCorr_vs_qT",          "u_{#parallel} vs q_{T}",       qTnumBins, qTbinning, 231, -500.0, +75.0);
  TH2* histogram_uParlDivQt_unclCorr_vs_qT     = dir.make<TH2D>("uParlDivQt_unclCorr_vs_qT",     "u_{#parallel}/q_{T} vs q_{T}", qTnumBins, qTbinning, 401,   -5.0,  +5.0);
  TH2* histogram_uParl_jets_vs_qT              = dir.make<TH2D>("uParl_jets_vs_qT",              "u_{#parallel} vs q_{T}",       qTnumBins, qTbinning, 231, -500.0, +75.0);
  TH2* histogram_uParlDivQt_jets_vs_qT         = dir.make<TH2D>("uParlDivQt_jets_vs_qT",         "u_{#parallel}/q_{T} vs q_{T}", qTnumBins, qTbinning, 401,   -5.0,  +5.0);
  TH2* histogram_uParl_jetCorr_vs_qT           = dir.make<TH2D>("uParl_jetCorr_vs_qT",           "u_{#parallel} vs q_{T}",       qTnumBins, qTbinning, 231, -500.0, +75.0);
  TH2* histogram_uParlDivQt_jetCorr_vs_qT      = dir.make<TH2D>("uParlDivQt_jetCorr_vs_qT",      "u_{#parallel}/q_{T} vs q_{T}", qTnumBins, qTbinning, 401,   -5.0,  +5.0);
  TH1* histogramQt                             = dir.make<TH1D>("qT",                            "q_{T}",                        600, 0., 300.);
  TH1* histogramNumPileUp                      = dir.make<TH1D>("numPileUp",                     "N_{PU} (mean)",                600, 0.,  60.);

  std::vector<double> qTbinning_forReweighting;
  double qT_forReweighting = 0.;
  while ( qT_forReweighting <= 1000. ) {
    qTbinning_forReweighting.push_back(qT_forReweighting);
    if      ( qT_forReweighting <  10. ) qT_forReweighting +=   1.;
    else if ( qT_forReweighting <  20. ) qT_forReweighting +=   2.;
    else if ( qT_forReweighting <  40. ) qT_forReweighting +=   5.;
    else if ( qT_forReweighting <  60. ) qT_forReweighting +=  10.;
    else if ( qT_forReweighting < 100. ) qT_forReweighting +=  20.;
    else if ( qT_forReweighting < 150. ) qT_forReweighting +=  50.;
    else if ( qT_forReweighting < 300. ) qT_forReweighting += 150.;
    else                                 qT_forReweighting += 700.;
  }
  TArrayD qTbinning_forReweighting_array = convertToArray(qTbinning_forReweighting);
  int qTnumBins_forReweighting = qTbinning_forReweighting_array.GetSize() - 1;
  assert(qTnumBins_forReweighting >= 1);
  
  TH1* histogramQt_forReweighting = dir.make<TH1D>("qT_forReweighting", "q_{T}", qTnumBins_forReweighting, qTbinning_forReweighting_array.GetArray());

  std::vector<double> sumEtBinning;
  double sumEtMin = 0.0;
  double sumEtMax = 500.0;
  int sumEtNumBins = 500;
  double sumEtBinEdgeLow = sumEtMin;
  double sumEtBinWidth = (sumEtMax - sumEtMin)/sumEtNumBins;
  double sumEtBinEdgeHigh = sumEtBinEdgeLow + sumEtBinWidth;
  while ( sumEtBinEdgeHigh <= (sumEtMax + 1.e-3) ) {
    if ( sumEtBinning.size() == 0 ) sumEtBinning.push_back(sumEtBinEdgeLow);
    sumEtBinning.push_back(sumEtBinEdgeHigh);
    sumEtBinEdgeLow = sumEtBinEdgeHigh;
    sumEtBinEdgeHigh = sumEtBinEdgeLow + sumEtBinWidth;
  }
  assert((sumEtBinning.size() - 1) == sumEtNumBins);
  TArrayD sumEtBinning_array = convertToArray(sumEtBinning);
  
  TH1* histogramSumEt        = dir.make<TH1D>("sumEt",        "#Sigma E_{T}",         250, 0., 2500.);
  TH2* histogramSumEt_vs_eta = dir.make<TH2D>("sumEt_vs_eta", "#Sigma E_{T} vs #eta", etaNumBins, etaBinning_array.GetArray(), sumEtNumBins, sumEtBinning_array.GetArray());

  std::vector<double> uParlBinning;
  double uParlMin = -500.0;
  double uParlMax = +75.0;
  int uParlNumBins = 231;
  double uParlBinEdgeLow = uParlMin;
  double uParlBinWidth = (uParlMax - uParlMin)/uParlNumBins;
  double uParlBinEdgeHigh = uParlBinEdgeLow + uParlBinWidth;
  while ( uParlBinEdgeHigh <= (uParlMax + 1.e-3) ) {
    if ( uParlBinning.size() == 0 ) uParlBinning.push_back(uParlBinEdgeLow);
    uParlBinning.push_back(uParlBinEdgeHigh);
    uParlBinEdgeLow = uParlBinEdgeHigh;
    uParlBinEdgeHigh = uParlBinEdgeLow + uParlBinWidth;
  }
  assert((uParlBinning.size() - 1) == uParlNumBins);
  TArrayD uParlBinning_array = convertToArray(uParlBinning);

  std::vector<histogramEntryType*> histograms_vs_numPileUp_vs_qT;
  int numPileUpBins = numPileUpBinning.size() - 1;
  for ( int iBin_numPileUp = 0; iBin_numPileUp < numPileUpBins; ++iBin_numPileUp ) {
    double numPileUpMin = numPileUpBinning[iBin_numPileUp];
    double numPileUpMax = numPileUpBinning[iBin_numPileUp + 1];
    for ( int iBin_qT = 0; iBin_qT < qTnumBins; ++iBin_qT ) {
      double qTmin = qTbinning[iBin_qT];
      double qTmax = qTbinning[iBin_qT + 1];
      histogramEntryType* histogramEntry = new histogramEntryType(numPileUpMin, numPileUpMax, qTmin, qTmax);
      histogramEntry->bookHistograms(dir, etaBinning_array, uParlBinning_array);
      histograms_vs_numPileUp_vs_qT.push_back(histogramEntry);
    }
  }

  TH1* histogramZvtx = dir.make<TH1D>("zVtx", "zVtx", 500, -25., +25.);

  Int_t run, event, lumi;
  tree->SetBranchAddress("run", &run);
  tree->SetBranchAddress("event", &event);
  tree->SetBranchAddress("lumi", &lumi);
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
  Float_t metPx, metPy, metPt;
  tree->SetBranchAddress(Form("%sPx", branchName_met.data()), &metPx);
  tree->SetBranchAddress(Form("%sPy", branchName_met.data()), &metPy);
  tree->SetBranchAddress(Form("%sPt", branchName_met.data()), &metPt);
  Int_t numVertices;
  tree->SetBranchAddress(branchName_numVertices.data(), &numVertices);
  Float_t zVtx;
  tree->SetBranchAddress(branchName_zVtx.data(), &zVtx);
  Float_t numPileUp;
  tree->SetBranchAddress(branchName_numPileUp.data(), &numPileUp);

  int numEvents_selected = 0;
  for ( int iEvent = 0; iEvent < numEvents && (maxEvents == -1 || iEvent < maxEvents); ++iEvent ) {
    if ( iEvent > 0 && (iEvent % 100000) == 0 ) {
      std::cout << "processing Event " << iEvent << std::endl;
    }

    tree->GetEntry(iEvent);

    //if ( !(qT > 80. && qT < 90.) ) continue; // CV: ONLY FOR TESTING !!

    if ( verbosity ) {
      std::cout << "processing Event " << iEvent << ": run# = " << run << ", ls# = " << lumi << ", event# = " << event << std::endl;
      std::cout << " qT = " << qT << " (qX = " << qX << ", qY = " << qY << "), numPileUp = " << numPileUp << std::endl;
      std::cout << " met: Pt = " << metPt << " (Px = " << metPx << ", Py = " << metPy << ")" << std::endl;
      bool subtract_qT_met = ( branchName_met.find("calo") != std::string::npos ) ?
	false : true;
      int errorFlag = 0;      
      std::pair<double, double> uProj = compMEtProjU(qX, qY, metPx, metPy, errorFlag, subtract_qT_met);
      if ( !errorFlag ) {
	double uParl = uProj.first;
	double uPerp = uProj.second;
	std::cout << "uParl(met) = " << uParl << ", uPerp = " << uPerp << std::endl;
      }
    }

    reco::Candidate::LorentzVector muPlusP4(muPlusPx, muPlusPy, muPlusPz, muPlusEn);
    reco::Candidate::LorentzVector muMinusP4(muMinusPx, muMinusPy, muMinusPz, muMinusEn);
    for ( std::vector<binningEntryType*>::iterator etaBinningEntry = etaBinningEntries.begin();
	  etaBinningEntry != etaBinningEntries.end(); ++etaBinningEntry ) {
      double residualCorrection = 1.0;
      if ( applyResidualCorr ) {
	double residualCorrParameter_data_offset = getResidualCorrection(residualCorrector_data_offset, (*etaBinningEntry)->binCenter_);
	double residualCorrParameter_data_slope  = getResidualCorrection(residualCorrector_data_slope, (*etaBinningEntry)->binCenter_);
	double response_data = residualCorrParameter_data_offset + residualCorrParameter_data_slope*numPileUp;
	double residualCorrParameter_mc_offset   = getResidualCorrection(residualCorrector_mc_offset, (*etaBinningEntry)->binCenter_);
	double residualCorrParameter_mc_slope    = getResidualCorrection(residualCorrector_mc_slope, (*etaBinningEntry)->binCenter_);
	double response_mc = residualCorrParameter_mc_offset + residualCorrParameter_mc_slope*numPileUp;
	if ( response_mc > 0. ) {
	  residualCorrection = response_data/response_mc; // CV: correct Monte Carlo to match response measured in Data
	  if ( verbosity ) {
	    std::cout << "response(eta = " << (*etaBinningEntry)->binCenter_ << "): data = " << response_data << ", mc = " << response_mc 
	  	      << " --> residualCorrection = " << residualCorrection << std::endl;
	  }
	}
      }
      (*etaBinningEntry)->update(qX, qY, muPlusP4, muMinusP4, residualCorrection);
    }

    double sumHadPx_total = 0.;
    double sumHadPy_total = 0.;
    double sumHadEt_total = 0.;    
    for ( std::vector<binningEntryType*>::iterator etaBinningEntry = etaBinningEntries.begin();
	  etaBinningEntry != etaBinningEntries.end(); ++etaBinningEntry ) {
      sumHadPx_total += (*etaBinningEntry)->sumHadPx();
      sumHadPy_total += (*etaBinningEntry)->sumHadPy();
      sumHadEt_total += (*etaBinningEntry)->sumHadEt();
    }

    double evtWeight = 1.0;
    for ( std::vector<weightEntryType*>::const_iterator weightEntry = weightEntries.begin();
	  weightEntry != weightEntries.end(); ++weightEntry ) {
      if ( verbosity ) {
	std::cout << "branchName = " << (*weightEntry)->branchName_ << ": weight = " << (*weightEntry)->weight_ << std::endl;
      }
      evtWeight *= (*weightEntry)->weight_;
    }
    if ( applyZvtxReweight ) {
      int bin = zVtxReweightHistogram->FindBin(zVtx);
      if ( bin <=  1                                 ) bin = 1;
      if ( bin >= zVtxReweightHistogram->GetNbinsX() ) bin = zVtxReweightHistogram->GetNbinsX();
      if ( verbosity ) {
	std::cout << "bin(zVtx) = " << bin << ": weight = " << zVtxReweightHistogram->GetBinContent(bin) << std::endl;
      }
      evtWeight *= zVtxReweightHistogram->GetBinContent(bin);
    }
    if ( applyQtReweight ) {
      int bin = qTreweightHistogram->FindBin(qT);
      if ( bin <=  1                               ) bin = 1;
      if ( bin >= qTreweightHistogram->GetNbinsX() ) bin = qTreweightHistogram->GetNbinsX();
      if ( verbosity ) {
	std::cout << "bin(qT) = " << bin << ": weight = " << qTreweightHistogram->GetBinContent(bin) << std::endl;
      }
      evtWeight *= qTreweightHistogram->GetBinContent(bin);
    }
    if ( verbosity ) {
      std::cout << "evtWeight = " << evtWeight << std::endl;
    }
    if ( evtWeight < 1.e-3 || evtWeight > 1.e+3 ) continue;

    if ( numVertices >= 1 ) {
      int errorFlag = 0;
      std::pair<double, double> uProj = compMEtProjU(qX, qY, -sumHadPx_total, -sumHadPy_total, errorFlag, false);
      if ( !errorFlag ) {
	double uParl = uProj.first;
	double uPerp = uProj.second;
	//std::cout << "uParl = " << uParl << ", uPerp = " << uPerp << std::endl;
	histogram_uParl_vs_qT->Fill(qT, uParl, evtWeight);
	if ( qT > 5. ) histogram_uParlDivQt_vs_qT->Fill(qT, uParl/qT, evtWeight);
	double uParl_unclustered = 0.;
	double uPerp_unclustered = 0.;
	double uParl_unclCorr    = 0.;
	double uPerp_unclCorr    = 0.;
	double uParl_jets        = 0.;
	double uPerp_jets        = 0.;
	double uParl_jetCorr     = 0.;
	double uPerp_jetCorr     = 0.;
	double uParl_sum         = 0.;
	double uPerp_sum         = 0.;
	double uPx_unclustered   = 0.;
	double uPy_unclustered   = 0.;
	double uPx_unclCorr      = 0.;
	double uPy_unclCorr      = 0.;
	double uPx_jets          = 0.;
	double uPy_jets          = 0.;
	double uPx_jetCorr       = 0.;
	double uPy_jetCorr       = 0.;
	double uPx_sum           = 0.;
	double uPy_sum           = 0.;
	for ( std::vector<binningEntryType*>::iterator etaBinningEntry = etaBinningEntries.begin();
	      etaBinningEntry != etaBinningEntries.end(); ++etaBinningEntry ) {
	  uParl_unclustered += (*etaBinningEntry)->uParl_unclustered();
	  uPerp_unclustered += (*etaBinningEntry)->uPerp_unclustered();
	  uParl_unclCorr    += (*etaBinningEntry)->uParl_unclCorr();
	  uPerp_unclCorr    += (*etaBinningEntry)->uPerp_unclCorr();
	  uParl_jets        += (*etaBinningEntry)->uParl_jets();
	  uPerp_jets        += (*etaBinningEntry)->uPerp_jets();
	  uParl_jetCorr     += (*etaBinningEntry)->uParl_jetCorr();
	  uPerp_jetCorr     += (*etaBinningEntry)->uPerp_jetCorr();
	  if ( verbosity ) {
	    uParl_sum       += (*etaBinningEntry)->uParl();
	    uPerp_sum       += (*etaBinningEntry)->uPerp();
	    uPx_unclustered += (*etaBinningEntry)->sumUnclEnPx_;
	    uPy_unclustered += (*etaBinningEntry)->sumUnclEnPy_;
	    uPx_unclCorr    += (*etaBinningEntry)->uPx_unclCorr_;
	    uPy_unclCorr    += (*etaBinningEntry)->uPy_unclCorr_;
	    uPx_jets        += (*etaBinningEntry)->sumJetPx_;
	    uPy_jets        += (*etaBinningEntry)->sumJetPy_;
	    uPx_jetCorr     += (*etaBinningEntry)->uPx_jetCorr_;
	    uPy_jetCorr     += (*etaBinningEntry)->uPy_jetCorr_;
	    uPx_sum         += (*etaBinningEntry)->sumHadPx();
	    uPy_sum         += (*etaBinningEntry)->sumHadPy();
	    //std::cout << "eta = " << (*etaBinningEntry)->binCenter_ << ":" << std::endl;
	    //std::cout << " uPx: unclustered = " << (*etaBinningEntry)->sumUnclEnPx_ << "," 
	    //	        << " unclCorr = " << (*etaBinningEntry)->uPx_unclCorr_ << "," 
	    //	        << " jets = " << (*etaBinningEntry)->sumJetPx_ << "," 
	    //	        << " jetCorr = " << (*etaBinningEntry)->uPx_jetCorr_ << std::endl;	  
	    //std::cout << " uPy: unclustered = " << (*etaBinningEntry)->sumUnclEnPy_ << "," 
	    //	        << " unclCorr = " << (*etaBinningEntry)->uPy_unclCorr_ << "," 
	    //          << " jets = " << (*etaBinningEntry)->sumJetPy_ << "," 
	    //          << " jetCorr = " << (*etaBinningEntry)->uPy_jetCorr_ << std::endl;
	  }
	}
	if ( verbosity ) {
	  std::cout << "sum:" << std::endl;
	  std::cout << " uParl: unclustered = " << uParl_unclustered << ", unclCorr = " << uParl_unclCorr << ", jets = " << uParl_jets << ", jetCorr = " << uParl_jetCorr << ","
		    << " sum-q = " << (uParl_unclustered + uParl_unclCorr + uParl_jets + uParl_jetCorr + qT) 
		    << " (" << uParl_sum << ")" << std::endl;
	  std::cout << " uPerp: unclustered = " << uPerp_unclustered << ", unclCorr = " << uPerp_unclCorr << ", jets = " << uPerp_jets << ", jetCorr = " << uPerp_jetCorr << ","
		    << " sum-q = " << (uPerp_unclustered + uPerp_unclCorr + uPerp_jets + uPerp_jetCorr) 
		    << " (" << uPerp_sum << ")" << std::endl;
	  std::cout << " uPx: unclustered = " << uPx_unclustered << ", unclCorr = " << uPx_unclCorr << ", jets = " << uPx_jets << ", jetCorr = " << uPx_jetCorr << ","
		    << " sum-q = " << (uPx_unclustered + uPx_unclCorr + uPx_jets + uPx_jetCorr - qX) 
		    << " (" << uPx_sum << ")" << std::endl;
	  std::cout << " uPy: unclustered = " << uPy_unclustered << ", unclCorr = " << uPy_unclCorr << ", jets = " << uPy_jets << ", jetCorr = " << uPy_jetCorr << ","
		    << " sum-q = " << (uPy_unclustered + uPy_unclCorr + uPy_jets + uPy_jetCorr - qY) 
		    << " (" << uPy_sum << ")" << std::endl;
	}
	histogram_uParl_unclustered_vs_qT->Fill(qT, uParl_unclustered, evtWeight);	
	histogram_uParl_unclCorr_vs_qT->Fill(qT, uParl_unclCorr, evtWeight);
	histogram_uParl_jets_vs_qT->Fill(qT, uParl_jets, evtWeight);
	histogram_uParl_jetCorr_vs_qT->Fill(qT, uParl_jetCorr, evtWeight);
	if ( qT > 5. ) {
	  histogram_uParlDivQt_unclustered_vs_qT->Fill(qT, uParl_unclustered/qT, evtWeight);
	  histogram_uParlDivQt_unclCorr_vs_qT->Fill(qT, uParl_unclCorr/qT, evtWeight);
	  histogram_uParlDivQt_jets_vs_qT->Fill(qT, uParl_jets/qT, evtWeight);
	  histogram_uParlDivQt_jetCorr_vs_qT->Fill(qT, uParl_jetCorr/qT, evtWeight);
	}
	histogramQt->Fill(qT, evtWeight); 
	histogramQt_forReweighting->Fill(qT, evtWeight); 
       	histogramSumEt->Fill(sumHadEt_total, evtWeight); 
	histogramNumPileUp->Fill(numPileUp, evtWeight); 
      }
      for ( std::vector<histogramEntryType*>::iterator histogramEntry = histograms_vs_numPileUp_vs_qT.begin();
	    histogramEntry != histograms_vs_numPileUp_vs_qT.end(); ++histogramEntry ) {
	(*histogramEntry)->fillHistograms(numPileUp, qT, etaBinningEntries, evtWeight);
      }
      for ( std::vector<binningEntryType*>::iterator etaBinningEntry = etaBinningEntries.begin();
	    etaBinningEntry != etaBinningEntries.end(); ++etaBinningEntry ) {
	histogramSumEt_vs_eta->Fill((*etaBinningEntry)->binCenter_, (*etaBinningEntry)->sumHadEt(), evtWeight); 
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

  for ( std::vector<binningEntryType*>::iterator it = etaBinningEntries.begin();
	it != etaBinningEntries.end(); ++it ) {
    delete (*it);
  }

  for ( std::vector<weightEntryType*>::iterator it = weightEntries.begin();
	it != weightEntries.end(); ++it ) {
    delete (*it);
  }

  for ( std::vector<histogramEntryType*>::iterator it = histograms_vs_numPileUp_vs_qT.begin();
	it != histograms_vs_numPileUp_vs_qT.end(); ++it ) {
    delete (*it);
  }
  
  delete zVtxReweightInputFile;
  //delete zVtxReweightHistogram;
  delete qTreweightInputFile;
  //delete qTreweightHistogram;

  clock.Show("FWLiteUnclusteredEnergyAnalyzer");

  return 0;
}
