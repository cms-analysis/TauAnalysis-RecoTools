#ifndef TauAnalysis_RecoTools_PATLeptonEfficiencyCorrectionProducer_h
#define TauAnalysis_RecoTools_PATLeptonEfficiencyCorrectionProducer_h

/** \class PATLeptonEfficiencyCorrectionProducer
 *
 * Produce downgrade weights to correct for differences in lepton reconstruction, id., 
 * isolation and trigger efficiencies between data and Monte Carlo simulation
 *
 * NOTE:
 *      o weight > 1: efficiency is underestimated by Monte Carlo simulation
 *      o weight = 1: efficiency is well modeled by Monte Carlo simulation
 *      o weight < 1: efficiency is overestimated by Monte Carlo simulation
 *
 * \authors Christian Veelken
 *
 * \version $Revision: 1.2 $
 *
 * $Id: PATLeptonEfficiencyCorrectionProducer.h,v 1.2 2010/09/28 11:23:36 jkolb Exp $
 *
 */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "CommonTools/Utils/interface/StringObjectFunction.h"

#include "TauAnalysis/CandidateTools/interface/FetchCollection.h"

#include <TFile.h>
#include <TH1.h>

#include <vector>
#include <string>

template<typename T>
class PATLeptonEfficiencyCorrectionProducer : public edm::EDProducer 
{
  typedef std::vector<edm::InputTag> vInputTag;

 public:

  explicit PATLeptonEfficiencyCorrectionProducer(const edm::ParameterSet& cfg)
    : inputFile_(0),
      lut_(0),
      stringFunctionX_(0),
      stringFunctionY_(0),
      stringFunctionZ_(0)
  {
    edm::FileInPath inputFileName = cfg.getParameter<edm::FileInPath>("inputFileName");
    std::string lutName = cfg.getParameter<std::string>("lutName");

    if ( inputFileName.isLocal()) {       
      inputFile_ = new TFile(inputFileName.fullPath().data());
      lut_ = (TH1*)inputFile_->Get(lutName.data());
    } else {
      throw cms::Exception("PATLeptonEfficiencyCorrectionProducer") 
	<< " Failed to find file = " << inputFileName << " !!";
    }

    edm::ParameterSet cfgParametrization = cfg.getParameter<edm::ParameterSet>("parametrization");
    srcLeptons_ = cfgParametrization.getParameter<vInputTag>("srcLeptons");

    stringFunctionX_ = new StringObjectFunction<T>(cfgParametrization.getParameter<std::string>("x"));
    if ( lut_->GetDimension() >= 2 ) stringFunctionY_ = new StringObjectFunction<T>(cfgParametrization.getParameter<std::string>("y"));
    if ( lut_->GetDimension() >= 3 ) stringFunctionZ_ = new StringObjectFunction<T>(cfgParametrization.getParameter<std::string>("z"));
    
    noLeptonSubstituteValue_ = cfg.getParameter<double>("noLeptonSubstituteValue");

    produces<double>("");
  }

  ~PATLeptonEfficiencyCorrectionProducer() 
  {
    delete inputFile_;

    delete stringFunctionX_;
    delete stringFunctionY_;
    delete stringFunctionZ_;
  }

  int getBin(TAxis* axis, double leptonX)
  {
    int bin = axis->FindFixBin(leptonX);
    int nBins = axis->GetNbins();
    if ( bin < 1     ) bin = 1;
    if ( bin > nBins ) bin = nBins;
    return bin;
  }

  double getEfficiencyCorrection(const T& lepton)
  {
    double leptonX, leptonY, leptonZ;
    int binX, binY, binZ;
    switch ( lut_->GetDimension() ) {
    case 1:
      leptonX = (*stringFunctionX_)(lepton);
      binX = getBin(lut_->GetXaxis(), leptonX);
      return lut_->GetBinContent(binX);
    case 2:
      leptonX = (*stringFunctionX_)(lepton);
      binX = getBin(lut_->GetXaxis(), leptonX);
      leptonY = (*stringFunctionY_)(lepton);
      binY = getBin(lut_->GetYaxis(), leptonY);
      return lut_->GetBinContent(binX, binY);
    case 3:
      leptonX = (*stringFunctionX_)(lepton);
      binX = getBin(lut_->GetXaxis(), leptonX);
      leptonY = (*stringFunctionY_)(lepton);
      binY = getBin(lut_->GetYaxis(), leptonY);
      leptonZ = (*stringFunctionZ_)(lepton);
      binZ = getBin(lut_->GetZaxis(), leptonZ);
      return lut_->GetBinContent(binX, binY, binZ);
    } 

    edm::LogError ("getEfficiencyCorrection") 
      << " Invalid dimensionality = " << lut_->GetDimension() << " of LUT histogram !!";
    return 1.;
  }

  void produce(edm::Event& evt, const edm::EventSetup& es)
  {
    //std::cout << "<PATLeptonEfficiencyCorrectionProducer::produce>:" << std::endl;

    double weight = 0.;

    typedef edm::View<T> TView;

    bool foundLepton = false;
    for ( vInputTag::const_iterator src = srcLeptons_.begin();
	  src != srcLeptons_.end(); ++src ) {
      edm::Handle<TView> leptonCollection;
      pf::fetchCollection(leptonCollection, *src, evt);

      for ( typename TView::const_iterator lepton = leptonCollection->begin();
	    lepton != leptonCollection->end(); ++lepton ) {
	weight = getEfficiencyCorrection(*lepton);
	foundLepton = true;
	break;
      }
    }

    if ( !foundLepton ) weight = noLeptonSubstituteValue_;
    //std::cout << "--> weight = " << weight << std::endl;
  
    std::auto_ptr<double> weightPtr(new double(weight));
    evt.put(weightPtr);
  }

 private:

  TFile* inputFile_;
  TH1* lut_;

  vInputTag srcLeptons_;

  StringObjectFunction<T>* stringFunctionX_;
  StringObjectFunction<T>* stringFunctionY_;
  StringObjectFunction<T>* stringFunctionZ_;

  double noLeptonSubstituteValue_;
};

#endif

