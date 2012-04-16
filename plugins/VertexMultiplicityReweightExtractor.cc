#include "TauAnalysis/RecoTools/plugins/VertexMultiplicityReweightExtractor.h"

#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "TauAnalysis/CandidateTools/interface/generalAuxFunctions.h"

#include <TMath.h>

#include <vector>

//-------------------------------------------------------------------------------
// CV: code to computed genPileUp reweights
//     copied & pasted from SimGeneral/MixingModule/python/mix_E7TeV_FlatDist10_2011EarlyData_inTimeOnly_cfi.py
//     cf. https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupReweighting
//    (extended up to 36 pile-up interactions by just adding zeros at the end)
namespace 
{
  std::vector<float> gen_pileup_flat10_summer11mc()
  {
    const int gen_pu_max = 36;
    const double npu_probs[gen_pu_max] = {
      0.0698146584, 0.0698146584, 0.0698146584, 0.0698146584, 0.0698146584,
      0.0698146584, 0.0698146584, 0.0698146584, 0.0698146584, 0.0698146584, 0.0698146584 /* <-- 10*/,
      0.0630151648, 0.0526654164, 0.0402754482, 0.0292988928, 0.0194384503, 0.0122016783, 0.007207042, 0.004003637, 0.0020278322,
      0.0010739954, 0.0004595759, 0.0002229748, 0.0001028162, 4.58337152809607E-05, /* <-- 24 */
      0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 
      0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,
      0.0000000000
    };
    std::vector<float> retVal(gen_pu_max);
    for ( int i = 0; i < gen_pu_max; ++i ) {
      retVal[i] = npu_probs[i];
    }
    return retVal;
  }

  int getBin(TAxis* axis, double x)
  {
    int bin = axis->FindFixBin(x);
    int nBins = axis->GetNbins();
    if ( bin < 1     ) bin = 1;
    if ( bin > nBins ) bin = nBins;
    return bin;
  }
}
//-------------------------------------------------------------------------------

VertexMultiplicityReweightExtractor::VertexMultiplicityReweightExtractor(const edm::ParameterSet& cfg)
  : src_(cfg.getParameter<edm::InputTag>("src")),
    inputFile_(0),
    genLumiReweight_(0),
    genLumiReweight3d_(0),
    recVtxMultiplicityReweight_(0),
    type_(kUndefined)
{
  edm::FileInPath inputFileName = cfg.getParameter<edm::FileInPath>("inputFileName");
  std::string lutName = cfg.getParameter<std::string>("lutName");
  if ( !inputFileName.isLocal()) throw cms::Exception("VertexMultiplicityReweightExtractor") 
    << " Failed to find File = " << inputFileName << " !!\n";

  std::string type_string = cfg.getParameter<std::string>("type");
  if      ( type_string == "gen"   ) type_ = kGenLevel;
  else if ( type_string == "gen3d" ) type_ = kGenLevel3d;
  else if ( type_string == "rec"   ) type_ = kRecLevel;
  else throw cms::Exception("VertexMultiplicityReweightExtractor") 
    << " Invalid Configuration Parameter 'type' = " << type_string << " !!\n";

  if ( type_ == kGenLevel || type_ == kGenLevel3d ) {
    inputFile_ = new TFile(inputFileName.fullPath().data());
    puDist_data_ = dynamic_cast<TH1*>(inputFile_->Get(lutName.data()));
    if ( !puDist_data_ ) 
      throw cms::Exception("VertexMultiplicityReweightExtractor") 
	<< " Failed to load LUT = " << lutName.data() << " from file = " << inputFileName.fullPath().data() << " !!\n";
    
    std::vector<float> gen_pileup_mc = gen_pileup_flat10_summer11mc();
    TH1* puDist_mc = new TH1D("MC_distr", "MC_distr", gen_pileup_mc.size(), 0., gen_pileup_mc.size());
    int numBins = puDist_mc->GetNbinsX();
    for ( int iBin = 1; iBin <= numBins; ++iBin ) {
      double binCenter = puDist_mc->GetBinCenter(iBin);
      int idx = TMath::FloorNint(binCenter);
      double binContent = ( idx >= 0 && idx < (int)gen_pileup_mc.size() ) ?
	gen_pileup_mc[idx] : 0.;
      puDist_mc->SetBinContent(iBin, binContent);
      puDist_mc->SetBinError(iBin, 0.);
    }
    std::string puFileName_mc = "puDist_mc.root";
    TFile* puFile_mc = new TFile(puFileName_mc.data(), "RECREATE");
    puDist_mc->Write();
    delete puFile_mc;
    
    if ( type_ == kGenLevel ) {
      genLumiReweight_ = new edm::LumiReWeighting(puFileName_mc.data(), inputFileName.fullPath().data(), "MC_distr", lutName.data());
    } else if ( type_ == kGenLevel3d ) {
      genLumiReweight3d_ = new edm::Lumi3DReWeighting(puFileName_mc.data(), inputFileName.fullPath().data(), "MC_distr", lutName.data());
      // CV: use pp inelastic cross-section of 73.5mb measured by TOTEM
      //     instead of CMS measurement of 68mb (default in Lumi3DReWeighting)
      genLumiReweight3d_->weight3D_init(73.5/68.);
    }
  } else {
    inputFile_ = new TFile(inputFileName.fullPath().data());
    recVtxMultiplicityReweight_ = dynamic_cast<TH1*>(inputFile_->Get(lutName.data()));
    if ( !recVtxMultiplicityReweight_ ) 
      throw cms::Exception("VertexMultiplicityReweightExtractor") 
	<< " Failed to load LUT = " << lutName.data() << " from file = " << inputFileName.fullPath().data() << " !!\n";
  }
}

VertexMultiplicityReweightExtractor::~VertexMultiplicityReweightExtractor() 
{
  delete inputFile_;
  delete genLumiReweight_;
  delete genLumiReweight3d_;
  delete recVtxMultiplicityReweight_;
}

double VertexMultiplicityReweightExtractor::operator()(const edm::Event& evt) const
{
  double weight = 1.;
  
  if ( type_ == kGenLevel || type_ == kGenLevel3d ) {
    typedef std::vector<PileupSummaryInfo> PileupSummaryInfoCollection;
    edm::Handle<PileupSummaryInfoCollection> genPileUpInfos;
    evt.getByLabel(src_, genPileUpInfos);

    int numPileUp_inTime     = -1;
    int numPileUp_bxPrevious = -1;
    int numPileUp_bxNext     = -1;
    for ( PileupSummaryInfoCollection::const_iterator genPileUpInfo = genPileUpInfos->begin();
	  genPileUpInfo != genPileUpInfos->end(); ++genPileUpInfo ) {
      // CV: in-time PU is stored in getBunchCrossing = 0, 
      //    cf. https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupInformation
      int bx = genPileUpInfo->getBunchCrossing();
      if      ( bx ==  0 ) numPileUp_inTime     = genPileUpInfo->getPU_NumInteractions();
      else if ( bx == -1 ) numPileUp_bxPrevious = genPileUpInfo->getPU_NumInteractions();
      else if ( bx == +1 ) numPileUp_bxNext     = genPileUpInfo->getPU_NumInteractions();
    }
    if ( numPileUp_bxPrevious == -1 || numPileUp_inTime == -1 || numPileUp_bxNext == -1 ) 
      throw cms::Exception("VertexMultiplicityReweightExtractor") 
	<< " Failed to find PileupSummaryInfo object for in-time Pile-up !!\n";
    if      ( type_ == kGenLevel   ) weight = genLumiReweight_->weight(numPileUp_inTime);
    else if ( type_ == kGenLevel3d ) weight = genLumiReweight3d_->weight3D(numPileUp_bxPrevious, numPileUp_inTime, numPileUp_bxNext);
    else assert(0);
  } else if ( type_ == kRecLevel ) {
    edm::Handle<reco::VertexCollection> vertices;
    evt.getByLabel(src_, vertices);
    int binIdx = getBin(recVtxMultiplicityReweight_->GetXaxis(), vertices->size());
    weight = recVtxMultiplicityReweight_->GetBinContent(binIdx);
  }

  //std::cout << " numPileUp = " << numPileUp << ": weight = " << weight << std::endl;

  return weight;
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(ObjValExtractorPluginFactory, VertexMultiplicityReweightExtractor, "VertexMultiplicityReweightExtractor");

