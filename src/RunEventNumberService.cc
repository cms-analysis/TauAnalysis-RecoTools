#include "TauAnalysis/RecoTools/interface/RunEventNumberService.h"

#include "TauAnalysis/DQMTools/interface/dqmAuxFunctions.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DQMServices/Core/interface/DQMStore.h"

#include <iostream>
#include <iomanip>

void checkKeyword(const std::string& filterCondition, const char* keyword, bool& isMatched, const std::string& dqmDirectory_store, 
		  bool& doSaveRunEventNumbers, std::string& dqmDirectory_saveRunEventNumbers)
{
  if ( filterCondition == keyword ) {
    doSaveRunEventNumbers = true;
    dqmDirectory_saveRunEventNumbers = dqmDirectoryName(dqmDirectory_store).append("events_").append(keyword);
    isMatched = true;
  }
}

RunEventNumberService::RunEventNumberService(const edm::ParameterSet& cfg)
{
  //std::cout << "<RunEventNumberService::RunEventNumberService>:" << std::endl; 

  cfgError_ = 0;

  name_ = cfg.getParameter<std::string>("name");
  //std::cout << " name = " << name_ << std::endl;

  dqmDirectory_store_ = cfg.getParameter<std::string>("dqmDirectory_store");
  //std::cout << " dqmDirectory_store = " << dqmDirectory_store_ << std::endl;

  typedef std::vector<edm::ParameterSet> vParameterSet;

  if ( cfg.exists("config") ) {
    vParameterSet cfgFilters = cfg.getParameter<vParameterSet>("config");
    for ( vParameterSet::const_iterator cfgFilter = cfgFilters.begin(); 
	  cfgFilter != cfgFilters.end(); ++cfgFilter ) {
      std::string filterName = cfgFilter->getParameter<std::string>("filterName");
      typedef std::vector<std::string> vstring;
      vstring saveRunEventNumbers = ( cfgFilter->exists("saveRunEventNumbers") ) ? 
	cfgFilter->getParameter<vstring>("saveRunEventNumbers") : vstring();

      std::string dqmDirectory_filter = dqmDirectoryName(dqmDirectory_store_).append(filterName);

      filterConfigEntry entry;
      entry.doSaveRunEventNumbers_passed_ = false;
      entry.dqmDirectory_passed_ = "";
      entry.doSaveRunEventNumbers_rejected_ = false;
      entry.dqmDirectory_rejected_ = "";
      entry.doSaveRunEventNumbers_exclRejected_ = false;
      entry.dqmDirectory_exclRejected_ = "";
      entry.doSaveRunEventNumbers_passed_cumulative_ = false;
      entry.dqmDirectory_passed_cumulative_ = "";
      entry.doSaveRunEventNumbers_rejected_cumulative_ = false;
      entry.dqmDirectory_rejected_cumulative_ = "";
      
      for ( vstring::const_iterator filterCondition = saveRunEventNumbers.begin();
	    filterCondition != saveRunEventNumbers.end(); ++filterCondition ) {
	bool isMatched = false;
	checkKeyword(*filterCondition, "passed", isMatched, dqmDirectory_filter,
		     entry.doSaveRunEventNumbers_passed_, entry.dqmDirectory_passed_);
	checkKeyword(*filterCondition, "rejected", isMatched, dqmDirectory_filter,
		     entry.doSaveRunEventNumbers_rejected_, entry.dqmDirectory_rejected_);
	checkKeyword(*filterCondition, "exclRejected", isMatched, dqmDirectory_filter,
		     entry.doSaveRunEventNumbers_exclRejected_, entry.dqmDirectory_exclRejected_);
	checkKeyword(*filterCondition, "passed_cumulative", isMatched, dqmDirectory_filter,
		     entry.doSaveRunEventNumbers_passed_cumulative_, entry.dqmDirectory_passed_cumulative_);
	checkKeyword(*filterCondition, "rejected_cumulative", isMatched, dqmDirectory_filter,
		     entry.doSaveRunEventNumbers_rejected_cumulative_, entry.dqmDirectory_rejected_cumulative_);
	if ( !isMatched ) {
	  edm::LogError ("RunEventNumberService") << " Undefined filterCondition = " << (*filterCondition) 
						  << " in Configuration Parameter saveRunEventNumbers !!";
	  cfgError_ = 1;
	}
      }

      filterConfigs_.insert(std::pair<std::string, filterConfigEntry>(filterName, entry));
    }
  }

  //for ( std::map<std::string, filterConfigEntry>::const_iterator filterConfig = filterConfigs_.begin();
  //	  filterConfig != filterConfigs_.end(); ++filterConfig ) {
  //  std::cout << " eventSelection = " << filterConfig->first << std::endl;
  //  std::cout << "  doSaveRunEventNumbers_passed = " << filterConfig->second.doSaveRunEventNumbers_passed_ << std::endl;
  //  std::cout << "   dqmDirectory_passed = " << filterConfig->second.dqmDirectory_passed_ << std::endl;
  //  std::cout << "  doSaveRunEventNumbers_rejected = " << filterConfig->second.doSaveRunEventNumbers_rejected_ << std::endl;
  //  std::cout << "   dqmDirectory_rejected = " << filterConfig->second.dqmDirectory_rejected_ << std::endl;
  //  std::cout << "  doSaveRunEventNumbers_exclRejected = " << filterConfig->second.doSaveRunEventNumbers_exclRejected_ << std::endl;
  //  std::cout << "   dqmDirectory_exclRejected = " << filterConfig->second.dqmDirectory_exclRejected_ << std::endl;
  //  std::cout << "  doSaveRunEventNumbers_passed_cumulative = " << filterConfig->second.doSaveRunEventNumbers_passed_cumulative_ << std::endl;
  //  std::cout << "   dqmDirectory_passed_cumulative = " << filterConfig->second.dqmDirectory_passed_cumulative_ << std::endl;
  //  std::cout << "  doSaveRunEventNumbers_rejected_cumulative = " << filterConfig->second.doSaveRunEventNumbers_rejected_cumulative_ << std::endl;
  //  std::cout << "   dqmDirectory_rejected_cumulative = " << filterConfig->second.dqmDirectory_rejected_cumulative_ << std::endl;
  //}
}

RunEventNumberService::~RunEventNumberService()
{
// nothing to be done yet...
}

void saveRunEventNumber(const std::string& dqmDirectory_store, edm::RunNumber_t runNumber, edm::EventNumber_t eventNumber, double eventWeight)
{
  if ( edm::Service<DQMStore>().isAvailable() ) {
    DQMStore& dqmStore = (*edm::Service<DQMStore>());
    //std::cout << "<store>:" << std::endl;
    if ( dqmDirectory_store != "" ) dqmStore.setCurrentFolder(dqmDirectory_store);
    std::ostringstream meName;
    meName << "r" << runNumber << "ev" << eventNumber;
    //std::cout << " meName = " << meName.str() << std::endl;
    std::ostringstream meValue;
    meValue << std::setprecision(3) << eventWeight;
    //std::cout << " meValue = " << meValue.str() << std::endl;
    dqmStore.bookString(meName.str(), meValue.str());
  } else {
    edm::LogError ("store") << " Failed to access dqmStore !!";
  }
}

void RunEventNumberService::update(const edm::RunNumber_t& runNumber, const edm::EventNumber_t& eventNumber, 
				   const filterResults_type& filterResults_cumulative, 
				   const filterResults_type& filterResults_individual, double eventWeight)
{
//--- check that configuration parameters contain no errors
  if ( cfgError_ ) {
    edm::LogError ("RunEventNumberService::update") << " Error in Configuration ParameterSet --> skipping !!";
    return;
  }  

//--- first pass through filterResults: 
//    count number of filters which passed/rejected the event
  unsigned numFiltersPassed_individual = 0;
  unsigned numFiltersRejected_individual = 0;
  for ( filterResults_type::const_iterator filterResult_individual = filterResults_individual.begin();
	filterResult_individual != filterResults_individual.end(); ++filterResult_individual ) {
    if (  filterResult_individual->second ) ++numFiltersPassed_individual;
    if ( !filterResult_individual->second ) ++numFiltersRejected_individual;
  }

//--- second pass through filterResults: 
//    save run & event numbers
  bool previousFiltersPassed = true;
  for ( filterResults_type::const_iterator filterResult_cumulative = filterResults_cumulative.begin();
	filterResult_cumulative != filterResults_cumulative.end(); ++filterResult_cumulative ) {
    const std::string& filterName = filterResult_cumulative->first;
    bool filterPassed_cumulative = filterResult_cumulative->second;

    filterResults_type::const_iterator filterResult_individual = filterResults_individual.end();
    for ( filterResults_type::const_iterator it = filterResults_individual.begin();
	  it != filterResults_individual.end(); ++it ) {
      if ( it->first == filterName ) {
	filterResult_individual = it;
	break;
      }
    }
    if ( filterResult_individual == filterResults_individual.end() ) {
      edm::LogError ("RunEventNumberService::update") << " Failed to find filterResult_individual for filterName = " << filterName
						      << " --> skipping !!";     
      continue;
    }

    bool filterPassed_individual = filterResult_individual->second;

    std::map<std::string, filterConfigEntry>::const_iterator it = filterConfigs_.find(filterName);
    if ( it == filterConfigs_.end() ) {
      edm::LogError ("RunEventNumberService::update") << " Failed to access filterConfigEntry for filterName = " << filterName
						      << " --> skipping !!";     
      continue;
    }
    
    const filterConfigEntry& entry = it->second;
    
    if ( entry.doSaveRunEventNumbers_passed_ && filterPassed_individual ) 
      saveRunEventNumber(entry.dqmDirectory_passed_, runNumber, eventNumber, eventWeight);
    if ( entry.doSaveRunEventNumbers_rejected_ && !filterPassed_individual ) 
      saveRunEventNumber(entry.dqmDirectory_rejected_, runNumber, eventNumber, eventWeight);
    if ( entry.doSaveRunEventNumbers_exclRejected_ && !filterPassed_individual && numFiltersRejected_individual == 1 ) 
      saveRunEventNumber(entry.dqmDirectory_exclRejected_, runNumber, eventNumber, eventWeight);
    if ( entry.doSaveRunEventNumbers_passed_cumulative_ && previousFiltersPassed && filterPassed_cumulative ) 
      saveRunEventNumber(entry.dqmDirectory_passed_cumulative_, runNumber, eventNumber, eventWeight);
    if ( entry.doSaveRunEventNumbers_rejected_cumulative_ && previousFiltersPassed && !filterPassed_cumulative ) 
      saveRunEventNumber(entry.dqmDirectory_rejected_cumulative_, runNumber, eventNumber, eventWeight);
    
    if ( !filterPassed_cumulative ) previousFiltersPassed = false;
  }
}
