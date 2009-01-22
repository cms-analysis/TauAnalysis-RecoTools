#ifndef TauAnalysis_RecoTools_RunEventNumberService_h  
#define TauAnalysis_RecoTools_RunEventNumberService_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Provenance/interface/RunID.h"
#include "DataFormats/Provenance/interface/EventID.h"

#include <string>
#include <map>

class RunEventNumberService
{
  struct filterConfigEntry
  {
    bool doSaveRunEventNumbers_passed_;
    std::string dqmDirectory_passed_;
    bool doSaveRunEventNumbers_rejected_;
    std::string dqmDirectory_rejected_;
    bool doSaveRunEventNumbers_exclRejected_;
    std::string dqmDirectory_exclRejected_;
    bool doSaveRunEventNumbers_passed_cumulative_;
    std::string dqmDirectory_passed_cumulative_;
    bool doSaveRunEventNumbers_rejected_cumulative_;
    std::string dqmDirectory_rejected_cumulative_;
  };

 public: 
  explicit RunEventNumberService(const edm::ParameterSet&);
  ~RunEventNumberService();

  typedef std::vector<std::pair<std::string, bool> > filterResults_type;
  void update(const edm::RunNumber_t&, const edm::EventNumber_t&, const filterResults_type&, const filterResults_type&, double);

 private:
  std::string name_;

  std::string dqmDirectory_store_;

  std::map<std::string, filterConfigEntry> filterConfigs_;

  int cfgError_;
};

#endif  


