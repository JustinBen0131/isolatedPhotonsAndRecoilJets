#ifndef RJ_DATA_SOURCE_ENTRY_CONTRACT_H
#define RJ_DATA_SOURCE_ENTRY_CONTRACT_H
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace rj_source_entry
{
struct Identity { std::int64_t entry, physical; int run; };
class Accounting
{
 public:
  Accounting(std::int64_t begin, std::int64_t total, int run)
    : begin_(begin), total_(total), run_(run)
  {
    require(begin >= 0 && total > 0 && run > 0 &&
            begin <= std::numeric_limits<std::int64_t>::max()-total,
            "invalid source-entry interval");
  }
  void input(std::int64_t jet, std::int64_t calo, int run, std::int64_t physical)
  {
    require(!pending_, "previous input has no classified outcome");
    if (!(attempted_ < total_ && jet == calo && jet == begin_+attempted_))
      throw std::runtime_error("paired input cursor mismatch, skip, repeat, or overrun: jet="+
        std::to_string(jet)+" calo="+std::to_string(calo)+" expected="+
        std::to_string(begin_+attempted_));
    require(run == run_ && physical >= 0, "original EventHeader identity absent or different");
    current_ = {jet,physical,run};
    ++attempted_; pending_=true; skim_seen_=false; writer_seen_=false;
  }
  void skim(bool rejected)
  {
    require(pending_ && !skim_seen_, "skimmer outcome without unique current input");
    skim_seen_=true; ++skimmed_inputs_;
    if (rejected) { rejected_.push_back(current_); pending_=false; }
  }
  std::int64_t offset(int run, std::int64_t physical, bool valid)
  {
    require(pending_ && skim_seen_, "retained event lacks input/skimmer witness");
    require(valid && run==current_.run && physical==current_.physical,
            "retained EventHeader differs from original input");
    writer_seen_=true;
    return current_.entry-begin_; // Runtime adds the configured source begin once.
  }
  std::int64_t retained(int run, std::int64_t physical, bool valid)
  {
    const auto relative = offset(run,physical,valid);
    ++retained_; pending_=false;
    return relative;
  }
  void finishEvent()
  {
    if (!pending_) return; // A documented skimmer rejection is already classified.
    require(skim_seen_ && writer_seen_, "input never reached the replay writer");
    ++retained_; pending_=false;
  }
  void complete() const
  {
    require(!pending_ && attempted_==total_ && skimmed_inputs_==attempted_ &&
            retained_+static_cast<std::int64_t>(rejected_.size())==attempted_,
            "input/retained/documented-rejection conservation failed");
  }
  std::int64_t begin() const { return begin_; }
  std::int64_t attempted() const { return attempted_; }
  std::int64_t retained() const { return retained_; }
  const std::vector<Identity>& rejected() const { return rejected_; }
 private:
  static void require(bool condition, const char* reason)
  { if (!condition) throw std::runtime_error(reason); }
  std::int64_t begin_,total_; int run_;
  std::int64_t attempted_=0,retained_=0,skimmed_inputs_=0;
  bool pending_=false,skim_seen_=false,writer_seen_=false;
  Identity current_{-1,-1,-1};
  std::vector<Identity> rejected_;
};
}
#endif
