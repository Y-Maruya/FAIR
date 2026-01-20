#ifndef RootRawHitReader_HPP
#define RootRawHitReader_HPP
#include "common/Hit.hpp"

#include <memory>
#include <string>
#include <vector>

class TFile;
class TTree;

class RootRawHitReader {
public:
    RootRawHitReader(std::string filename,
                    std::string treename = "Raw_Hit");

    // returns false when no more events
    bool next(std::vector<AHCALRawHit>& out_hits);

    // optional
    long long entry() const { return m_entry; }
    long long entries() const { return m_entries; }

private:
    std::unique_ptr<TFile> m_file;
    TTree* m_tree = nullptr;

    long long m_entry = -1;
    long long m_entries = 0;
    // branch buffers (simple branches) 
    int                     b_runNo;
    int                     b_triggerID;
    int                     b_cycleID;
    int            b_Event_Time;
    int b_timestamp;
    int b_bc_id_tlu;
    std::vector<int>* b_inputs = nullptr;
    std::vector<int>* b_fine_timestamps = nullptr;
  // branch buffers (vector branches)
    std::vector<int>*      b_cellID = nullptr;
    std::vector<unsigned short>* b_hg = nullptr;
    std::vector<unsigned short>* b_lg = nullptr;
    std::vector<unsigned short>* b_bcid = nullptr; 
    std::vector<unsigned short>* b_hitTag = nullptr; // erase no meaning data

  void bind_branches_();
};

#endif // RootRawHitReader_HPP