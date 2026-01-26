#include "RootRawHitReader.hpp"

#include <stdexcept>

#include "TFile.h"
#include "TTree.h"

RootRawHitReader::RootRawHitReader(std::string filename,
                                   std::string treename)
{
  m_file.reset(TFile::Open(filename.c_str(), "READ"));
  if (!m_file || m_file->IsZombie()) {
    throw std::runtime_error("RootRawHitReader: cannot open file: " + filename);
  }

  m_tree = dynamic_cast<TTree*>(m_file->Get(treename.c_str()));
  if (!m_tree) {
    throw std::runtime_error("RootRawHitReader: cannot find TTree: " + treename);
  }

  m_entries = m_tree->GetEntries();
  bind_branches_();
  // Read first branch to initialize the enviromental variables
  m_tree->GetEntry(0);
}

void RootRawHitReader::bind_branches_() {
    //TLU branches
    m_tree->SetBranchAddress("Timestamp", &b_timestamp);
    m_tree->SetBranchAddress("BCID_TLU", &b_bc_id_tlu);
    m_tree->SetBranchAddress("Inputs", &b_inputs);
    m_tree->SetBranchAddress("FineTimestamps", &b_fine_timestamps);
    //AHCAL branches
    m_tree->SetBranchAddress("Run_No", &b_runNo);
    m_tree->SetBranchAddress("CycleID", &b_cycleID);
    if (m_tree->GetBranch("TriggerID_TLU"))
    {
      m_tree->SetBranchAddress("TriggerID_TLU", &b_triggerID);
    }
    else{
      m_tree->SetBranchAddress("TriggerID", &b_triggerID);
    }
    // m_tree->SetBranchAddress("TriggerID_TLU", &b_triggerID);
    m_tree->SetBranchAddress("Event_Time", &b_Event_Time);
    m_tree->SetBranchAddress("CellID", &b_cellID);
    m_tree->SetBranchAddress("BCID", &b_bcid);
    m_tree->SetBranchAddress("HitTag", &b_hitTag);
    m_tree->SetBranchAddress("HG_Charge", &b_hg);
    m_tree->SetBranchAddress("LG_Charge", &b_lg);
}

bool RootRawHitReader::next(std::vector<AHCALRawHit>& out_hits) {
  ++m_entry;
  if (m_entry >= m_entries) return false;

  m_tree->GetEntry(m_entry);

  if (!b_cellID || !b_hg || !b_lg) {
    throw std::runtime_error("RootRawHitReader: branch buffers not bound");
  }
  if (b_cellID->size() != b_hg->size() || b_cellID->size() != b_lg->size()) {
    throw std::runtime_error("RootRawHitReader: branch vector size mismatch");
  }

  out_hits.clear();
  out_hits.reserve(b_cellID->size());

  for (size_t i = 0; i < b_cellID->size(); ++i) {
    AHCALRawHit h;
    h.index = static_cast<int>(i);
    h.cellID = (*b_cellID)[i];
    h.hg_adc = (*b_hg)[i];
    h.lg_adc = (*b_lg)[i];
    h.bcid   = (b_bcid ? (*b_bcid)[i] : 0);
    h.hittag = (b_hitTag ? (*b_hitTag)[i] : 0);
    out_hits.push_back(h);
  }
  return true;
}
