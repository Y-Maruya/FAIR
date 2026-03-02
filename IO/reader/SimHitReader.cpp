#include "SimHitReader.hpp"

#include <stdexcept>

#include "TFile.h"
#include "TTree.h"

SimHitReader::SimHitReader(std::string filename,
                                   std::string treename)
{
    m_file.reset(TFile::Open(filename.c_str(), "READ"));
    if (!m_file || m_file->IsZombie()) {
        throw std::runtime_error("SimHitReader: cannot open file: " + filename);
    }

    m_tree = dynamic_cast<TTree*>(m_file->Get(treename.c_str()));
    if (!m_tree) {
        throw std::runtime_error("SimHitReader: cannot find TTree: " + treename);
    }

    m_entries = m_tree->GetEntries();
    LOG_INFO("SimHitReader: opened file '{}', found TTree '{}', total entries: {}", filename, treename, m_entries);
    // Read first branch to initialize the enviromental variables
    bind_branches_();
    m_tree->GetEntry(0);
}
SimHitReader::~SimHitReader() {
    cleanup_buffers_();
}
void SimHitReader::cleanup_buffers_() {
    delete vecHcalCellID; vecHcalCellID = nullptr;
    delete vecHcalStepsCell; vecHcalStepsCell = nullptr;
    delete vecHcalEdepCell; vecHcalEdepCell = nullptr;
    delete vecHcalVisibleEdepCell; vecHcalVisibleEdepCell = nullptr;
    delete vecHcalHitTimeCell; vecHcalHitTimeCell = nullptr;
    delete vecHcalToaCell; vecHcalToaCell = nullptr;

    delete vecTruth_pdgID; vecTruth_pdgID = nullptr;
    delete vecTruth_px; vecTruth_px = nullptr;
    delete vecTruth_py; vecTruth_py = nullptr;
    delete vecTruth_pz; vecTruth_pz = nullptr;
    delete vecTruth_x; vecTruth_x = nullptr;
    delete vecTruth_y; vecTruth_y = nullptr;
    delete vecTruth_z; vecTruth_z = nullptr;
    delete vecTruth_energy; vecTruth_energy = nullptr;
    delete vecTruth_vertexIndex; vecTruth_vertexIndex = nullptr;
    delete vecTruth_trackID; vecTruth_trackID = nullptr;

    // delete vecPlane_pdgID; vecPlane_pdgID = nullptr;
    // delete vecPlane_charge; vecPlane_charge = nullptr;
    // delete vecPlane_globalTime; vecPlane_globalTime = nullptr;
    // delete vecPlane_px; vecPlane_px = nullptr;
    // delete vecPlane_py; vecPlane_py = nullptr;
    // delete vecPlane_pz; vecPlane_pz = nullptr;
    // delete vecPlane_x; vecPlane_x = nullptr;
    // delete vecPlane_y; vecPlane_y = nullptr;
    // delete vecPlane_z; vecPlane_z = nullptr;
    // delete vecPlane_energy; vecPlane_energy = nullptr;
    // delete vecPlane_vertexIndex; vecPlane_vertexIndex = nullptr;
    // delete vecPlane_trackID; vecPlane_trackID = nullptr;
    // delete vecPlane_primary_Dmeson; vecPlane_primary_Dmeson = nullptr;
}
void SimHitReader::clear_vectors_() {
    if (vecHcalCellID) vecHcalCellID->clear();
    if (vecHcalStepsCell) vecHcalStepsCell->clear();
    if (vecHcalEdepCell) vecHcalEdepCell->clear();
    if (vecHcalVisibleEdepCell) vecHcalVisibleEdepCell->clear();
    if (vecHcalHitTimeCell) vecHcalHitTimeCell->clear();
    if (vecHcalToaCell) vecHcalToaCell->clear();

    if (vecTruth_pdgID) vecTruth_pdgID->clear();
    if (vecTruth_px) vecTruth_px->clear();
    if (vecTruth_py) vecTruth_py->clear();
    if (vecTruth_pz) vecTruth_pz->clear();
    if (vecTruth_x) vecTruth_x->clear();
    if (vecTruth_y) vecTruth_y->clear();
    if (vecTruth_z) vecTruth_z->clear();
    if (vecTruth_energy) vecTruth_energy->clear();
    if (vecTruth_vertexIndex) vecTruth_vertexIndex->clear();
    if (vecTruth_trackID) vecTruth_trackID->clear();

    // if (vecPlane_pdgID) vecPlane_pdgID->clear();
    // if (vecPlane_charge) vecPlane_charge->clear();
    // if (vecPlane_globalTime) vecPlane_globalTime->clear();
    // if (vecPlane_px) vecPlane_px->clear();
    // if (vecPlane_py) vecPlane_py->clear();
    // if (vecPlane_pz) vecPlane_pz->clear();
    // if (vecPlane_x) vecPlane_x->clear();
    // if (vecPlane_y) vecPlane_y->clear();
    // if (vecPlane_z) vecPlane_z->clear();
    // if (vecPlane_energy) vecPlane_energy->clear();
    // if (vecPlane_vertexIndex) vecPlane_vertexIndex->clear();
    // if (vecPlane_trackID) vecPlane_trackID->clear();
    // if (vecPlane_primary_Dmeson) vecPlane_primary_Dmeson->clear();
}
void SimHitReader::bind_branches_() {
    if (!vecHcalCellID) vecHcalCellID = new std::vector<int>();
    if (!vecHcalStepsCell) vecHcalStepsCell = new std::vector<int>();
    if (!vecHcalEdepCell) vecHcalEdepCell = new std::vector<double>();
    if (!vecHcalVisibleEdepCell) vecHcalVisibleEdepCell = new std::vector<double>();
    if (!vecHcalHitTimeCell) vecHcalHitTimeCell = new std::vector<double>();
    if (!vecHcalToaCell) vecHcalToaCell = new std::vector<double>();
    if (!vecTruth_pdgID) vecTruth_pdgID = new std::vector<int>();
    if (!vecTruth_px) vecTruth_px = new std::vector<double>();
    if (!vecTruth_py) vecTruth_py = new std::vector<double>();
    if (!vecTruth_pz) vecTruth_pz = new std::vector<double>();
    if (!vecTruth_x) vecTruth_x = new std::vector<double>();
    if (!vecTruth_y) vecTruth_y = new std::vector<double>();
    if (!vecTruth_z) vecTruth_z = new std::vector<double>();
    if (!vecTruth_energy) vecTruth_energy = new std::vector<double>();
    if (!vecTruth_vertexIndex) vecTruth_vertexIndex = new std::vector<int>();
    if (!vecTruth_trackID) vecTruth_trackID = new std::vector<int>();
    // if (!vecPlane_pdgID) vecPlane_pdgID = new std::vector<int>();
    // if (!vecPlane_charge) vecPlane_charge = new std::vector<double>();
    // if (!vecPlane_globalTime) vecPlane_globalTime = new std::vector<double>();
    // if (!vecPlane_px) vecPlane_px = new std::vector<double>();
    // if (!vecPlane_py) vecPlane_py = new std::vector<double>();
    // if (!vecPlane_pz) vecPlane_pz = new std::vector<double>();
    // if (!vecPlane_x) vecPlane_x = new std::vector<double>();
    // if (!vecPlane_y) vecPlane_y = new std::vector<double>();
    // if (!vecPlane_z) vecPlane_z = new std::vector<double>();
    // if (!vecPlane_energy) vecPlane_energy = new std::vector<double>();
    // if (!vecPlane_vertexIndex) vecPlane_vertexIndex = new std::vector<int>();
    // if (!vecPlane_trackID) vecPlane_trackID = new std::vector<int>();
    // if (!vecPlane_primary_Dmeson) vecPlane_primary_Dmeson = new std::vector<int>();
    m_tree->SetBranchStatus("*", 0); // disable all branches
    // enable only the branches we need

    // m_tree->SetBranchAddress("EvtID", &EvtID);
    m_tree->SetBranchStatus("ParticleEnergy", 1);
    m_tree->SetBranchAddress("ParticleEnergy", &ParticleEnergy);
    m_tree->SetBranchStatus("ftagNulabel", 1);
    m_tree->SetBranchAddress("ftagNulabel", &ftagNulabel);
    m_tree->SetBranchStatus("Interaction_x", 1);
    m_tree->SetBranchStatus("Interaction_y", 1);
    m_tree->SetBranchStatus("Interaction_z", 1);
    m_tree->SetBranchAddress("Interaction_x", &Interaction_x);
    m_tree->SetBranchAddress("Interaction_y", &Interaction_y);
    m_tree->SetBranchAddress("Interaction_z", &Interaction_z);

    if (m_tree->GetBranch("SecondaryEnergy"))
    {
        SecondaryBranchExist = true;
        m_tree->SetBranchStatus("SecondaryEnergy", 1);
        m_tree->SetBranchStatus("Secondarypdgid", 1);
        m_tree->SetBranchStatus("SecondaryMomentum_px", 1);
        m_tree->SetBranchStatus("SecondaryMomentum_py", 1);
        m_tree->SetBranchStatus("SecondaryMomentum_pz", 1);
        m_tree->SetBranchAddress("SecondaryEnergy", &SecondaryEnergy);
        m_tree->SetBranchAddress("Secondarypdgid", &Secondarypdgid);
        m_tree->SetBranchAddress("SecondaryMomentum_px", &SecondaryMomentum_px);
        m_tree->SetBranchAddress("SecondaryMomentum_py", &SecondaryMomentum_py);
        m_tree->SetBranchAddress("SecondaryMomentum_pz", &SecondaryMomentum_pz);
    }

    
    // m_tree->SetBranchAddress("CaloEdepSum", &CaloEdepSum);
    // m_tree->SetBranchAddress("CaloVisibleEdepSum", &CaloVisibleEdepSum);
    // m_tree->SetBranchAddress("EcalEdepSum", &EcalEdepSum);
    // m_tree->SetBranchAddress("EcalVisibleEdepSum", &EcalVisibleEdepSum);
    // m_tree->SetBranchAddress("HcalEdepSum", &HcalEdepSum);

    // m_tree->SetBranchAddress("EcalMaxEdepCell", &EcalMaxEdepCell);
    // m_tree->SetBranchAddress("HcalMaxEdepCell", &HcalMaxEdepCell);
    m_tree->SetBranchStatus("vecHcalCellID", 1);
    m_tree->SetBranchAddress("vecHcalCellID", &vecHcalCellID);
    // m_tree->SetBranchAddress("vecHcalStepsCell", &vecHcalStepsCell);
    // m_tree->SetBranchAddress("vecHcalEdepCell", &vecHcalEdepCell);
    m_tree->SetBranchStatus("vecHcalVisibleEdepCell", 1);
    m_tree->SetBranchAddress("vecHcalVisibleEdepCell", &vecHcalVisibleEdepCell);
    m_tree->SetBranchStatus("vecHcalHitTimeCell", 1);
    m_tree->SetBranchAddress("vecHcalHitTimeCell", &vecHcalHitTimeCell);
    m_tree->SetBranchStatus("vecHcalToaCell", 1);
    m_tree->SetBranchAddress("vecHcalToaCell", &vecHcalToaCell);
    if (m_tree->GetBranch("nstoredTruthParticles"))
    {
        vecTruthBranchExist = true;
        m_tree->SetBranchStatus("nstoredTruthParticles", 1);
        m_tree->SetBranchStatus("vecTruth_pdgID", 1);
        m_tree->SetBranchStatus("vecTruth_px", 1);
        m_tree->SetBranchStatus("vecTruth_py", 1);
        m_tree->SetBranchStatus("vecTruth_pz", 1);
        // m_tree->SetBranchStatus("vecTruth_x", 1);
        // m_tree->SetBranchStatus("vecTruth_y", 1);
        // m_tree->SetBranchStatus("vecTruth_z", 1);
        m_tree->SetBranchStatus("vecTruth_energy", 1);
        // m_tree->SetBranchStatus("vecTruth_vertexIndex", 1);
        // m_tree->SetBranchStatus("vecTruth_trackID", 1);
        m_tree->SetBranchAddress("nstoredTruthParticles", &nstoredTruthParticles);
        m_tree->SetBranchAddress("vecTruth_pdgID", &vecTruth_pdgID);
        m_tree->SetBranchAddress("vecTruth_px", &vecTruth_px);
        m_tree->SetBranchAddress("vecTruth_py", &vecTruth_py);
        m_tree->SetBranchAddress("vecTruth_pz", &vecTruth_pz);
        // m_tree->SetBranchAddress("vecTruth_x", &vecTruth_x);
        // m_tree->SetBranchAddress("vecTruth_y", &vecTruth_y);
        // m_tree->SetBranchAddress("vecTruth_z", &vecTruth_z);
        m_tree->SetBranchAddress("vecTruth_energy", &vecTruth_energy);
        // m_tree->SetBranchAddress("vecTruth_vertexIndex", &vecTruth_vertexIndex);
        // m_tree->SetBranchAddress("vecTruth_trackID", &vecTruth_trackID);
    }
    // m_tree->SetBranchAddress("nstoredPlaneParticles", &nstoredPlaneParticles);
    // m_tree->SetBranchAddress("vecPlane_pdgID", &vecPlane_pdgID);
    // m_tree->SetBranchAddress("vecPlane_charge", &vecPlane_charge);
    // m_tree->SetBranchAddress("vecPlane_globalTime", &vecPlane_globalTime);
    // m_tree->SetBranchAddress("vecPlane_px", &vecPlane_px);
    // m_tree->SetBranchAddress("vecPlane_py", &vecPlane_py);
    // m_tree->SetBranchAddress("vecPlane_pz", &vecPlane_pz);
    // m_tree->SetBranchAddress("vecPlane_x", &vecPlane_x);
    // m_tree->SetBranchAddress("vecPlane_y", &vecPlane_y);
    // m_tree->SetBranchAddress("vecPlane_z", &vecPlane_z);
    // m_tree->SetBranchAddress("vecPlane_energy", &vecPlane_energy);
    // m_tree->SetBranchAddress("vecPlane_vertexIndex", &vecPlane_vertexIndex);
    // m_tree->SetBranchAddress("vecPlane_trackID", &vecPlane_trackID);
    // m_tree->SetBranchAddress("vecPlane_primary_Dmeson", &vecPlane_primary_Dmeson);
}

bool SimHitReader::next(std::vector<AHCALSimHit>& out_hits, SimData& out_simdata) {
    if (m_entry >= m_entries) return false;
    clear_vectors_();
    m_tree->GetEntry(m_entry);
    LOG_DEBUG("SimHitReader: reading entry {}/{}", m_entry, m_entries);
    LOG_DEBUG("SimHitReader: ParticleEnergy={}, ftagNulabel={}, Interaction_x={}, Interaction_y={}, Interaction_z={}", 
              ParticleEnergy, ftagNulabel, Interaction_x, Interaction_y, Interaction_z);
    if (!vecHcalCellID || !vecHcalVisibleEdepCell) {
        LOG_ERROR("SimHitReader: missing essential branch data {}", 
                  vecHcalCellID ? "vecHcalVisibleEdepCell" : "vecHcalCellID");
        LOG_ERROR("SimHitReader: missing essential branch data");
        throw std::runtime_error("SimHitReader: missing branch data");
    }
    if (vecTruthBranchExist){
      out_simdata.injected_pdgId = vecTruth_pdgID->empty() ? 0 : vecTruth_pdgID->at(0);
      out_simdata.injected_energy = ParticleEnergy;
      if (ParticleEnergy != vecTruth_energy->at(0)) {
          LOG_WARN("SimHitReader: ParticleEnergy does not match vecTruth_energy[0]");
      }
      out_simdata.injected_time = 0; // not available in current TTree, set to -1 as default
      out_simdata.injected_z = -999999; // not available in current TTree, set to -1 as default
      out_simdata.injected_x = -999999; // not available in current TTree, set to -1 as default
      out_simdata.injected_y = -999999; // not available in current TTree, set to -1 as default
      out_simdata.injected_px = vecTruth_px->empty() ? 0 : vecTruth_px->at(0);
      out_simdata.injected_py = vecTruth_py->empty() ? 0 : vecTruth_py->at(0);
      out_simdata.injected_pz = vecTruth_pz->empty() ? 0 : vecTruth_pz->at(0);
    }else{
      out_simdata.injected_pdgId = 0;
      out_simdata.injected_energy = ParticleEnergy;
      out_simdata.injected_time = 0;
      out_simdata.injected_z = -999999;
      out_simdata.injected_x = -999999;
      out_simdata.injected_y = -999999;
      out_simdata.injected_px = -999999;
      out_simdata.injected_py = -999999;
      out_simdata.injected_pz = -999999;
    }
    out_simdata.isNeutrinoEvent = ftagNulabel >= 0 && ftagNulabel <= 3;
    out_simdata.ftagNulabel = ftagNulabel;
    out_simdata.isCC = (ftagNulabel == 0 || ftagNulabel == 1 || ftagNulabel == 2);
    out_simdata.Interaction_x = Interaction_x;
    out_simdata.Interaction_y = Interaction_y;
    out_simdata.Interaction_z = Interaction_z;
    if (SecondaryBranchExist){
      out_simdata.SecondaryEnergy = SecondaryEnergy;
      out_simdata.Secondarypdgid = Secondarypdgid;
      out_simdata.Secondary_px = SecondaryMomentum_px;
      out_simdata.Secondary_py = SecondaryMomentum_py;
      out_simdata.Secondary_pz = SecondaryMomentum_pz;
    }else{
      out_simdata.SecondaryEnergy = -999999;
      out_simdata.Secondarypdgid = 0;
      out_simdata.Secondary_px = -999999;
      out_simdata.Secondary_py = -999999;
      out_simdata.Secondary_pz = -999999;
    }
    out_hits.clear();
    out_hits.reserve(vecHcalCellID->size());
    for (size_t i = 0; i < vecHcalCellID->size(); ++i) {
        AHCALSimHit hit;
        hit.cellID = cellid_conversion(vecHcalCellID->at(i));
        hit.Edep = vecHcalVisibleEdepCell->at(i); // after birks saturation, in MeV
        hit.Nmip = hit.Edep / AHCALGeometry::MIPEnergy; // convert energy deposition to number of MIPs, using the most probable energy deposition of a MIP in AHCAL as the conversion factor
        hit.HitTime = vecHcalHitTimeCell->at(i); // the time of maximum energy deposition in the cell, in ns
        hit.TimeOfArrival = vecHcalToaCell->at(i); //Time of arrival, in ns
        out_hits.push_back(hit);
    }
    ++m_entry;
    return true;
}

