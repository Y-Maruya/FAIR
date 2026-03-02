#pragma once 
#include "common/edm/SimHit.hpp"
#include "common/edm/SimData.hpp"

class SimHitReader{
public:
    SimHitReader(std::string filename,
                    std::string treename = "treeEvt");
    ~SimHitReader();

    bool next(std::vector<AHCALSimHit>& out_hits, SimData& out_simdata);

    long long entry() const { return m_entry; }
    long long entries() const { return m_entries; }

private:
    std::unique_ptr<TFile> m_file;
    TTree* m_tree = nullptr;
    long long m_entries = 0;
    long long m_entry = 0;
        
    
    int EvtID = 0;
    double ParticleEnergy = 0;
    int ftagNulabel = 0;
    double Interaction_x = 0;
    double Interaction_y = 0;
    double Interaction_z = 0;

    bool SecondaryBranchExist = false;
    double SecondaryEnergy = 0;
    int Secondarypdgid = 0;
    double SecondaryMomentum_px = 0;
    double SecondaryMomentum_py = 0;
    double SecondaryMomentum_pz = 0;
    
    double CaloEdepSum = 0;
    double CaloVisibleEdepSum = 0;
    double HcalEdepSum = 0;
    double HcalVisibleEdepSum = 0;
    double HcalMaxEdepCell = 0;
    std::vector<int>* vecHcalCellID = nullptr;
    std::vector<int>* vecHcalStepsCell = nullptr;
    std::vector<double>* vecHcalEdepCell = nullptr;
    std::vector<double>* vecHcalVisibleEdepCell = nullptr;
    std::vector<double>* vecHcalHitTimeCell = nullptr;
    std::vector<double>* vecHcalToaCell = nullptr;
    bool vecTruthBranchExist;
    int nstoredTruthParticles;
    std::vector<int>* vecTruth_pdgID = nullptr;
    std::vector<double>* vecTruth_px = nullptr;
    std::vector<double>* vecTruth_py = nullptr;
    std::vector<double>* vecTruth_pz = nullptr;
    std::vector<double>* vecTruth_x = nullptr;
    std::vector<double>* vecTruth_y = nullptr;
    std::vector<double>* vecTruth_z = nullptr;
    std::vector<double>* vecTruth_energy = nullptr;
    std::vector<int>* vecTruth_vertexIndex = nullptr;
    std::vector<int>* vecTruth_trackID = nullptr;
    // int nstoredPlaneParticles;
    // std::vector<int>* vecPlane_pdgID = nullptr;
    // std::vector<double>* vecPlane_charge = nullptr;
    // std::vector<double>* vecPlane_globalTime = nullptr;
    // std::vector<double>* vecPlane_px = nullptr;
    // std::vector<double>* vecPlane_py = nullptr;
    // std::vector<double>* vecPlane_pz = nullptr;
    // std::vector<double>* vecPlane_x = nullptr;
    // std::vector<double>* vecPlane_y = nullptr;
    // std::vector<double>* vecPlane_z = nullptr;
    // std::vector<double>* vecPlane_energy = nullptr;
    // std::vector<int>* vecPlane_vertexIndex = nullptr;
    // std::vector<int>* vecPlane_trackID = nullptr;
    // std::vector<int>* vecPlane_primary_Dmeson = nullptr;

    void bind_branches_();
    void cleanup_buffers_();
    void clear_vectors_();
    int cellid_conversion(int sim_cellid){
        int ID_X = sim_cellid % 100 - 1;
        int ID_Y = (sim_cellid % 10000) / 100 - 1;
        int ID_Z = sim_cellid / 10000 -1;
        int chipID = ID_X / 6 + ID_Y  / 6 * 3;
        int X_inchip = 5 - ID_X % 6;
        int Y_inchip = ID_Y % 6;
        int ChannelID = AHCALGeometry::_Channel[Y_inchip][X_inchip];
        if (chipID%3 != 0){
            if (ChannelID == 2) ChannelID = 0;
            else if (ChannelID == 0) ChannelID = 2;
            if (ChannelID == 33) ChannelID = 35;
            else if (ChannelID == 35) ChannelID = 33;
        }
        int new_cellid = ID_Z * 100000 + chipID * 10000 + ChannelID;
        return new_cellid;
    }
};