#pragma once
#include "calibration/RefValues.hpp"
#include "common/geometry/Geometry.hpp"
#include "common/EventStore.hpp"
#include "common/IAlg.hpp"
#include <string>
#include <vector>
#include <memory>
#include <TFile.h>
#include <TTree.h>
#include <cstddef>

namespace AHCALRecoAlg {
    struct AdcToEnergyReadTTreeAlgCfg{
        std::string in_rawhit_key = "RawHits";
        std::string out_recohit_key = "RecoHits";
        int mip_cellid_version = 1;
        int ped_cellid_version = 1;
        int dac_cellid_version = 1;
        std::string mip_file = "";
        std::string ped_file = "";
        std::string dac_file = "";
        std::string mip_cut_string = "";
        std::string ped_cut_string = "";
        std::string dac_cut_string = "";
    };
    class AdcToEnergyReadTTreeAlg final : public IAlg {
    public:
        AdcToEnergyReadTTreeAlg(RunContext& ctx, std::string name)
            : IAlg(ctx, name){ }
        ~AdcToEnergyReadTTreeAlg();
        bool initialize_mip();
        bool initialize_ped();
        bool initialize_dac();
        void execute(EventStore& evt) override;
        void parse_cfg(const YAML::Node& n);
        void initialize() override {
            const std::size_t n = total_channel_count();
            mip_values_.assign(n, AHCALRefValues::ref_MIP);
            hg_ped_values_.assign(n, AHCALRefValues::ref_ped_highgain);
            lg_ped_values_.assign(n, AHCALRefValues::ref_ped_lowgain);
            gainratio_values_.assign(n, AHCALRefValues::ref_gain_ratio);
            gainplat_values_.assign(n, AHCALRefValues::lowgain_plat);
            initialize_mip();
            initialize_ped();
            initialize_dac();
        }
    private:
        std::string m_in_rawhit_key;
        std::string m_out_recohit_key;
        std::unique_ptr<TFile> m_in_file;
        int file_cellid_version = 1; 
        TTree* m_in_tree = nullptr;
        std::vector<double> mip_values_;
        std::vector<double> hg_ped_values_;
        std::vector<double> lg_ped_values_;
        std::vector<double> gainratio_values_;
        std::vector<int> gainplat_values_;
        int cellid_conversion(int input_cellid);
        static std::size_t total_channel_count();
        static std::size_t channel_index_from_cellid(int cellid);
        AdcToEnergyReadTTreeAlgCfg m_cfg;
    };
}
