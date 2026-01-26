#ifndef ADC_TO_ENERGY_READ_TTree_HPP
#define ADC_TO_ENERGY_READ_TTree_HPP
#include "common/EventStore.hpp"
#include "common/IAlg.hpp"
#include <string>
#include <vector>
#include <memory>
#include <TFile.h>
#include <TTree.h>
#include <map>
#include <tuple>

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

    private:
        std::string m_in_rawhit_key;
        std::string m_out_recohit_key;
        std::unique_ptr<TFile> m_in_file;
        int file_cellid_version = 1; 
        TTree* m_in_tree = nullptr;
        std::string cut_string;
        std::map<int, double> mip_map; // cellID to MPV
        std::map<int, double> hg_ped_map; // cellID to pedestal
        std::map<int, double> lg_ped_map; // cellID to pedestal
        std::map<int, double> gainratio_map; // cellID to gain ratio
        std::map<int, int> gainplat_map; // cellID to gain plat
        int cellid_conversion(int input_cellid);
        AdcToEnergyReadTTreeAlgCfg m_cfg;
    };
}
#endif // ADC_TO_ENERGY_READ_MIP_TTree_HPP