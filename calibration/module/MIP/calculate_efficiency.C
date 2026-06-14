/**
 * Calculate per-channel efficiency directly from a pre-fit MIP ROOT file.
 *
 * efficiency = hMIP_<cellid>::GetEntries()
 *              / Ntracks_pass_through_channel_<cellid>
 *
 * Usage:
 *   root -b -q 'calculate_efficiency.C("mip_nofit.root", "mip_efficiency.root")'
 */

#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>
#include <TKey.h>
#include <TParameter.h>
#include <TTree.h>

#include <iostream>
#include <string>
#include <unordered_map>

struct EfficiencyInput {
    std::unordered_map<int, Long64_t> entries;
    std::unordered_map<int, int> ntrack_pass_through;
};

void collectEfficiencyInput(TDirectory* dir, EfficiencyInput& input) {
    if (!dir) return;

    TIter next(dir->GetListOfKeys());
    TKey* key = nullptr;
    while ((key = static_cast<TKey*>(next()))) {
        TObject* obj = dir->Get(key->GetName());
        if (!obj) continue;

        if (obj->InheritsFrom(TDirectory::Class())) {
            collectEfficiencyInput(static_cast<TDirectory*>(obj), input);
        } else if (obj->InheritsFrom(TH1::Class())) {
            const std::string name = obj->GetName();
            const std::string prefix = "hMIP_";
            if (name.find(prefix) == 0) {
                const int cellid = std::stoi(name.substr(prefix.size()));
                input.entries[cellid] =
                    static_cast<TH1*>(obj)->GetEntries();
            }
        } else if (obj->InheritsFrom(TParameter<int>::Class())) {
            const std::string name = obj->GetName();
            const std::string prefix = "Ntracks_pass_through_channel_";
            if (name.find(prefix) == 0) {
                const int cellid = std::stoi(name.substr(prefix.size()));
                input.ntrack_pass_through[cellid] =
                    static_cast<TParameter<int>*>(obj)->GetVal();
            }
        }

    }
}

void calculate_efficiency(const char* input_file = "/afs/cern.ch/user/y/ymaruya/private/FASERlink/AHCAL/FAIR/out_calibration/mip_calib/21987-22133/mip_neighborcheck_nofit.root",
                          const char* output_file = "/afs/cern.ch/user/y/ymaruya/private/FASERlink/AHCAL/FAIR/out_calibration/mip_calib/21987-22133/mip_neighborcheck_efficiency.root") {
    TFile* input_file_ptr = TFile::Open(input_file, "READ");
    if (!input_file_ptr || input_file_ptr->IsZombie()) {
        std::cerr << "Error: cannot open input file: " << input_file
                  << std::endl;
        delete input_file_ptr;
        return;
    }

    EfficiencyInput input;
    TDirectory* mip_dir = input_file_ptr->GetDirectory("MIP");
    collectEfficiencyInput(mip_dir ? mip_dir : input_file_ptr, input);

    TFile* output_file_ptr = TFile::Open(output_file, "RECREATE");
    if (!output_file_ptr || output_file_ptr->IsZombie()) {
        std::cerr << "Error: cannot create output file: " << output_file
                  << std::endl;
        delete output_file_ptr;
        input_file_ptr->Close();
        delete input_file_ptr;
        return;
    }

    int cellid = -1;
    int layer = -1;
    int chip = -1;
    int channel = -1;
    Long64_t entries = 0;
    int ntrack_pass_through_channel = -1;
    double efficiency = -1.0;

    TTree* tree =
        new TTree("mip_efficiency", "MIP efficiency before fitting");
    tree->Branch("cellid", &cellid);
    tree->Branch("layer", &layer);
    tree->Branch("chip", &chip);
    tree->Branch("channel", &channel);
    tree->Branch("entries", &entries);
    tree->Branch("ntrack_pass_through_channel",
                 &ntrack_pass_through_channel);
    tree->Branch("efficiency", &efficiency);

    int missing_ntrack_count = 0;
    for (const auto& item : input.entries) {
        cellid = item.first;
        layer = cellid / 100000;
        chip = (cellid / 10000) % 10;
        channel = cellid % 10000;
        entries = item.second;

        const auto ntrack_it = input.ntrack_pass_through.find(cellid);
        if (ntrack_it != input.ntrack_pass_through.end()) {
            ntrack_pass_through_channel = ntrack_it->second;
        } else {
            ntrack_pass_through_channel = -1;
        }

        if (ntrack_pass_through_channel > 0) {
            efficiency =
                static_cast<double>(entries) / ntrack_pass_through_channel;
        } else {
            efficiency = -1.0;
            ++missing_ntrack_count;
        }

        tree->Fill();
    }

    output_file_ptr->cd();
    tree->Write();
    delete tree;
    output_file_ptr->Close();
    delete output_file_ptr;

    input_file_ptr->Close();
    delete input_file_ptr;

    std::cout << "Wrote efficiencies for " << input.entries.size()
              << " channels to " << output_file << std::endl;
    std::cout << "Channels without a positive NtrackPassThrough value: "
              << missing_ntrack_count << std::endl;
}
