#include "InputEffAlg.hpp"
#include "common/AHCALGeometry.hpp"
#include "common/Logger.hpp"
#include "common/config/YAMLUtil.hpp"
#include "common/AlgRegistry.hpp"
#include "common/edm/SimpleFittedTrack.hpp"
#include "common/edm/Track.hpp"

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TEfficiency.h>
#include <TDirectory.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLatex.h>
#include <TLine.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TLegend.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <fstream>

AHCAL_REGISTER_ALG(AHCALRecoAlg::InputEffAlg, "InputEffAlg")
const double XYMIN = -500;
const double XYMAX = +500;
constexpr int NBIN_XY = 100;
namespace AHCALRecoAlg{

    struct InputEffAlg::Impl {
        explicit Impl(InputEffAlgCfg cfg) : cfg_(std::move(cfg)) {}
    
        void fill_inputs(const std::vector<int> inputs) {
            h_input_sum->Fill(std::accumulate(inputs.begin(), inputs.end(), 0));
        }
        void fill_single(int input_index, int nHits, bool is_has_hit) {
            h_full->Fill(input_index);
            if (is_has_hit) {
                h_passed->Fill(input_index);
                h_nHits_passed->Fill(nHits);
            } else {
                h_nHits_failed->Fill(nHits);
            }
        }
        void fill_single_input1(int input_index, int nHits, unsigned long long event_number, const RunContext& ctx, bool is_has_hit, std::ofstream* out_file) {
            h_inputsumeqq1->Fill(input_index);
            h_nHits_sum_inputs_eq_1->Fill(nHits);
            if (cfg_.write_to_txt && out_file) {
                (*out_file) << ctx.config.runNumber << " " << ctx.config.poolIndex << " " << event_number << " " << input_index << " " << is_has_hit << " " << nHits << "\n";
            }
            if (nHits > 20){
                LOG_WARN("InputEffAlg::Impl::fill_single_input1: event {} has {} hits, which is unusually high. Check if this is expected.", event_number, nHits);
                std::cout << ctx.config.runNumber << " " << ctx.config.poolIndex << " " << event_number << " " << input_index << " " << nHits << "\n";
            }
        }
        void fill_difBCID(int TLUBCID, int HitsBCID, bool is_passed) {
            int difBCID = TLUBCID - HitsBCID;
            if (is_passed) {
                h_difBCID_passed->Fill(difBCID);
                h2_BCID_2D->Fill(TLUBCID, HitsBCID);
                h_clock_passed->Fill((TLUBCID % 8));
            } else {
                h_difBCID_failed->Fill(difBCID);
                h_clock_failed->Fill((TLUBCID % 8));
            }
        }
        void fill_loss_category(std::vector<bool> lost_inputs) {
            bool lost_all = false;
            bool lost_1 = false;
            bool lost_2 = false;
            bool lost_3 = false;
            bool nolost = false;
            int num = lost_inputs.size();
            int lost_num = std::count(lost_inputs.begin(), lost_inputs.end(), true);
            if (lost_num == num) {
                lost_all = true;
                lost_1 = false;
                lost_2 = false;
                lost_3 = false;
                nolost = false;
            } else if (lost_num == 0) {
                lost_all = false;
                lost_1 = false;
                lost_2 = false;
                lost_3 = false;
                nolost = true;
            } else if (lost_num == num-1) {
                lost_all = false;
                lost_1 = true;
                lost_2 = false;
                lost_3 = false;
                nolost = false;
             } else if (lost_num == num-2) {
                lost_all = false;
                lost_1 = false;
                lost_2 = true;
                lost_3 = false;
                nolost = false;
            } else if (lost_num == num-3) {
                lost_all = false;
                lost_1 = false;
                lost_2 = false;
                lost_3 = true;
                nolost = false;
            } else {
                LOG_WARN("InputEffAlg::Impl::fill_loss_category: unexpected number of lost inputs: {}, which is larger than 3. Check if this is expected.", lost_num);
            }
            if (lost_all) {
                h_loss->Fill(0);
            } else if (lost_1) {
                h_loss->Fill(1);
            } else if (lost_2) {
                h_loss->Fill(2);
            } else if (lost_3) {
                h_loss->Fill(3);
            } else if (nolost) {
                h_loss->Fill(4); // nolost
            }
        }

        void fill_pattern( int layer1, int layer2, bool isdoublelost, bool issinglelost, int lost_layer, std::vector<int> fineTimeStamps) {
            if (isdoublelost){
                h_double_lost[layer1]->Fill(layer2);
                h_double_lost[layer2]->Fill(layer1);
                h_timing_diff_double_lost[layer1*4+layer2]->Fill(fineTimeStamps[layer2] - fineTimeStamps[layer1]);
                h_timing_diff_double_lost[layer2*4+layer1]->Fill(fineTimeStamps[layer1] - fineTimeStamps[layer2]);
            } else if (issinglelost) {
                int other_layer = -1;
                for (int i = 0; i < cfg_.ntrigger_layer; ++i) {
                    if (i != lost_layer && (i == layer1 || i == layer2)) {
                        other_layer = i;
                        break;
                    }
                }
                h_single_lost[other_layer]->Fill(lost_layer);
                h_timing_diff_single_lost[other_layer*4+lost_layer]->Fill(fineTimeStamps[lost_layer] - fineTimeStamps[other_layer]);
            } else {
                h_timing_diff_no_lost[layer1*4+layer2]->Fill(fineTimeStamps[layer2] - fineTimeStamps[layer1]);
                h_timing_diff_no_lost[layer2*4+layer1]->Fill(fineTimeStamps[layer1] - fineTimeStamps[layer2]);
             }
        }
        void write() {
            if (written_) return;
            written_ = true;

            auto fout = std::unique_ptr<TFile>(TFile::Open(cfg_.out_filename.c_str(), "RECREATE"));
            if (!fout || fout->IsZombie()) {
                LOG_ERROR("InputEffAlg: cannot create output file: {}", cfg_.out_filename);
                return;
            }
            LOG_INFO("InputEffAlg: writing output to {}", cfg_.out_filename);
            fout->cd();
            h_full->Write();
            h_passed->Write();
            h_input_sum->Write();
            h_inputsumeqq1->Write();
            h_nHits_sum_inputs_eq_1->Write();
            h_difBCID_passed->Write();
            h_difBCID_failed->Write();
            h2_BCID_2D->Write();
            h_clock_passed->Write();
            h_clock_failed->Write();
            h_nHits_passed->Write();
            h_nHits_failed->Write();
            h_deltaBCID->Write();
            h_deltaCycleID->Write();
            h_deltaTimestamp->Write();
            h_deltaTriggerID->Write();
            h_cycle_duration->Write();
            for (int i = 0; i < cfg_.ntrigger_layer; ++i) {
                for (int j = 1; j < cfg_.ntrigger_layer; ++j) {
                    h_timing_diff_double_lost[i*4+j]->Write();
                    h_timing_diff_single_lost[i*4+j]->Write();
                    h_timing_diff_no_lost[i*4+j]->Write();
                }
             }
            h_loss->GetXaxis()->SetBinLabel(1, "All Lost");
            h_loss->GetXaxis()->SetBinLabel(2, "1 Lost");
            h_loss->GetXaxis()->SetBinLabel(3, "2 Lost");
            h_loss->GetXaxis()->SetBinLabel(4, "3 Lost"); 
            h_loss->GetXaxis()->SetBinLabel(5, "No Loss");
            h_loss->Write();
            for (int i = 0; i < cfg_.ntrigger_layer; ++i) {
                h_double_lost[i]->Write();
                h_single_lost[i]->Write();
            }

            std::unique_ptr<TEfficiency> eff = std::make_unique<TEfficiency>(*h_passed, *h_full);
            eff->SetName("input_efficiency");
            eff->SetTitle("HitTag Efficiency for Inputs;Input;Efficiency");
            eff->Draw("AP");
            eff->Write();
            gPad->Update();
            auto graph = eff->GetPaintedGraph();
            if (graph) {
                graph->GetXaxis()->SetBinLabel(1, Form("Input0 L%d", cfg_.trigger_layer[0]));
                graph->GetXaxis()->SetBinLabel(2, Form("Input1 L%d", cfg_.trigger_layer[1]));
                graph->GetXaxis()->SetBinLabel(3, Form("Input2 L%d", cfg_.trigger_layer[2]));
                graph->GetXaxis()->SetBinLabel(4, Form("Input3 L%d", cfg_.trigger_layer[3]));
                graph->GetXaxis()->SetBinLabel(5, "Input4 Veto1");
                graph->GetXaxis()->SetBinLabel(6, "Input5 Veto2");
            }

            std::vector<std::unique_ptr<TEfficiency>> effs;
            for (int i = 0; i < cfg_.ntrigger_layer; ++i) {
                std::unique_ptr<TEfficiency> eff_single = std::make_unique<TEfficiency>(*h_single_lost[i], *h_double_lost[i]);
                eff_single->SetName(Form("eff_single_lost_input%d", i));
                eff_single->SetTitle(Form("Accidental Hit Efficiency for Input %d;Input;Efficiency", i));
                effs.push_back(std::move(eff_single));
            }

            if (cfg_.write_to_png) {
                if (gSystem->AccessPathName(cfg_.out_png_dir.c_str())) {
                    if (gSystem->mkdir(cfg_.out_png_dir.c_str(), true) != 0) {
                        LOG_ERROR("InputEffAlg: cannot create output PNG directory: {}", cfg_.out_png_dir);
                        return;
                    }
                }
                TCanvas c("c_input_efficiency", "HitTag Efficiency for Inputs", 800, 600);
                TH1D* h_dummy = new TH1D("h_dummy", eff->GetTitle(), 6, -0.5, 5.5);
                h_dummy->GetXaxis()->SetBinLabel(1, Form("Input0 L%d", cfg_.trigger_layer[0]));
                h_dummy->GetXaxis()->SetBinLabel(2, Form("Input1 L%d", cfg_.trigger_layer[1]));
                h_dummy->GetXaxis()->SetBinLabel(3, Form("Input2 L%d", cfg_.trigger_layer[2]));
                h_dummy->GetXaxis()->SetBinLabel(4, Form("Input3 L%d", cfg_.trigger_layer[3]));
                h_dummy->GetXaxis()->SetBinLabel(5, "Input4 Veto1");
                h_dummy->GetXaxis()->SetBinLabel(6, "Input5 Veto2");
                h_dummy->Draw();
                eff->Draw("Psame");
                c.SaveAs(Form("%s/input_efficiency.png", cfg_.out_png_dir.c_str()));
                h_full->Draw();
                h_full->GetXaxis()->SetBinLabel(1, Form("Input0 L%d", cfg_.trigger_layer[0]));
                h_full->GetXaxis()->SetBinLabel(2, Form("Input1 L%d", cfg_.trigger_layer[1]));
                h_full->GetXaxis()->SetBinLabel(3, Form("Input2 L%d", cfg_.trigger_layer[2]));
                h_full->GetXaxis()->SetBinLabel(4, Form("Input3 L%d", cfg_.trigger_layer[3]));
                h_full->GetXaxis()->SetBinLabel(5, "Input4 Veto1");
                h_full->GetXaxis()->SetBinLabel(6, "Input5 Veto2");
                h_full->Draw();
                c.SaveAs(Form("%s/input_full.png", cfg_.out_png_dir.c_str()));
                h_input_sum->Draw();
                c.SaveAs(Form("%s/input_sum.png", cfg_.out_png_dir.c_str()));
                h_inputsumeqq1->Draw();
                c.SaveAs(Form("%s/input_sum_eq_1.png", cfg_.out_png_dir.c_str()));
                h_nHits_sum_inputs_eq_1->Draw();
                c.SaveAs(Form("%s/nHits_sum_input_eq_1.png", cfg_.out_png_dir.c_str()));
                h_difBCID_passed->Draw();
                c.SaveAs(Form("%s/input_difBCID_passed.png", cfg_.out_png_dir.c_str()));
                h_difBCID_failed->Draw();
                c.SaveAs(Form("%s/input_difBCID_failed.png", cfg_.out_png_dir.c_str()));
                h2_BCID_2D->Draw("COLZ");
                c.SaveAs(Form("%s/input_BCID_2D.png", cfg_.out_png_dir.c_str()));
                h_clock_passed->Draw("histe1");
                h_clock_passed->SetMinimum(0);
                c.SaveAs(Form("%s/input_clock_passed.png", cfg_.out_png_dir.c_str()));
                h_clock_failed->Draw("histe1");
                h_clock_failed->SetMinimum(0);
                c.SaveAs(Form("%s/input_clock_failed.png", cfg_.out_png_dir.c_str()));
                h_nHits_passed->Draw("histe1");
                h_nHits_passed->SetMinimum(0);
                c.SaveAs(Form("%s/input_nHits_passed.png", cfg_.out_png_dir.c_str()));
                h_nHits_failed->Draw("histe1");
                h_nHits_failed->SetMinimum(0);
                c.SaveAs(Form("%s/input_nHits_failed.png", cfg_.out_png_dir.c_str()));
                h_deltaBCID->Draw("histe1");
                h_deltaBCID->SetMinimum(0);
                c.SaveAs(Form("%s/input_deltaBCID.png", cfg_.out_png_dir.c_str()));
                h_deltaCycleID->Draw("histe1");
                h_deltaCycleID->SetMinimum(0);
                c.SaveAs(Form("%s/input_deltaCycleID.png", cfg_.out_png_dir.c_str()));
                h_deltaTimestamp->Draw("histe1");
                h_deltaTimestamp->SetMinimum(0);
                c.SaveAs(Form("%s/input_deltaTimestamp.png", cfg_.out_png_dir.c_str()));
                h_deltaTriggerID->Draw("histe1");
                h_deltaTriggerID->SetMinimum(0);
                c.SaveAs(Form("%s/input_deltaTriggerID.png", cfg_.out_png_dir.c_str()));
                h_cycle_duration->Draw("histe1");
                h_cycle_duration->SetMinimum(0);
                c.SaveAs(Form("%s/input_cycle_duration.png", cfg_.out_png_dir.c_str()));
                h_loss->Draw("histe1");
                h_loss->SetMinimum(0);
                c.SaveAs(Form("%s/input_loss.png", cfg_.out_png_dir.c_str()));
                for (auto& eff : effs) {
                    eff->Draw("AP");
                    c.SaveAs(Form("%s/%s.png", cfg_.out_png_dir.c_str(), eff->GetName()));
                }
                for (int i = 0; i < cfg_.ntrigger_layer; ++i) {
                    h_double_lost[i]->Draw();
                    c.SaveAs(Form("%s/double_lost_input%d.png", cfg_.out_png_dir.c_str(), i));
                    h_single_lost[i]->Draw();
                    c.SaveAs(Form("%s/single_lost_input%d.png", cfg_.out_png_dir.c_str(), i));
                }
                for (int i = 0; i < cfg_.ntrigger_layer; ++i) {
                    for (int j = 0; j < cfg_.ntrigger_layer; ++j) {
                        if (i >= j) continue; // only plot upper triangle to avoid duplication
                        h_timing_diff_double_lost[i*cfg_.ntrigger_layer+j]->Draw("hist");
                        // h_timing_diff_double_lost[i]->SetMinimum(0);
                        c.SaveAs(Form("%s/timing_diff_double_lost_%d_%d.png", cfg_.out_png_dir.c_str(), i, j));
                        h_timing_diff_single_lost[i*cfg_.ntrigger_layer+j]->Draw("hist");
                        // h_timing_diff_single_lost[i]->SetMinimum(0);
                        c.SaveAs(Form("%s/timing_diff_single_lost_%d_%d.png", cfg_.out_png_dir.c_str(), i, j));
                        h_timing_diff_no_lost[i*cfg_.ntrigger_layer+j]->Draw("hist");
                        // h_timing_diff_no_lost[i]->SetMinimum(0);
                        c.SaveAs(Form("%s/timing_diff_no_lost_%d_%d.png", cfg_.out_png_dir.c_str(), i, j));
                        TLegend*legend = new TLegend(0.6, 0.7, 0.9, 0.9);
                        legend->AddEntry(h_timing_diff_double_lost[i*cfg_.ntrigger_layer+j].get(), "Double Lost", "l");
                        legend->AddEntry(h_timing_diff_single_lost[i*cfg_.ntrigger_layer+j].get(), Form("Single Lost (Lost layer %d)", j), "l");
                        legend->AddEntry(h_timing_diff_no_lost[i*cfg_.ntrigger_layer+j].get(), "No Lost", "l");
                        TCanvas c_comparison(Form("c_timing_diff_comparison_%d_%d", i, j), Form("Timing Difference Comparison for Layers %d and %d", i, j), 800, 600);
                        h_timing_diff_double_lost[i*cfg_.ntrigger_layer+j]->SetLineColor(kRed);
                        h_timing_diff_single_lost[i*cfg_.ntrigger_layer+j]->SetLineColor(kBlue);
                        h_timing_diff_no_lost[i*cfg_.ntrigger_layer+j]->SetLineColor(kGreen);
                        h_timing_diff_no_lost[i*cfg_.ntrigger_layer+j]->SetTitle(Form("Timing Difference Comparison for Layers %d - %d;Timing Difference (0.625 ns);Events", j, i));
                        h_timing_diff_no_lost[i*cfg_.ntrigger_layer+j]->Draw("hist");
                        h_timing_diff_double_lost[i*cfg_.ntrigger_layer+j]->Draw("histsame");
                        h_timing_diff_single_lost[i*cfg_.ntrigger_layer+j]->Draw("histsame");
                        legend->Draw();
                        c_comparison.SaveAs(Form("%s/timing_diff_comparison_%d_%d.png", cfg_.out_png_dir.c_str(), i, j));
                    }
                }
                delete h_dummy;
            }
        
            fout->Close();
            LOG_INFO("InputEffAlg: wrote {}", cfg_.out_filename);
        }

        InputEffAlgCfg cfg_;
        bool written_ = false;
        void initialize_histograms() {
            h_full = std::make_unique<TH1D>("h_full", "Denominator for HitTag efficiency;Input;Events", 6, -0.5, 5.5);
            h_passed = std::make_unique<TH1D>("h_passed", "Numerator for HitTag efficiency;Input;Events", 6, -0.5, 5.5);
            h_input_sum = std::make_unique<TH1D>("h_input_sum", "Distribution of sum of Inputs;Sum of Inputs;Events", 6, -0.5, 5.5);
            h_inputsumeqq1 = std::make_unique<TH1D>("h_inputsumeqq1", "Distribution of HitTag for events with sum of Inputs == 1;Input;Events", 6, -0.5, 5.5);
            h_nHits_sum_inputs_eq_1 = std::make_unique<TH1D>("h_nHits_sum_inputs_eq_1", "Distribution of nHits for events with sum of Inputs == 1;nHits;Events", 100, 0, 100);
            h_difBCID_passed = std::make_unique<TH1D>("h_difBCID_passed", "Distribution of BCID difference for events that passed the HitTag;BCID difference;Events", 1024, -1024, 1024);
            h_difBCID_failed = std::make_unique<TH1D>("h_difBCID_failed", "Distribution of BCID difference for events that failed the HitTag;BCID difference;Events", 1024, -1024, 1024);
            h2_BCID_2D = std::make_unique<TH2D>("h2_BCID_2D", "2D distribution of BCID; TLU BCID; Hits BCID", 1024, 0, 65516, 1024, 0, 1024);
            h_clock_passed = std::make_unique<TH1D>("h_clock_passed", "Distribution of clock for events that passed the HitTag;Clock;Events", 8, -0.5, 7.5);
            h_clock_failed = std::make_unique<TH1D>("h_clock_failed", "Distribution of clock for events that failed the HitTag;Clock;Events", 8, -0.5, 7.5);
            h_nHits_passed = std::make_unique<TH1D>("h_nHits_passed", "Distribution of nHits for events that passed the HitTag;nHits;Events", 100, 0, 100);
            h_nHits_failed = std::make_unique<TH1D>("h_nHits_failed", "Distribution of nHits for events that failed the HitTag;nHits;Events", 100, 0, 100);
            h_deltaBCID = std::make_unique<TH1D>("h_deltaBCID", "Distribution of BCID difference between AHCAL hit and TLU;#Delta BCID;Events", 1024, -1024, 1024);
            h_deltaCycleID = std::make_unique<TH1D>("h_deltaCycleID", "Distribution of CycleID difference between AHCAL hit and TLU;#Delta CycleID;Events", 1024, -1024, 1024);
            h_deltaTimestamp = std::make_unique<TH1D>("h_deltaTimestamp", "Distribution of Timestamp difference between AHCAL hit and TLU;#Delta Timestamp;Events", 1024, -1024, 1024);
            h_deltaTriggerID = std::make_unique<TH1D>("h_deltaTriggerID", "Distribution of TriggerID difference between AHCAL hit and TLU;#Delta TriggerID;Events", 1024, 0, 1024);
            h_cycle_duration = std::make_unique<TH1D>("h_cycle_duration", "Distribution of cycle duration calculated from TLU timestamp and BCID difference;Cycle duration (ns);Events", 100, 0, 10000000);
            h_loss = std::make_unique<TH1D>("h_loss", "Distribution of syncronization of the lost trigger;Loss category;Events", 5, -0.5, 4.5);
            h_full->SetDirectory(nullptr);
            h_passed->SetDirectory(nullptr);
            h_input_sum->SetDirectory(nullptr);
            h_inputsumeqq1->SetDirectory(nullptr);
            h_nHits_sum_inputs_eq_1->SetDirectory(nullptr);
            h_difBCID_passed->SetDirectory(nullptr);
            h_difBCID_failed->SetDirectory(nullptr);
            h2_BCID_2D->SetDirectory(nullptr);
            h_clock_passed->SetDirectory(nullptr);
            h_clock_failed->SetDirectory(nullptr);
            h_nHits_passed->SetDirectory(nullptr);
            h_nHits_failed->SetDirectory(nullptr);
            h_deltaBCID->SetDirectory(nullptr);
            h_deltaCycleID->SetDirectory(nullptr);
            h_deltaTimestamp->SetDirectory(nullptr);
            h_deltaTriggerID->SetDirectory(nullptr);
            h_cycle_duration->SetDirectory(nullptr);
            h_loss->SetDirectory(nullptr);

            //
            h_double_lost.clear();
            for (int i = 0; i < 4; ++i) {
                h_double_lost.push_back(std::make_unique<TH1D>(Form("h_double_lost_%d", i), Form("Distribution of double lost pattern for layer %d;Input;Events", i), 4, -0.5, 3.5));
                h_double_lost[i]->SetDirectory(nullptr);
            }
            h_single_lost.clear();
            for (int i = 0; i < 4; ++i) {
                h_single_lost.push_back(std::make_unique<TH1D>(Form("h_single_lost_%d", i), Form("Distribution of single lost for layer %d;Input;Events", i), 4, -0.5, 3.5));
                h_single_lost[i]->SetDirectory(nullptr);
            }
            h_timing_diff_double_lost.clear();
            h_timing_diff_single_lost.clear();
            h_timing_diff_no_lost.clear();
            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 4; ++j) {
                    h_timing_diff_double_lost.push_back(std::make_unique<TH1D>(Form("h_timing_diff_double_lost_%d_%d", i, j), Form("Timing difference between layers %d - %d for double lost pattern;Time difference [ns];Events", i, j), 64, -256., 256));
                    h_timing_diff_double_lost.back()->SetDirectory(nullptr);
                    h_timing_diff_single_lost.push_back(std::make_unique<TH1D>(Form("h_timing_diff_single_lost_%d_%d", i, j), Form("Timing difference between lost layer and other layer for single lost pattern: lost layer %d, other layer %d;Time difference [ns];Events", i, j), 64, -256., 256));
                    h_timing_diff_single_lost.back()->SetDirectory(nullptr);
                    h_timing_diff_no_lost.push_back(std::make_unique<TH1D>(Form("h_timing_diff_no_lost_%d_%d", i, j), Form("Timing difference between layers %d and %d for no lost pattern;Time difference [ns];Events", i, j), 64, -256., 256));
                    h_timing_diff_no_lost.back()->SetDirectory(nullptr);
                }
            }
        }

        std::unique_ptr<TH1D> h_full; // HitTag efficiency denominator 
        std::unique_ptr<TH1D> h_passed; // HitTag efficiency numerator 
        std::unique_ptr<TH1D> h_input_sum; // distribution of sum of inputs, for sanity check
        std::unique_ptr<TH1D> h_inputsumeqq1; // distribution of sum of inputs == 1, for sanity check
        std::unique_ptr<TH1D> h_nHits_sum_inputs_eq_1; // distribution of nHits for events with sum of inputs == 1, for sanity check
        std::unique_ptr<TH1D> h_difBCID_passed; // distribution of BCID difference between hit and trigger for events that passed the HitTag, for sanity check
        std::unique_ptr<TH1D> h_difBCID_failed; // distribution of BCID difference between hit and trigger for events that failed the HitTag, for sanity check
        std::unique_ptr<TH2D> h2_BCID_2D; // 2D distribution of BCID vs sum of inputs, for sanity check
        std::unique_ptr<TH1D> h_clock_passed; // distribution of clock for events that passed the HitTag, for sanity check
        std::unique_ptr<TH1D> h_clock_failed; // distribution of clock for events that failed the HitTag, for sanity check
        std::unique_ptr<TH1D> h_nHits_passed; // distribution of nHits for events that passed the HitTag, for sanity check
        std::unique_ptr<TH1D> h_nHits_failed; // distribution of nHits for events that failed the HitTag, for sanity check
        std::unique_ptr<TH1D> h_deltaBCID; 
        std::unique_ptr<TH1D> h_deltaCycleID;
        std::unique_ptr<TH1D> h_deltaTimestamp;
        std::unique_ptr<TH1D> h_deltaTriggerID;
        std::unique_ptr<TH1D> h_cycle_duration;
        std::unique_ptr<TH1D> h_loss;
        std::vector<std::unique_ptr<TH1D>> h_double_lost; // distribution of single/double lost pattern, for sanity check
        std::vector<std::unique_ptr<TH1D>> h_single_lost; // distribution of single/double lost pattern, for sanity check
        std::vector<std::unique_ptr<TH1D> > h_timing_diff_double_lost; // distribution of timing difference between the two layers for double lost pattern
        std::vector<std::unique_ptr<TH1D> > h_timing_diff_single_lost; // distribution of timing difference between the lost layer and the other layer for single lost pattern
        std::vector<std::unique_ptr<TH1D> > h_timing_diff_no_lost; // distribution of timing difference between the two layers for events with no lost input, for comparison with the above two categories
    };

    void InputEffAlg::ImplDeleter::operator()(Impl* p) const {
        delete p;
    }
    InputEffAlg::~InputEffAlg(){
        if (impl_) impl_->write();
    }
    bool InputEffAlg::is_hittag1_exist(std::vector<AHCALRawHit>& rawHits, int layer) {
        for (const auto& hit : rawHits) {
            if (hit.layer() == layer && hit.hittag == 1) {
                return true;
            }
        }
        return false;
    } 
    void InputEffAlg::execute(EventStore& evt){
        if (!impl_){
            impl_.reset(new Impl(cfg_));
            impl_->initialize_histograms();
        }
        auto rawdata = evt.get<AHCALTLURawData>(cfg_.in_data_key);
        auto rawHits = evt.get<std::vector<AHCALRawHit>>(cfg_.in_raw_hit_key);
        impl_->fill_inputs(rawdata.Inputs);
        int current_cycle_id = rawdata.CycleID;
        int current_AHCAL_bcid = rawHits.size() > 0 ? rawHits[0].bcid : -1; // assuming rawHits are sorted by time and the first one is the trigger hit, which should be checked in real data
        int current_TLU_timestamp = rawdata.Timestamp;
        int current_trigger_id = rawdata.TriggerID;
        if (last_cycle_id != -1 && last_AHCAL_bcid != -1 && last_TLU_timestamp != -1) {
            int delta_cycle_id = current_cycle_id - last_cycle_id;
            int delta_bcid = current_AHCAL_bcid - last_AHCAL_bcid;
            int delta_timestamp = current_TLU_timestamp - last_TLU_timestamp;
            int delta_trigger_id = current_trigger_id - last_trigger_id;
            impl_->h_deltaCycleID->Fill(delta_cycle_id);
            impl_->h_deltaBCID->Fill(delta_bcid);
            impl_->h_deltaTimestamp->Fill(delta_timestamp);
            impl_->h_deltaTriggerID->Fill(delta_trigger_id);
            if (delta_cycle_id > 1 && delta_timestamp > 0 && delta_trigger_id == 1) { // sanity check to avoid unphysical values, which should be the case for real data but might not be true for some MC samples
                double cycle_duration = (delta_timestamp*25 - delta_bcid*200) / delta_cycle_id; 
                impl_->h_cycle_duration->Fill(cycle_duration);
            }
        }
        last_cycle_id = current_cycle_id;
        last_AHCAL_bcid = current_AHCAL_bcid;
        last_TLU_timestamp = current_TLU_timestamp;
        last_trigger_id = current_trigger_id;
        // categorize the inputs 
        std::vector<bool> lost_inputs;
        for (int l = 0; l < cfg_.ntrigger_layer; ++l){
            if (rawdata.Inputs[l]) {
                bool is_has_hit = is_hittag1_exist(rawHits, cfg_.trigger_layer[l]);
                int nHits = 0;
                for (const auto& hit : rawHits) {
                    if (hit.hittag == 1) {
                        nHits++;
                    }
                }
                impl_->fill_single(l, nHits, is_has_hit);
                if (rawHits.size() > 0) // check if there is at least one hit to avoid out-of-range access, which should be the case for real data but might not be true for some MC samples
                    impl_->fill_difBCID(rawdata.BCID_TLU,rawHits[0].bcid, is_has_hit); // assuming rawHits are sorted by time and the first one is the trigger hit, which should be checked in real data
                if(!is_has_hit){
                    lost_inputs.push_back(true);
                    LOG_DEBUG("InputEffAlg: event {} has input {} but no hittag1 in layer {}, total nHits = {}", evt.event_counter(), l, cfg_.trigger_layer[l], nHits);
                } else {
                    lost_inputs.push_back(false);
                }
            }            
        }
        impl_->fill_loss_category(lost_inputs);
        // single lost / double lost with several pattern 
        int layer1 = -1;
        int layer2 = -1;
        bool isdoublelost = false;
        bool issinglelost = false;
        int lost_layer = -1;

        if (std::accumulate(rawdata.Inputs.begin(), rawdata.Inputs.end(), 0) == 2) {
            for (int l1 = 0; l1 < cfg_.ntrigger_layer; ++l1){
                for (int l2 = l1+1; l2 < cfg_.ntrigger_layer; ++l2){
                    if (rawdata.Inputs[l1] && rawdata.Inputs[l2]) {
                        layer1 = l1;
                        layer2 = l2;
                        if (!is_hittag1_exist(rawHits, cfg_.trigger_layer[l1]) && !is_hittag1_exist(rawHits, cfg_.trigger_layer[l2])) {
                            isdoublelost = true;
                            LOG_DEBUG("InputEffAlg: event {} has double lost inputs {} and {}, total nHits = {}", evt.event_counter(), l1, l2, rawHits.size());
                        } else if (!is_hittag1_exist(rawHits, cfg_.trigger_layer[l1]) && is_hittag1_exist(rawHits, cfg_.trigger_layer[l2])) {
                            issinglelost = true;
                            lost_layer = l1;
                        } else if (is_hittag1_exist(rawHits, cfg_.trigger_layer[l1]) && !is_hittag1_exist(rawHits, cfg_.trigger_layer[l2])) {
                            issinglelost = true;
                            lost_layer = l2;
                        }
                    }
                }
            }
            impl_->fill_pattern(layer1, layer2, isdoublelost, issinglelost, lost_layer, rawdata.FineTimestamps);
        }
        if (std::accumulate(rawdata.Inputs.begin(), rawdata.Inputs.end(), 0) == 0) {
            LOG_DEBUG("InputEffAlg: event {} has no input, total nHits = {}", evt.event_counter(), rawHits.size());
        }
        // sum_inputs == 1 analysis
        if (std::accumulate(rawdata.Inputs.begin(), rawdata.Inputs.end(), 0) == 1) {
            for (int l = 0; l < cfg_.ntrigger_layer; ++l){
                if (rawdata.Inputs[l]) {
                    bool is_has_hit = is_hittag1_exist(rawHits, cfg_.trigger_layer[l]);
                    int nHits = 0;
                    for (const auto& hit : rawHits) {
                        if (hit.layer() == cfg_.trigger_layer[l] && hit.hittag == 1) {
                            nHits++;
                        }
                    }
                    impl_->fill_single_input1(l, nHits, evt.event_counter(), this->ctx(), is_has_hit, out_file_.get());
                }
            }
        }
    }
    void InputEffAlg::parse_cfg(const YAML::Node& cfg){
        cfg_.in_data_key = get_or<std::string>(cfg, "in_data_key", cfg_.in_data_key);
        cfg_.in_raw_hit_key = get_or<std::string>(cfg, "in_raw_hit_key", cfg_.in_raw_hit_key);
        // cfg_.string_track_struct = get_or<std::string>(cfg, "string_track_struct", cfg_.string_track_struct);
        // cfg_.track_selection_string = get_or<std::string>(cfg, "track_selection_string", cfg_.track_selection_string);
        cfg_.out_filename = get_or<std::string>(cfg, "out_filename", cfg_.out_filename);
        cfg_.write_to_png = get_or<bool>(cfg, "write_to_png", cfg_.write_to_png);
        cfg_.out_png_dir = get_or<std::string>(cfg, "out_png_dir", cfg_.out_png_dir);
        cfg_.write_to_txt = get_or<bool>(cfg, "write_to_txt", cfg_.write_to_txt);
        cfg_.out_txt_name = get_or<std::string>(cfg, "out_txt_name", cfg_.out_txt_name);
        cfg_.ntrigger_layer = get_or<int>(cfg, "ntrigger_layer", cfg_.ntrigger_layer);
        cfg_.trigger_layer.clear(); //{1:9, 2:19, 3:29, 4:38}
        std::map<int, int> default_trigger_layer = {{0,9}, {1,19}, {2,29}, {3,38}};
        default_trigger_layer = get_or<std::map<int,int>>(cfg, "trigger_layer", default_trigger_layer);
        for (int i = 0; i < cfg_.ntrigger_layer; ++i) {
            if (default_trigger_layer.find(i) != default_trigger_layer.end()) {
                cfg_.trigger_layer.push_back(default_trigger_layer[i]);
            } else {
                LOG_ERROR("InputEffAlg: trigger_layer config is missing entry for input index {}, skipping this input", i);
            }
        }
    }
    void InputEffAlg::initialize() {
        if (cfg_.write_to_txt) {
            out_file_ = std::make_unique<std::ofstream>();
            out_file_->open(cfg_.out_txt_name, std::ios::out);
            if (!out_file_->is_open()) {
                LOG_ERROR("InputEffAlg: cannot open output text file: {}", cfg_.out_txt_name);
                out_file_.reset();
            } else {
                LOG_INFO("InputEffAlg: writing detailed event info to {}", cfg_.out_txt_name);
            }
        }
    }
} // namespace AHCALRecoAlg