#include "fft_helper.h"

#include <TAxis.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TF1.h>
#include <TFitResultPtr.h>
#include <TH1.h>
#include <TH1D.h>
#include <TKey.h>
#include <TObject.h>
#include <TObjArray.h>
#include <TSpectrum.h>
#include <TString.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace {

struct SpeHistRef {
  std::string path;
  int cellid = -1;
  int layer = -1;
  int chip = -1;
  int channel = -1;
};

struct GausPeakFit {
  int nPeaks = 0;
  int nFitPeaks = 0;
  int nDiffs = 0;
  int fitStatus = -1;
  double peakX = -1.0;
  double peakAmp = -1.0;
  double fitMin = -1.0;
  double fitMax = -1.0;
  double spe = -1.0;
  double speErr = -1.0;
  double sigma = -1.0;
  double sigmaErr = -1.0;
  double chi2ndf = -1.0;
  std::vector<double> peakXs;
  std::vector<double> peakAmps;
  std::vector<double> fitMeans;
  std::vector<double> fitMeanErrs;
  std::vector<double> fitSigmas;
  std::vector<double> fitSigmaErrs;
  std::vector<double> fitChi2Ndfs;
  std::vector<double> peakDiffs;
  std::vector<double> peakRawDiffs;
  std::vector<int> peakDiffOrders;
  std::vector<int> fitStatuses;
};

bool parseSpeHistName(const std::string &name, SpeHistRef &ref) {
  int layer = -1;
  int chip = -1;
  int channel = -1;
  if (std::sscanf(name.c_str(), "hSPE_L%d_C%d_ch%d", &layer, &chip, &channel) == 3) {
    ref.layer = layer;
    ref.chip = chip;
    ref.channel = channel;
    ref.cellid = layer * 100000 + chip * 10000 + channel;
    return true;
  }

  int cellid = -1;
  int consumed = 0;
  if (std::sscanf(name.c_str(), "hSPE_%d%n", &cellid, &consumed) == 1 &&
      consumed == static_cast<int>(name.size())) {
    ref.cellid = cellid;
    ref.layer = cellid / 100000;
    ref.chip = (cellid / 10000) % 10;
    ref.channel = cellid % 10000;
    return true;
  }

  return false;
}

std::string joinRootPath(const std::string &dir, const std::string &name) {
  if (dir.empty()) return name;
  return dir + "/" + name;
}

void collectSpeHists(TDirectory *dir, const std::string &dirPath, std::vector<SpeHistRef> &refs) {
  if (!dir) return;

  TIter next(dir->GetListOfKeys());
  while (TKey *key = static_cast<TKey *>(next())) {
    std::unique_ptr<TObject> obj(key->ReadObj());
    if (!obj) continue;

    const std::string name = key->GetName();
    const std::string path = joinRootPath(dirPath, name);

    if (obj->InheritsFrom(TDirectory::Class())) {
      collectSpeHists(static_cast<TDirectory *>(obj.get()), path, refs);
      continue;
    }

    if (!obj->InheritsFrom(TH1::Class())) continue;

    SpeHistRef ref;
    if (!parseSpeHistName(name, ref)) continue;
    ref.path = path;
    refs.push_back(ref);
  }
}

std::vector<double> histogramBins(const TH1D *hist) {
  std::vector<double> values(hist ? hist->GetNbinsX() : 0, 0.0);
  for (int bin = 1; hist && bin <= hist->GetNbinsX(); ++bin) {
    values[bin - 1] = hist->GetBinContent(bin);
  }

  double mean = 0.0;
  for (const double value : values) mean += value;
  if (!values.empty()) mean /= static_cast<double>(values.size());
  for (double &value : values) value -= mean;
  return values;
}

GausPeakFit fitPhotonPeaksWithGaussian(TH1D *hist,
                                       double searchMin,
                                       double searchMax,
                                       double minRelHeight,
                                       double fitHalfWidth,
                                       double rootPeakSigma,
                                       double diffMin,
                                       double diffMax,
                                       double expectedSpe,
                                       const char *fitNamePrefix) {
  GausPeakFit result;
  if (!hist || hist->GetEntries() <= 0) return result;

  TAxis *xaxis = hist->GetXaxis();
  const int nbins = hist->GetNbinsX();
  searchMin = std::max(searchMin, xaxis->GetXmin());
  searchMax = std::min(searchMax, xaxis->GetXmax());
  if (searchMax <= searchMin) return result;

  const int binMin = std::max(1, xaxis->FindBin(searchMin));
  const int binMax = std::min(nbins, xaxis->FindBin(searchMax));
  if (binMax < binMin) return result;

  minRelHeight = std::max(0.0, minRelHeight);
  rootPeakSigma = std::max(1.0, rootPeakSigma);
  if (fitHalfWidth <= 0.0) {
    const double expectedHalfWidth = (expectedSpe > 0.0 ? 0.25 * expectedSpe : 5.0);
    fitHalfWidth = std::min(12.0, std::max(5.0, expectedHalfWidth));
  }
  diffMin = std::max(1.0, diffMin);
  if (diffMax <= diffMin) diffMax = std::max(diffMin, searchMax - searchMin);

  std::unique_ptr<TH1D> hSearch(static_cast<TH1D *>(
      hist->Clone(Form("%s_search", hist->GetName()))));
  hSearch->SetDirectory(nullptr);
  for (int bin = 1; bin <= nbins; ++bin) {
    if (bin < binMin || bin > binMax) hSearch->SetBinContent(bin, 0.0);
  }

  const int maxPeaks = std::max(8, std::min(100, binMax - binMin + 1));
  TSpectrum spectrum(maxPeaks);
  const int nFound = spectrum.Search(hSearch.get(), rootPeakSigma, "nodraw", minRelHeight);

  std::vector<std::pair<double, double>> candidates;
  double *peakXs = spectrum.GetPositionX();
  double *peakYs = spectrum.GetPositionY();
  for (int i = 0; i < nFound; ++i) {
    const double x = peakXs ? peakXs[i] : -1.0;
    if (x < searchMin || x > searchMax) continue;
    const int bin = xaxis->FindBin(x);
    if (bin < binMin || bin > binMax) continue;
    const double amp = peakYs ? peakYs[i] : hist->GetBinContent(bin);
    candidates.emplace_back(x, amp);
  }

  std::sort(candidates.begin(), candidates.end(),
            [](const auto &a, const auto &b) { return a.first < b.first; });

  const double binWidth = hist->GetBinWidth(std::max(1, binMin));
  std::vector<std::pair<double, double>> uniqueCandidates;
  for (const auto &candidate : candidates) {
    if (!uniqueCandidates.empty() &&
        std::abs(candidate.first - uniqueCandidates.back().first) < binWidth) {
      if (candidate.second > uniqueCandidates.back().second) {
        uniqueCandidates.back() = candidate;
      }
      continue;
    }
    uniqueCandidates.push_back(candidate);
  }
  candidates.swap(uniqueCandidates);
  result.nPeaks = static_cast<int>(candidates.size());

  std::vector<double> goodMeans;
  for (int i = 0; i < static_cast<int>(candidates.size()); ++i) {
    const double peakX = candidates[i].first;
    const int peakBin = xaxis->FindBin(peakX);
    const double peakAmp = hist->GetBinContent(peakBin);

    const double fitMin = std::max(xaxis->GetXmin(), peakX - fitHalfWidth);
    const double fitMax = std::min(xaxis->GetXmax(), peakX + fitHalfWidth);
    if (fitMax <= fitMin) continue;

    const TString fitName = Form("%s_peak%d", fitNamePrefix, i);
    TF1 gaus(fitName, "gaus", fitMin, fitMax);
    const double sigma0 = std::max(binWidth, fitHalfWidth / 3.0);
    const double sigmaMax = std::max(fitHalfWidth + binWidth, 1.5 * fitHalfWidth);
    gaus.SetParameters(std::max(peakAmp, 1.0), peakX, sigma0);
    gaus.SetParLimits(1, fitMin, fitMax);
    gaus.SetParLimits(2, 0.25 * binWidth, std::max(sigmaMax, binWidth));

    TFitResultPtr fitResult = hist->Fit(&gaus, "RQ0SN");
    const int fitStatus = static_cast<int>(fitResult);
    const double mean = gaus.GetParameter(1);
    const double meanErr = gaus.GetParError(1);
    const double sigma = gaus.GetParameter(2);
    const double sigmaErr = gaus.GetParError(2);
    const double chi2ndf = (gaus.GetNDF() > 0 ? gaus.GetChisquare() / gaus.GetNDF() : -1.0);

    result.peakXs.push_back(peakX);
    result.peakAmps.push_back(peakAmp);
    result.fitStatuses.push_back(fitStatus);
    result.fitMeans.push_back(mean);
    result.fitMeanErrs.push_back(meanErr);
    result.fitSigmas.push_back(sigma);
    result.fitSigmaErrs.push_back(sigmaErr);
    result.fitChi2Ndfs.push_back(chi2ndf);
    hist->GetListOfFunctions()->Add(gaus.Clone(fitName));

    if (!std::isfinite(mean) || !std::isfinite(sigma)) continue;
    if (mean < searchMin || mean > searchMax || sigma <= 0.0) continue;
    if (sigma > sigmaMax + 0.5 * binWidth) continue;
    goodMeans.push_back(mean);
  }

  std::sort(goodMeans.begin(), goodMeans.end());
  result.nFitPeaks = static_cast<int>(goodMeans.size());

  std::vector<double> spacingMeans = goodMeans;
  if (spacingMeans.size() < 2 && candidates.size() >= 2) {
    spacingMeans.clear();
    for (const auto &candidate : candidates) spacingMeans.push_back(candidate.first);
    std::sort(spacingMeans.begin(), spacingMeans.end());
  }

  struct PairDiff {
    double diff = 0.0;
    int i = -1;
    int j = -1;
  };
  std::vector<PairDiff> pairDiffs;
  for (int i = 0; i < static_cast<int>(spacingMeans.size()); ++i) {
    for (int j = i + 1; j < static_cast<int>(spacingMeans.size()); ++j) {
      const double diff = spacingMeans[j] - spacingMeans[i];
      if (diff > 0.0) pairDiffs.push_back({diff, i, j});
    }
  }

  constexpr int maxPeakOrder = 16;
  std::vector<double> speCandidates;
  for (const auto &pd : pairDiffs) {
    const int maxOrder = std::max(1, std::min(maxPeakOrder, static_cast<int>(std::floor(pd.diff / diffMin))));
    for (int order = 1; order <= maxOrder; ++order) {
      const double unit = pd.diff / static_cast<double>(order);
      if (diffMin > 0.0 && unit < diffMin) continue;
      if (diffMax > 0.0 && unit > diffMax) continue;
      speCandidates.push_back(unit);
    }
  }

  double bestScore = -1.0;
  double bestSpread = 1e99;
  std::vector<double> bestUnits;
  std::vector<double> bestRawDiffs;
  std::vector<int> bestOrders;

  for (const double candidate : speCandidates) {
    const double tolerance = std::max(1.0, 0.20 * candidate);
    double score = 0.0;
    std::vector<double> units;
    std::vector<double> rawDiffs;
    std::vector<int> orders;

    for (const auto &pd : pairDiffs) {
      const int order = std::max(1, static_cast<int>(std::llround(pd.diff / candidate)));
      if (order > maxPeakOrder) continue;
      const double unit = pd.diff / static_cast<double>(order);
      if (diffMin > 0.0 && unit < diffMin) continue;
      if (diffMax > 0.0 && unit > diffMax) continue;
      if (std::abs(unit - candidate) > tolerance) continue;
      units.push_back(unit);
      rawDiffs.push_back(pd.diff);
      orders.push_back(order);
      score += 1.0;
    }

    if (expectedSpe > 0.0) {
      const double priorTolerance = std::max(1.0, 0.30 * expectedSpe);
      if (std::abs(candidate - expectedSpe) <= priorTolerance) score += 3.0;
    }

    if (units.empty()) continue;
    double mean = 0.0;
    for (const double unit : units) mean += unit;
    mean /= static_cast<double>(units.size());
    double spread = 0.0;
    for (const double unit : units) {
      const double delta = unit - mean;
      spread += delta * delta;
    }
    spread = std::sqrt(spread / static_cast<double>(units.size()));

    if (score > bestScore || (score == bestScore && spread < bestSpread)) {
      bestScore = score;
      bestSpread = spread;
      bestUnits = std::move(units);
      bestRawDiffs = std::move(rawDiffs);
      bestOrders = std::move(orders);
    }
  }

  result.peakDiffs = std::move(bestUnits);
  result.peakRawDiffs = std::move(bestRawDiffs);
  result.peakDiffOrders = std::move(bestOrders);
  result.nDiffs = static_cast<int>(result.peakDiffs.size());

  if (!result.peakXs.empty()) {
    result.peakX = result.peakXs.front();
    result.peakAmp = result.peakAmps.front();
  }
  if (!result.fitMeans.empty()) {
    result.fitStatus = result.fitStatuses.front();
    result.fitMin = result.fitMeans.front() - fitHalfWidth;
    result.fitMax = result.fitMeans.front() + fitHalfWidth;
    result.sigma = result.fitSigmas.front();
    result.sigmaErr = result.fitSigmaErrs.front();
    result.chi2ndf = result.fitChi2Ndfs.front();
  }

  if (result.peakDiffs.empty()) return result;

  double sum = 0.0;
  for (const double diff : result.peakDiffs) sum += diff;
  result.spe = sum / static_cast<double>(result.peakDiffs.size());

  if (result.peakDiffs.size() == 1) {
    result.speErr = 0.0;
  } else {
    double var = 0.0;
    for (const double diff : result.peakDiffs) {
      const double delta = diff - result.spe;
      var += delta * delta;
    }
    var /= static_cast<double>(result.peakDiffs.size() - 1);
    result.speErr = std::sqrt(var / static_cast<double>(result.peakDiffs.size()));
  }
  return result;
}

TDirectory *mkdirPath(TDirectory *top, const char *path) {
  if (!top) return nullptr;
  TDirectory *dir = top;
  TString rest(path);
  TObjArray *tokens = rest.Tokenize("/");
  for (int i = 0; tokens && i < tokens->GetEntriesFast(); ++i) {
    const char *token = tokens->At(i)->GetName();
    if (!token || !*token) continue;
    TDirectory *next = dynamic_cast<TDirectory *>(dir->Get(token));
    if (!next) next = dir->mkdir(token);
    dir = next;
  }
  delete tokens;
  return dir;
}

} // namespace

void refit_spe_from_hists(const char *inputFile = "spe_analysis.root",
                          const char *outputFile = "spe_refit.root",
                          int npadReq = 8192,
                          double gainMin = 5.0,
                          double gainMax = 30.0,
                          int minEntFft = 500,
                          double edgeGainTol = 0.01,
                          int kedgeSkip = 4,
                          bool writeFftHists = true,
                          bool writeAdcHists = true,
                          bool doGausPeakFit = true,
                          double peakSearchMin = 50.0,
                          double peakSearchMax = 600.0,
                          double peakMinRelHeight = 0.10,
                          double gausFitHalfWidth = -1.0,
                          double rootPeakSigma = 2.0,
                          double gausDiffMin = -1.0,
                          double gausDiffMax = -1.0) {
  std::unique_ptr<TFile> fin(TFile::Open(inputFile, "READ"));
  if (!fin || fin->IsZombie()) {
    std::cerr << "[refit_spe_from_hists] cannot open input file: " << inputFile << std::endl;
    return;
  }

  TDirectory *speDir = dynamic_cast<TDirectory *>(fin->Get("SPE"));
  if (!speDir) {
    std::cerr << "[refit_spe_from_hists] directory 'SPE' not found in " << inputFile << std::endl;
    return;
  }

  std::vector<SpeHistRef> refs;
  collectSpeHists(speDir, "SPE", refs);
  std::sort(refs.begin(), refs.end(), [](const SpeHistRef &a, const SpeHistRef &b) {
    return a.cellid < b.cellid;
  });

  std::unique_ptr<TFile> fout(TFile::Open(outputFile, "RECREATE"));
  if (!fout || fout->IsZombie()) {
    std::cerr << "[refit_spe_from_hists] cannot create output file: " << outputFile << std::endl;
    return;
  }

  TTree tree("spe", "SPE refit from saved histograms");
  int cellid = -1;
  int layer = -1;
  int chip = -1;
  int channel = -1;
  int entries = 0;
  double x_mm = -999.0;
  double y_mm = -999.0;
  int npad = 0;
  double gain = -1.0;
  double gainErr = -1.0;
  double snr = -1.0;
  int kp = -1;
  double k0 = -1.0;
  double peakAmp = -1.0;
  double noiseMed = -1.0;
  double spe_gaus = -1.0;
  double spe_gaus_err = -1.0;
  int gaus_n_peaks = 0;
  int gaus_n_fit_peaks = 0;
  int gaus_n_diffs = 0;
  int gaus_fit_status = -1;
  double gaus_peak_x = -1.0;
  double gaus_peak_amp = -1.0;
  double gaus_sigma = -1.0;
  double gaus_sigma_err = -1.0;
  double gaus_chi2ndf = -1.0;
  double gaus_fit_min = -1.0;
  double gaus_fit_max = -1.0;
  std::vector<double> photon_peak_xs;
  std::vector<double> photon_peak_amps;
  std::vector<double> photon_peak_means;
  std::vector<double> photon_peak_mean_errs;
  std::vector<double> photon_peak_sigmas;
  std::vector<double> photon_peak_sigma_errs;
  std::vector<double> photon_peak_chi2ndfs;
  std::vector<double> photon_peak_diffs;
  std::vector<double> photon_peak_raw_diffs;
  std::vector<int> photon_peak_diff_orders;
  std::vector<int> photon_peak_fit_statuses;

  tree.Branch("cellid", &cellid);
  tree.Branch("layer", &layer);
  tree.Branch("chip", &chip);
  tree.Branch("channel", &channel);
  tree.Branch("entries", &entries);
  tree.Branch("x_mm", &x_mm);
  tree.Branch("y_mm", &y_mm);
  tree.Branch("npad", &npad);
  tree.Branch("gain", &gain);
  tree.Branch("gainErr", &gainErr);
  tree.Branch("snr", &snr);
  tree.Branch("kp", &kp);
  tree.Branch("k0", &k0);
  tree.Branch("peakAmp", &peakAmp);
  tree.Branch("noiseMed", &noiseMed);
  tree.Branch("spe_gaus", &spe_gaus);
  tree.Branch("spe_gaus_err", &spe_gaus_err);
  tree.Branch("gaus_n_peaks", &gaus_n_peaks);
  tree.Branch("gaus_n_fit_peaks", &gaus_n_fit_peaks);
  tree.Branch("gaus_n_diffs", &gaus_n_diffs);
  tree.Branch("gaus_fit_status", &gaus_fit_status);
  tree.Branch("gaus_peak_x", &gaus_peak_x);
  tree.Branch("gaus_peak_amp", &gaus_peak_amp);
  tree.Branch("gaus_sigma", &gaus_sigma);
  tree.Branch("gaus_sigma_err", &gaus_sigma_err);
  tree.Branch("gaus_chi2ndf", &gaus_chi2ndf);
  tree.Branch("gaus_fit_min", &gaus_fit_min);
  tree.Branch("gaus_fit_max", &gaus_fit_max);
  tree.Branch("photon_peak_xs", &photon_peak_xs);
  tree.Branch("photon_peak_amps", &photon_peak_amps);
  tree.Branch("photon_peak_means", &photon_peak_means);
  tree.Branch("photon_peak_mean_errs", &photon_peak_mean_errs);
  tree.Branch("photon_peak_sigmas", &photon_peak_sigmas);
  tree.Branch("photon_peak_sigma_errs", &photon_peak_sigma_errs);
  tree.Branch("photon_peak_chi2ndfs", &photon_peak_chi2ndfs);
  tree.Branch("photon_peak_diffs", &photon_peak_diffs);
  tree.Branch("photon_peak_raw_diffs", &photon_peak_raw_diffs);
  tree.Branch("photon_peak_diff_orders", &photon_peak_diff_orders);
  tree.Branch("photon_peak_fit_statuses", &photon_peak_fit_statuses);

  TDirectory *fftTop = writeFftHists ? mkdirPath(fout.get(), "FFT") : nullptr;
  TDirectory *adcTop = writeAdcHists ? mkdirPath(fout.get(), "ADC") : nullptr;
  int nFit = 0;
  int nGausFit = 0;
  int i = 0;
  int success_fft_inlayer[40] = {0};
  int success_gaus_inlayer[40] = {0};
  int total_inlayer[40] = {0};
  for (const SpeHistRef &ref : refs) {
    TH1D *hist = dynamic_cast<TH1D *>(fin->Get(ref.path.c_str()));
    if (i % 100 == 0) {
      std::cout << "[refit_spe_from_hists] processing " << i << "/" << refs.size()
                << " histograms, fitted " << nFit << std::endl;
    }
    i++;
    if (!hist) continue;
    cellid = ref.cellid;
    layer = ref.layer;
    ++total_inlayer[layer];
    chip = ref.chip;
    channel = ref.channel;
    entries = static_cast<int>(hist->GetEntries());
    x_mm = -999.0;
    y_mm = -999.0;
    npad = 0;
    gain = -1.0;
    gainErr = -1.0;
    snr = -1.0;
    kp = -1;
    k0 = -1.0;
    peakAmp = -1.0;
    noiseMed = -1.0;
    spe_gaus = -1.0;
    spe_gaus_err = -1.0;
    gaus_n_peaks = 0;
    gaus_n_fit_peaks = 0;
    gaus_n_diffs = 0;
    gaus_fit_status = -1;
    gaus_peak_x = -1.0;
    gaus_peak_amp = -1.0;
    gaus_sigma = -1.0;
    gaus_sigma_err = -1.0;
    gaus_chi2ndf = -1.0;
    gaus_fit_min = -1.0;
    gaus_fit_max = -1.0;
    photon_peak_xs.clear();
    photon_peak_amps.clear();
    photon_peak_means.clear();
    photon_peak_mean_errs.clear();
    photon_peak_sigmas.clear();
    photon_peak_sigma_errs.clear();
    photon_peak_chi2ndfs.clear();
    photon_peak_diffs.clear();
    photon_peak_raw_diffs.clear();
    photon_peak_diff_orders.clear();
    photon_peak_fit_statuses.clear();

    std::unique_ptr<TH1D> hAdc(static_cast<TH1D *>(
        hist->Clone(Form("hADC_L%d_C%d_ch%d", layer, chip, channel))));
    hAdc->SetDirectory(nullptr);

    if (entries >= minEntFft) {
      const std::vector<double> bins = histogramBins(hist);
      npad = AHCALRecoAlg::chooseNpad(static_cast<int>(bins.size()), npadReq);
      std::unique_ptr<TH1D> hFFT(AHCALRecoAlg::doFFT(
          bins, npad, Form("hFFT_k_L%d_C%d_ch%d", layer, chip, channel)));

      if (hFFT) {
        const AHCALRecoAlg::KPick fit = AHCALRecoAlg::pickK0(
            hFFT.get(), npad, hist->GetBinWidth(1), gainMin, gainMax, edgeGainTol, kedgeSkip);
        gain = fit.gain;
        gainErr = fit.gainErr;
        snr = fit.snr;
        kp = fit.kp;
        k0 = fit.k0;
        peakAmp = fit.peakAmp;
        noiseMed = fit.noiseMed;
        ++nFit;
        if (snr > 3.0 && gain > 0.0) {
          ++success_fft_inlayer[layer];
        }

        if (fftTop) {
          TDirectory *dir = mkdirPath(
              fftTop, Form("Layer%02d/Chip%d", std::max(layer, 0), std::max(chip, 0)));
          if (dir) dir->cd();
          hFFT->Write();
          std::unique_ptr<TH1D> hGain(AHCALRecoAlg::makeFFTgainHist(
              hFFT.get(), npad, hist->GetBinWidth(1), gainMin, gainMax,
              Form("hFFT_gain_L%d_C%d_ch%d", layer, chip, channel)));
          hGain->Write();
        }
      }
    }

    if (doGausPeakFit && entries >= minEntFft) {
      const double searchMin = (peakSearchMin >= 0.0 ? peakSearchMin : hAdc->GetXaxis()->GetXmin());
      const double searchMax = (peakSearchMax >= 0.0 ? peakSearchMax : hAdc->GetXaxis()->GetXmax());
      const double expectedSpe = (gain > 0.0 ? gain : -1.0);
      const double gausUnitMin = (gausDiffMin > 0.0 ? gausDiffMin : gainMin);
      const double gausUnitMax = (gausDiffMax > 0.0 ? gausDiffMax : 1.5 * gainMax);
      const GausPeakFit gausFit = fitPhotonPeaksWithGaussian(
          hAdc.get(), searchMin, searchMax, peakMinRelHeight, gausFitHalfWidth, rootPeakSigma,
          gausUnitMin, gausUnitMax, expectedSpe,
          Form("fGausPeak_L%d_C%d_ch%d", layer, chip, channel));

      spe_gaus = gausFit.spe;
      spe_gaus_err = gausFit.speErr;
      gaus_n_peaks = gausFit.nPeaks;
      gaus_n_fit_peaks = gausFit.nFitPeaks;
      gaus_n_diffs = gausFit.nDiffs;
      gaus_fit_status = gausFit.fitStatus;
      gaus_peak_x = gausFit.peakX;
      gaus_peak_amp = gausFit.peakAmp;
      gaus_sigma = gausFit.sigma;
      gaus_sigma_err = gausFit.sigmaErr;
      gaus_chi2ndf = gausFit.chi2ndf;
      gaus_fit_min = gausFit.fitMin;
      gaus_fit_max = gausFit.fitMax;
      photon_peak_xs = gausFit.peakXs;
      photon_peak_amps = gausFit.peakAmps;
      photon_peak_means = gausFit.fitMeans;
      photon_peak_mean_errs = gausFit.fitMeanErrs;
      photon_peak_sigmas = gausFit.fitSigmas;
      photon_peak_sigma_errs = gausFit.fitSigmaErrs;
      photon_peak_chi2ndfs = gausFit.fitChi2Ndfs;
      photon_peak_diffs = gausFit.peakDiffs;
      photon_peak_raw_diffs = gausFit.peakRawDiffs;
      photon_peak_diff_orders = gausFit.peakDiffOrders;
      photon_peak_fit_statuses = gausFit.fitStatuses;
      if (gaus_n_diffs > 0 && spe_gaus > 0.0) {
        ++nGausFit;
        ++success_gaus_inlayer[layer];
      }
    }

    if (adcTop) {
      TDirectory *dir = mkdirPath(
          adcTop, Form("Layer%02d/Chip%d", std::max(layer, 0), std::max(chip, 0)));
      if (dir) dir->cd();
      hAdc->Write();
    }

    tree.Fill();
  }
  TEfficiency* eff_fft = new TEfficiency("efficiency_fft", "FFT fit efficiency", 40, 0.0, 40);
  TEfficiency* eff_gaus = new TEfficiency("efficiency_gaus", "Gaussian peak fit efficiency", 40, 0.0, 40);
  for (int l = 0; l < 40; ++l) {
    eff_fft->SetTotalEvents(l, total_inlayer[l]);
    eff_fft->SetPassedEvents(l, success_fft_inlayer[l]);
    eff_gaus->SetTotalEvents(l, total_inlayer[l]);
    eff_gaus->SetPassedEvents(l, success_gaus_inlayer[l]);
  }
  // if (fftTop) fftTop->Write();
  // if (adcTop) adcTop->Write();
  fout->cd();
  tree.Write();
  eff_fft->Write();
  eff_gaus->Write();
  fout->Close();

  std::cout << "[refit_spe_from_hists] scanned " << refs.size()
            << " saved SPE histograms, FFT-fitted " << nFit
            << ", Gaussian-fitted " << nGausFit
            << ", wrote " << outputFile << std::endl;
}
