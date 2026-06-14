#pragma once
#include <vector>
#include <algorithm>
#include <cmath>
#include "TH1D.h"
#include "TVirtualFFT.h"

namespace AHCALRecoAlg {

// ============ Utility Functions ============

/// Compute median of a vector (destructive - will reorder)
template<class V>
inline double median(V v) {
    if (v.empty()) return 0.0;
    size_t n = v.size() / 2;
    std::nth_element(v.begin(), v.begin() + n, v.end());
    double m = v[n];
    if (v.size() % 2 == 0) {
        std::nth_element(v.begin(), v.begin() + n - 1, v.end());
        m = 0.5 * (m + v[n - 1]);
    }
    return m;
}

/// Compute Npad for FFT (power of 2, at least Npad_req)
inline int chooseNpad(int nb, int Npad_req = 8192) {
    if (nb <= Npad_req) return Npad_req;
    int p = 1;
    while (p < nb) p <<= 1;
    return p;
}

// ============ FFT Processing ============

/// Perform real-to-complex FFT on zero-padded input
/// Returns TH1D* with magnitude spectrum, or nullptr on error
TH1D* doFFT(const std::vector<double>& vin, int Npad, const char* name) {
    if (vin.empty()) return nullptr;

    std::vector<double> vpad(Npad, 0.0);
    for (size_t i = 0; i < vin.size() && i < (size_t)Npad; ++i) {
        vpad[i] = vin[i];
    }

    Int_t n = Npad;
    TVirtualFFT* fft = TVirtualFFT::FFT(1, &n, "R2C ES K");
    if (!fft) return nullptr;

    fft->SetPoints(vpad.data());
    fft->Transform();

    std::vector<double> re(Npad / 2 + 1), im(Npad / 2 + 1);
    fft->GetPointsComplex(re.data(), im.data());
    delete fft;

    // ROOT bin q + 1 stores FFT frequency index q.
    TH1D* h = new TH1D(name, name, Npad / 2, -0.5, Npad / 2 - 0.5);
    h->SetDirectory(nullptr);
    for (int k = 0; k < Npad / 2; ++k) {
        h->SetBinContent(k + 1, hypot(re[k], im[k]));
    }
    return h;
}

/// Convert FFT spectrum to gain-space histogram
TH1D* makeFFTgainHist(const TH1D* hFFT, int Npad, double adc_bin_width,
                      double gain_min, double gain_max, const char* name) {
    TH1D* hG = new TH1D(name, "FFT spectrum vs gain;gain [ADC/p.e.];amplitude",
                        250, gain_min, gain_max);
    hG->SetDirectory(nullptr);
    if (!hFFT) return hG;

    for (int bin = 2; bin <= hFFT->GetNbinsX(); ++bin) {
        const int q = bin - 1;
        double a = hFFT->GetBinContent(bin);
        double g = (double)Npad * adc_bin_width / (double)q;
        if (g < gain_min || g > gain_max) continue;
        int bg = hG->FindBin(g);
        hG->SetBinContent(bg, hG->GetBinContent(bg) + a);
    }
    return hG;
}

/// Gain error estimation (simplified)
inline double gainErr_simple(double Npad, double adc_bin_width, double q0) {
    const double sigq = 0.5 / std::sqrt(12.0);
    return (q0 > 0 ? (Npad * adc_bin_width * sigq) / (q0 * q0) : -1.0);
}

// ============ Peak Detection ============

struct KPick {
    int kp = -1;
    double k0 = -1;
    double gain = -1;
    double gainErr = -1;
    double snr = -1;
    double peakAmp = -1;
    double noiseMed = -1;
};

/// Find FFT peak and extract gain + SNR
KPick pickK0(const TH1D* h, int Npad, double adc_bin_width,
             double gain_min, double gain_max, double edge_gain_tol = 0.01,
             int kedge_skip = 4) {
    KPick r;
    if (!h || Npad <= 0 || adc_bin_width <= 0 || gain_min <= 0 ||
        gain_max <= gain_min) {
        return r;
    }

    const double adc_span = (double)Npad * adc_bin_width;
    int kmin = std::max(1, (int)std::floor(adc_span / gain_max));
    int kmax = std::min(h->GetNbinsX() - 1, (int)std::ceil(adc_span / gain_min));
    if (kmax < kmin) return r;

    auto amp = [&](int q) {
        return (q >= 0 && q < h->GetNbinsX()) ? h->GetBinContent(q + 1) : 0.0;
    };

    auto argmax = [&](int a, int b) {
        int kp = a;
        double pk = amp(kp);
        for (int k = a + 1; k <= b; ++k) {
            double v = amp(k);
            if (v > pk) {
                pk = v;
                kp = k;
            }
        }
        return kp;
    };

    auto bestLocalMax = [&](int a, int b) {
        int best = -1;
        double bestA = -1;
        a = std::max(a, 2);
        b = std::min(b, h->GetNbinsX() - 2);
        for (int k = a; k <= b; ++k) {
            double c = amp(k);
            if (c > amp(k - 1) && c >= amp(k + 1)) {
                if (c > bestA) {
                    bestA = c;
                    best = k;
                }
            }
        }
        return best;
    };

    int kp = argmax(kmin, kmax);
    double pk = amp(kp);

    double gain0 = (kp > 0 ? adc_span / kp : -1);
    bool edgeLike = (kp <= kmin + 1) || 
                    (gain0 >= gain_max * (1.0 - edge_gain_tol));

    if (edgeLike) {
        int kstart = std::min(kmin + kedge_skip, kmax);
        int kp2 = bestLocalMax(kstart, kmax);
        if (kp2 < 0) kp2 = argmax(kstart, kmax);
        if (kp2 > 0 && amp(kp2) > 0) {
            kp = kp2;
            pk = amp(kp);
        }
    }

    // Amplitude-weighted interpolation around the peak
    double w = 0, kw = 0;
    for (int dk = -1; dk <= +1; ++dk) {
        int kk = kp + dk;
        if (kk < 1 || kk >= h->GetNbinsX()) continue;
        double a = amp(kk);
        w += a;
        kw += kk * a;
    }
    double k0 = (w > 0 ? kw / w : kp);

    // Compute noise median (avoid peak region)
    std::vector<double> noise;
    noise.reserve((size_t)(kmax - kmin + 1));
    for (int k = kmin; k <= kmax; ++k) {
        if (std::abs(k - kp) <= 2) continue;
        noise.push_back(amp(k));
    }
    double med = median(noise);

    r.kp = kp;
    r.k0 = k0;
    r.gain = (k0 > 0 ? adc_span / k0 : -1);
    r.gainErr = (k0 > 0 ? gainErr_simple((double)Npad, adc_bin_width, k0) : -1);
    r.peakAmp = pk;
    r.noiseMed = med;
    r.snr = (med > 0 ? pk / med : -1);
    return r;
}

} // namespace AHCALRecoAlg
