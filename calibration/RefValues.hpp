#pragma once
namespace AHCALRefValues {
    const double MIP_E=0.461;//MeV
    const double SwitchPoint=500;
    const double ref_ped_highgain=377.8;
    const double ref_ped_highgain_error=26.0;
    const double ref_ped_lowgain=373.8;
    const double ref_ped_lowgain_error=26.0;
    const double ref_MIP=344.3;
    const double ref_gain_ratio=30;
    const double ref_gain_ratio_error=1.0;
    const int lowgain_plat=2000;

    enum HGPedestalStatus {
        HGPedestal_OK = 0,
        HGPedestal_entries_not_enough = -999,
        HGPedestal_cut_by_rms_per_sigma = -998,
        HGPedestal_cut_by_sigma_vs_chi2 = -997,
        HGPedestal_cut_by_both_rms_per_sigma_and_sigma_vs_chi2 = -996,
        HGPedestal_cut_by_max_sigma_and_max_rms = -995,
        HGPedestal_cut_by_both_max_sigma_and_max_rms_and_rms_per_sigma = -994,
        HGPedestal_fit_status_covariance_matrix_not_pos_def_but_modified = 1,
        HGPedestal_fit_status_invalid_hesse = 2,
        HGPedestal_fit_status_EDM_too_large = 3,
        HGPedestal_fit_status_above_max_iterations = 4,
        HGPedestal_fit_status_covariance_matrix_not_pos_def = 5
    };
    enum LGPedestalStatus{
        LGPedestal_OK = 0,
        LGPedestal_entries_not_enough = -999,
        LGPedestal_cut_by_max_sigma = -998,
        LGPedestal_fit_status_covariance_matrix_not_pos_def_but_modified = 1,
        LGPedestal_fit_status_invalid_hesse = 2,
        LGPedestal_fit_status_EDM_too_large = 3,
        LGPedestal_fit_status_above_max_iterations = 4,
        LGPedestal_fit_status_covariance_matrix_not_pos_def = 5
    };

    inline bool HGPedestalStatus_is_ok(int status) {
        return status == HGPedestal_OK || status == HGPedestal_fit_status_covariance_matrix_not_pos_def_but_modified;
    }
    inline bool LGPedestalStatus_is_ok(int status) {
        return status == LGPedestal_OK || status == LGPedestal_fit_status_covariance_matrix_not_pos_def_but_modified;
    }
}