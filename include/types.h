#pragma once

#include <limits>

class Delphes; // forward declaration

namespace hlfv {

struct Parameters {
    // FCC-hh specific parameters
    // name convention: if object: <object>_<variable>_<cuttype>
    // if function (deltaR, deltaPhi): <cuttype>_<func>_<objects>

    double mode; // 1=etaumu, 0=mutaue; use for collinear mass calculation

    // Lepton selection
    double mu_pt_min;
    double e_pt_min;
    double min_dr_e_mu;
    double min_dphi_e_mu;

    // Jet cut
    double jet_pt_min;
    double n_jet;

    // Lepton pT cut
    double mu_pt_cut;
    double e_pt_cut;

    // DeltaPhi cut
    double max_dphi_lep_met;



    
};

struct Event {
    Delphes* d = nullptr; // bound Delphes tree object
};

struct Meta {
    int e_idx = -1;
    int mu_idx = -1;
    double e_pt = std::numeric_limits<double>::quiet_NaN();
    double mu_pt = std::numeric_limits<double>::quiet_NaN();
    double m_collinear = std::numeric_limits<double>::quiet_NaN();
    double MET = std::numeric_limits<double>::quiet_NaN();
    int n_jet = -1;
};

} // namespace hlfv
