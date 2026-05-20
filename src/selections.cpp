#include "../include/selections.h"

#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>
#include <TLorentzVector.h>
#include <TVector3.h>

namespace hlfv {

// --- Constants ---
constexpr double MASS_E   = 0.000511;
constexpr double MASS_MU  = 0.10565837;
constexpr double ECM      = 240.0; // Center of mass energy
constexpr double Z_MASS   = 91.1876;

enum Flavor { ELECTRON = 0, MUON = 1 };

struct LeptonObj {
    int index;
    int flavor; 
    double pt;
    int charge;
};

// --- Helper Functions ---
namespace {

    // Generic P4 retriever
    TLorentzVector get_p4(const Event& evt, int idx, int flavor) {
        TLorentzVector v;
        if (flavor == ELECTRON) {
            v.SetPtEtaPhiM(evt.d->Electron_PT[idx], evt.d->Electron_Eta[idx], 
                           evt.d->Electron_Phi[idx], MASS_E);
        } else {
            v.SetPtEtaPhiM(evt.d->Muon_PT[idx], evt.d->Muon_Eta[idx], 
                           evt.d->Muon_Phi[idx], MASS_MU);
        }
        return v;
    }

    double delta_phi(double phi1, double phi2) {
        double dphi = std::abs(phi1 - phi2);
        if (dphi > M_PI) dphi = 2 * M_PI - dphi;
        return dphi;
    }

    // Generic Collinear Mass Calculation
    // Calculates mass of (visible_lep + invisible_tau_products)
    // Assumes neutrinos go in direction of visible tau product (approx)
    double compute_collinear_mass(const Event& evt, const TLorentzVector& v_tau_vis, const TLorentzVector& v_other) {
        if (!evt.d || evt.d->MissingET_size <= 0) return std::numeric_limits<double>::quiet_NaN();
        if (v_tau_vis.Pt() <= 0) return std::numeric_limits<double>::quiet_NaN();

        double met = evt.d->MissingET_MET[0];
        double met_phi = evt.d->MissingET_Phi[0];
        double met_px = met * std::cos(met_phi);
        double met_py = met * std::sin(met_phi);

        // Project MET onto the visible tau direction
        double dir_x = v_tau_vis.Px() / v_tau_vis.Pt();
        double dir_y = v_tau_vis.Py() / v_tau_vis.Pt();
        double proj = met_px * dir_x + met_py * dir_y;

        if (proj <= 0) return std::numeric_limits<double>::quiet_NaN();

        // Visible fraction of tau momentum
        double x_tau_vis = v_tau_vis.Pt() / (v_tau_vis.Pt() + proj);

        if (x_tau_vis <= 0 || x_tau_vis > 1) return std::numeric_limits<double>::quiet_NaN();

        double m_vis = (v_tau_vis + v_other).M();
        return m_vis / std::sqrt(x_tau_vis);
    }

    double compute_transverse_mass(double pt, double phi, double met, double met_phi) {
        double dphi = delta_phi(phi, met_phi);
        return std::sqrt(2 * pt * met * (1 - std::cos(dphi)));
    }
    
} // anonymous namespace

// --- Selections ---

std::string EmptySelection::name() const { return "EmptySelection"; }
bool EmptySelection::apply(const Event& evt, Meta& meta, const Parameters& cfg) {
    // Example test
    return true;
}

std::string FinalState_NoCut::name() const { return "FinalState_NoCut"; }
bool FinalState_NoCut::apply(const Event& evt, Meta& meta, const Parameters& cfg) {
    return true;
}


// TODO: Plan the appropriate order of variable plots, as this FCC-hh run extremely slow.
// Show collinear mass distribution at the very first cut? as we capable to do so.
// Context: cfg.mode: 1=etaumu, 0=mutaue; use for collinear mass calculation.

// FCC-hh selection
// 1. Lepton veto
//     1. Filter muon with pT > 10 and electron with pT > 10, all |eta| < 6.0
//     2. require exactly one pair with charge conservation
//     3. dR > 0.3
//     4. DeltaPhi > 2.2
std::string LeptonSelection::name() const { return "LeptonSelection"; }
bool LeptonSelection::apply(const Event& evt, Meta& meta, const Parameters& cfg) {
    if (!evt.d) return false;

    // save met
    if (evt.d->MissingET_size > 0) {
        meta.MET = evt.d->MissingET_MET[0];
    }

    // Helper to extract passing leptons    
    auto extract = [&](int n, const float* pts, const float* etas, const float* phis, const int* charges, int flav, float pt_min) {
        std::vector<LeptonObj> leps;
        for (int i = 0; i < n; ++i) {
            if (pts[i] > pt_min && std::abs(etas[i]) < 6.0) {
                leps.push_back({i, flav, (double)pts[i], charges[i]});
            }
        }
        return leps;
    };
    // Retrive pt from config; cfg.e_pt_min and cfg.mu_pt_min
    auto electrons = extract(evt.d->Electron_size, evt.d->Electron_PT, evt.d->Electron_Eta, evt.d->Electron_Phi, evt.d->Electron_Charge, ELECTRON, cfg.e_pt_min);
    auto muons     = extract(evt.d->Muon_size, evt.d->Muon_PT, evt.d->Muon_Eta, evt.d->Muon_Phi, evt.d->Muon_Charge, MUON, cfg.mu_pt_min);

    // Exactly 2 leptons
    bool is_1e_1mu = (electrons.size() == 1 && muons.size() == 1);
    if (!is_1e_1mu) return false;

    // Net Charge Check
    int net_charge = electrons[0].charge + muons[0].charge;
    if (net_charge != 0) return false;

    // dR Check
    int idx_e = electrons[0].index;
    int idx_mu = muons[0].index;

    TLorentzVector p4_e  = get_p4(evt, idx_e, ELECTRON);
    TLorentzVector p4_mu = get_p4(evt, idx_mu, MUON);
    double deltaR = p4_e.DeltaR(p4_mu);
    if (deltaR <= cfg.min_dr_e_mu) return false;

    // DeltaPhi Check
    double dphi = std::abs(p4_e.DeltaPhi(p4_mu));
    if (dphi <= cfg.min_dphi_e_mu) return false;

    // Assign variables 
    meta.mu_idx = idx_mu;
    meta.e_idx = idx_e;

    meta.e_pt = p4_e.Pt();
    meta.mu_pt = p4_mu.Pt();

    // collinear mass calculation
    TLorentzVector v_tau_vis = (cfg.mode == 0) ? p4_mu : p4_e; // if mode=0 (etaumu), tau_vis is mu; if mode=1 (mutaue), tau_vis is e
    TLorentzVector v_other = (cfg.mode == 0) ? p4_e : p4_mu; // the other lepton
    meta.m_collinear = compute_collinear_mass(evt, v_tau_vis, v_other);


    return true;
}

// 2. Jet cut
//     1. Filter jets with pT > 30 and |eta| < 6.0
//         1. if any passing jet appear to be b-tagging, reject event
// Example bitwise btagging: (below is example to accept, if reject do opposite)
        // for (int j=0; j<indelphes->Jet_size; j++)
        // {
        //     if (indelphes->Jet_PT[j] < 40) continue;
        //     if (TMath::Abs(indelphes->Jet_Eta[j]) > 6.0) continue;
        //     passed_jets.push_back(j);
        //     if (indelphes->Jet_BTag[j] & 0b111) passed_b_jets.push_back(j);
        // }
std::string JetVeto::name() const { return "JetVeto"; }
bool JetVeto::apply(const Event& evt, Meta& meta, const Parameters& cfg) {
    if (!evt.d) return false;

    int jet_count = 0;
    for (int j = 0; j < evt.d->Jet_size; ++j) {
        if (evt.d->Jet_PT[j] > cfg.jet_pt_min && std::abs(evt.d->Jet_Eta[j]) < 6.0) {
            ++jet_count;
            // Check b-tagging; assuming b-tagging info is in Jet_BTag and uses bitwise flags
            if (evt.d->Jet_BTag[j] & 0b111) { // Example: if any of the first three bits are set, it's b-tagged
                return false; // Reject event if a b-tagged jet is found
            }
        }
    }
    // Assign variables (if any)
    meta.n_jet = jet_count;

    // debug, dump all jets info 
    // cout << "Event has " << evt.d->Jet_size << " jets, " << jet_count << " pass the pT and eta cuts." << endl;
    // for (int j = 0; j < evt.d->Jet_size; ++j) {
    //     cout << "Jet " << j << ": PT=" << evt.d->Jet_PT[j] << ", Eta=" << evt.d->Jet_Eta[j] 
    //          << ", BTag=" << evt.d->Jet_BTag[j] << endl;
    // }

    return true; // No jets or only non-b-tagged jets passed the criteria
}

// 3. Jet categorize: divide into two categories: 0 jet and 1 jet
//     1. configured via `n_jet` variable
std::string JetCategorization::name() const { return "JetCategorization"; }
bool JetCategorization::apply(const Event& evt, Meta& meta, const Parameters& cfg) {
    if (!evt.d) return false;

    // debug: print jet count and category info
    // cout << "JetCategorization: Event has " << meta.n_jet << " jets, required category is " << cfg.n_jet << endl;

    if (meta.n_jet != cfg.n_jet) {
        return false; // Reject event if it doesn't match the configured jet category
    }
    return true;
}

// 4. prompt lepton cut
//     1. (LFV h to mu tau_e)
//         1. Low-mass cut: pT_mu > 60
//         2. High-mass cut: pT_mu > 150
//         3. All: pT_e > 10
//     2. (LFV h to e tau_mu)
//         1. Low-mass cut: pT_e > 60
//         2. High-mass cut: pT_e > 150
//         3. All: pT_mu > 10
// NOTE: we dont care the high, low mass category here, only consider input configuration
// cfg.e_pt and cfg.mu_pt can be set to corresponding channels
std::string LeptonPTCut::name() const { return "LeptonPTCut"; }
bool LeptonPTCut::apply(const Event& evt, Meta& meta, const Parameters& cfg) {
    if (!evt.d) return false;

    int idx_e = meta.e_idx;
    int idx_mu = meta.mu_idx;

    double e_pt = evt.d->Electron_PT[idx_e];
    double mu_pt = evt.d->Muon_PT[idx_mu];

    if (e_pt < cfg.e_pt_cut || mu_pt < cfg.mu_pt_cut) {
        return false; // Reject event if either lepton fails the pT cut
    }
    return true; // Both leptons pass the pT cuts
}

// 5. DeltaPhi(MET, tau_vis): 
//     1. (LFV h to mu tau_e): DeltaPhi(MET, e) < 0.7
//     2. (LFV h to e tau_mu): DeltaPhi(MET, mu) < 0.7
// Read from Mode
std::string DeltaPhiMETLepton::name() const { return "DeltaPhiMETLepton"; }
bool DeltaPhiMETLepton::apply(const Event& evt, Meta& meta, const Parameters& cfg) {
    if (!evt.d) return false;
    if (evt.d->MissingET_size <= 0) return false;

    // if mode = 0, tau_vis is mu; if mode = 1, tau_vis is e
    int idx_lep = (cfg.mode == 0) ? meta.mu_idx : meta.e_idx;
    double lep_phi = (cfg.mode == 0) ? evt.d->Muon_Phi[idx_lep] : evt.d->Electron_Phi[idx_lep];
    double met_phi = evt.d->MissingET_Phi[0];

    double dphi = delta_phi(lep_phi, met_phi);

    if (dphi >= cfg.max_dphi_lep_met) {
        return false;
    }
    return true;
}

} // namespace hlfv
