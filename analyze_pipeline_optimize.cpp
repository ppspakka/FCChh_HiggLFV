#include <stdio.h>
#include <stdlib.h>
#include <memory>
#include <string>
#include <vector>
#include <functional>
#include <map>
#include <limits>
#include <cmath>

#include "TRef.h"
#include "TRefArray.h"
#include "TLorentzVector.h"
#include "Delphes.C"
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TString.h>
#include <TSystem.h>
#include <fstream>


using namespace std;

#include "include/types.h"
#include "include/utils.h"
#include "include/histman.h"
#include "include/selections.h"
#include "include/pipeline_config.h"

using namespace hlfv;

// Implementation units (to make ROOT macro self-contained without a separate build)
#include "src/utils.cpp"
#include "src/histman.cpp"
#include "src/selections.cpp"
#include "src/pipeline_config.cpp"

// ------------------------------------------------------------
// Main analysis entry point (ROOT-callable)
// ------------------------------------------------------------
void analyze_pipeline_optimize(const char* inputPath = "samples/HMuTauE_LFV_125.root",
                      const char* outputPath = "out_hist_demo.root",
                      const char* configPath = "pipeline.json")
{
    gErrorIgnoreLevel = kFatal;

    // Build input chain
    TChain* chain = new TChain("Delphes");
    vector<TString> files = getFileList(inputPath);
    if (files.empty()) {
        printf("No input files found matching pattern: %s\n", inputPath);
        return;
    }
    for (auto& f : files) chain->Add(f);

    // print list of files
    printf("Input files:\n");
    for (auto& f : files) {
        printf("  %s\n", f.Data());
    }
    printf("Total files: %zu\n", files.size());

    // Bind Delphes class and enable only needed branches
    Delphes* delphes = new Delphes(chain);
    chain->SetBranchStatus("*", 0);
    chain->SetBranchStatus("Muon*", 1);
    chain->SetBranchStatus("Electron*", 1);
    chain->SetBranchStatus("MissingET*", 1);
    // Enable Event Number for debugging
    chain->SetBranchStatus("Event", 1);
    chain->SetBranchStatus("Jet*", 1);
    
    // Load pipeline config (JSON). Fallback to defaults if missing.
    PipelineConfig cfg;
    bool haveCfg = loadPipelineConfig(configPath, cfg);
    
    // Build active selections from config
    std::vector<std::pair<std::string, std::unique_ptr<ISelection>>> selections;
    selections.reserve(cfg.selections.size());
    for (const auto& s : cfg.selections) {
        if (!s.enabled) continue; // skip disabled
        auto ptr = makeSelectionByName(s.name);
        if (ptr) selections.emplace_back(s.name, std::move(ptr));
        else fprintf(stderr, "Unknown selection name in config: %s\n", s.name.c_str());
    }
    if (selections.empty()) {
        fprintf(stderr, "No enabled selections. Nothing to do.\n");
        return;
    }

    // Steps naming: include initial (preselection) as step 0
    std::vector<std::string> stepNames;
    stepNames.push_back("00_Initial");
    for (size_t i = 0; i < selections.size(); ++i) {
        char buf[64];
        snprintf(buf, sizeof(buf), "%02zu_%s", i+1, selections[i].first.c_str());
        stepNames.emplace_back(buf);
    }
    
    // Variables registry (1D & 2D)
    std::vector<HistogramManager::VarSpec> variables = Variables::getDefault();
    std::vector<HistogramManager::Var2DSpec> variables2D = Variables::getDefault2D();
    
    HistogramManager hman(stepNames, variables, variables2D);
    Parameters params = cfg.params; // from config or defaults
    
    // Print all Parameters used (in loop)
    printf("==== Analysis Parameters ====\n");
    for (const auto& p : {
        std::make_pair("mode", params.mode),
        std::make_pair("mu_pt_min", params.mu_pt_min),
        std::make_pair("e_pt_min", params.e_pt_min),
        std::make_pair("min_dr_e_mu", params.min_dr_e_mu),
        std::make_pair("min_dphi_e_mu", params.min_dphi_e_mu),
        std::make_pair("jet_pt_min", params.jet_pt_min),
        std::make_pair("n_jet", params.n_jet),
        std::make_pair("e_pt_cut", params.e_pt_cut),
        std::make_pair("mu_pt_cut", params.mu_pt_cut),
        std::make_pair("max_dphi_lep_met", params.max_dphi_lep_met)
    }) {
        printf("%-20s : %g\n", p.first, p.second);
    }
    printf("=============================\n");
    

    // --- Cache Setup ---
    // TString cacheFile = TString(outputPath).ReplaceAll(".root", "_passed.txt");
    // Use the inputPath (IS DIR NOT THE .root FILE) and append passed_event.txt to build cache path.
    TString cacheFile;
    if (TString(inputPath).EndsWith(".root")) {
        cacheFile = TString(inputPath).ReplaceAll(".root", "_passed.txt");
    } else { // Is directory, build cache file inside it
        if (TString(inputPath).EndsWith("/")) {
            cacheFile = TString(inputPath) + "passed_events.txt";
        } else {
            cacheFile = TString(inputPath) + "/passed_events.txt";
        }
    }
    std::unordered_set<Long64_t> passedEvents;
    bool useCache = false;
    std::ifstream fin(cacheFile.Data());
    if (fin.is_open()) {
        useCache = true;
        Long64_t ev;
        
        // This will now safely read integers. 
        // Make sure your file ONLY contains numbers!
        while (fin >> ev) {
            passedEvents.insert(ev);
        }
        printf("Loaded %zu passed events from cache.\n", passedEvents.size());
        fin.close(); // Good practice to close when done reading
    } else {
        printf("No cache found. Running full selection and building cache...\n");
    }
    
    std::ofstream pfout;
    if (!useCache) {
        pfout.open(cacheFile.Data());
    }
    const size_t target_cache_step = 2; // e.g., cache after lepton(0) and jet(1)
    // Event loop
    const Long64_t nEntries = chain->GetEntries();
    // Dynamic cutflow: index 0 = Initial, then one per selection
    std::vector<Long64_t> cutflow(stepNames.size(), 0);

    Event evt{delphes};
    Meta meta;

    // for (Long64_t i = 0; i < nEntries; ++i) {
    //     delphes->GetEntry(i);

    //     // Step 0: initial
    //     hman.fill(0, evt, meta, 1.0);
    //     cutflow[0]++;

    //     Meta mcur = meta;
    //     bool passAll = true;
    for (Long64_t i = 0; i < nEntries; ++i) {
        
        // 1. ALWAYS increment the initial cutflow counter. 
        // This guarantees your "Total" event count remains exactly nEntries.
        cutflow[0]++; 

        // 2. Fast Bypass Check
        if (useCache && passedEvents.find(i) == passedEvents.end()) {
            // Event failed in a previous run. 
            // Skip reading branch data and processing selections.
            continue; 
        }

        // 3. For cached events (or all events on the first run), read the real data
        delphes->GetEntry(i);

        // Fill Step 0 kinematic histograms ONLY with real data
        hman.fill(0, evt, meta, 1.0);

        Meta mcur = meta;
        bool passAll = true;
        for (size_t si = 0; si < selections.size(); ++si) {
            Meta next = mcur;
            if (!selections[si].second->apply(evt, next, params)) { passAll = false; break; }
            hman.fill(si+1, evt, next, 1.0);
            cutflow[si+1]++;
            mcur = std::move(next);
            // --- WRITE TO CACHE ---
            // If it passed our target step (e.g., si == 1 is the 2nd cut) and we are building cache
            if (!useCache && si == target_cache_step - 1) {
                pfout << i << "\n";
            }
        }
    }
    // close cache file if we opened it
    if (pfout.is_open()) {
        pfout.close();
        printf("Cache file written: %s\n", cacheFile.Data());
    }

    // Write outputs
    TFile* fout = TFile::Open(outputPath, "RECREATE");
    if (!fout || fout->IsZombie()) {
        printf("Failed to create output file: %s\n", outputPath);
        return;
    }
    hman.writeAll(fout);

    // Small cutflow summary as a TNamed or histogram
    TH1F cutflowH("cutflow", "Cutflow;Step;Events", (int)cutflow.size(), 0.5, (double)cutflow.size()+0.5);
    for (size_t i = 0; i < cutflow.size(); ++i) {
        cutflowH.GetXaxis()->SetBinLabel((int)i+1, stepNames[i].c_str());
        cutflowH.SetBinContent((int)i+1, cutflow[i]);
    }
    cutflowH.Write();

    fout->Close();
    // for efficiency calculation
    int initial = cutflow[0];
    printf("\n==== Pipeline summary ====\n");
    printf("Total events:          %lld\n", nEntries);
    for (size_t i = 0; i < cutflow.size(); ++i) {
        // calculate efficiency
        double efficiency = static_cast<double>(cutflow[i]) / initial;
        printf("After %-30s %lld (Efficiency: %.6f%%)\n", stepNames[i].c_str(), cutflow[i], efficiency * 100);
    }
    printf("Output written to: %s\n", outputPath);
    printf("===============================\n");
    gSystem->Exit(0); // ensure ROOT macro exits cleanly
}
