# FCChh_HiggLFV

### How to run

```bash
root -l -b -q 'analyze_pipeline.cpp("PATH_TO_ROOT_FILE_OR_PATTERN","OUTPUT_ROOT_FILE","CONFIG_JSON_FILE")'
```
- `PATH_TO_ROOT_FILE_OR_PATTERN`: e.g. `samples/HMuTauE_LFV_125.root` or `samples/` (to process all ROOT files in the directory)
- `OUTPUT_ROOT_FILE`: e.g. `out_hist_cfg.root` contains the histograms after applying the pipeline selections
- `CONFIG_JSON_FILE`: e.g. `pipeline.json` defines the selection pipeline and parameters

### How to add new selections

1) Implement the new selection class in `selections.cpp` (e.g. `HToMuTauESelection`), following the structure of existing selections.
2) Add the new selection class declaration in `selection.h`.
3) In `pipeline_config.cpp`, add a case in the selection factory function to create an instance of your new selection class when the corresponding name is encountered in the config.
4) In `pipeline.json`, add a new entry in the "selections" list with the name that matches the case you added in `pipeline_config.cpp` and specify any parameters needed for that selection.

### How to add new parameters

1) Define the new parameter in `type.h` (e.g. `double my_new_param;`) under the `Parameters` struct.
2) In `pipeline_config.cpp`, update the code that reads the config file to parse the new parameter and store it in the `Parameters` struct.
3) In `pipeline.json`, add the new parameter under the "parameters" section with an appropriate value.
4) (Optional) Add print statements in `analyze_pipeline.cpp` to verify that the new parameter is being read correctly during execution.


### Metadata
The `Meta` struct is used to store intermediate information about the event that can be accessed by different selections in the pipeline. For example, it can store the indices of the selected leptons, the flavor of the Z boson leptons, and any kinematic variables that are calculated during the selections. 

#### How to use Meta
1) Define new variables in the `Meta` struct in `type.h` (e.g. `double deltaR_mu_e;`).
2) In your selection class (e.g. `HToMuTauESelection`), calculate the desired variable and store it in the `Meta` struct (e.g. `meta.deltaR_mu_e = deltaR_mu_e;`).
3) In subsequent selections in the pipeline, you can access this variable from the `Meta` struct (e.g. `double deltaR = meta.deltaR_mu_e;`).