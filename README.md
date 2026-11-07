# RestingState-Prep

MATLAB/EEGLAB pipeline for resting-state EEG. It covers import into BIDS with event
reconstruction, preprocessing, spectral analysis and sensitivity analyses.

## Project

The dataset contains resting-state EEG from 37 healthy adults, recorded with a BioSemi
ActiveTwo system (128 scalp channels plus 8 external channels). Each recording has 4
eyes-open and 4 eyes-closed periods of about 63 s each, presented in one of two
counterbalanced sequences. The pipeline organises the data in BIDS and preprocesses
them for spectral measures: individual alpha frequency, theta, alpha and beta power,
band ratios, the aperiodic (1/f) component and the Berger effect. It then tests how
robust the results are.

## Pipeline

Run the scripts in this order. Each one reads the output of the previous step.

1. **`import_bdf_to_bids.m`**: BDF import, event reconstruction and BIDS conversion.
   *Input:* raw BioSemi `.bdf` files and PsychoPy logs. *Output:* a BIDS raw dataset, plus
   the base EEGLAB datasets in `derivatives/eeglab-import`.
2. **`preprocess_rs_eeg.m`**: preprocessing (filtering, bad channels, electrode bridges,
   ICA with ICLabel, interpolation, re-referencing, epoching).
   *Input:* `derivatives/eeglab-import`. *Output:* `derivatives/preproc-spectral`.
3. **`spectral_analysis_rs.m`** (version 2.1): spectral analysis (IAF, band power,
   ratios, aperiodic component, Berger effect, per-period values).
   *Input:* `derivatives/preproc-spectral`. *Output:* `derivatives/spectral`.
4. **`sensitivity_rs.m`**: sensitivity analyses (screening, influence of single subjects,
   covariates, subject and analysis multiverse, case studies).
   *Input:* `derivatives/spectral`. *Output:* `derivatives/spectral/sensitivity`.

**Quality control.** `qc/qc_dataset_overview.m` is a read-only quality check of the EEG
datasets (bad and bridged channels, saturation, events, spectra). It can be run on the
data before or after preprocessing.

`resources/128BIOSEMI_EXG.xyz` contains the template electrode positions (128 scalp and
8 external channels).

## Requirements

- MATLAB R2026a
- EEGLAB with the clean_rawdata and ICLabel plugins
- Signal Processing Toolbox
- Step 4 (`sensitivity_rs.m`) needs base MATLAB only.
- The quality-control script also uses `kurtosis` and `corr` from the Statistics and
  Machine Learning Toolbox.

## How to run

Each script has a **CONFIGURATION** section at the top. Set the `base` path there to the
folder that contains your data, review the other parameters, then run the script in
MATLAB. The quality-control script sets its input and output folders directly
(`cfg.input_folder`, `cfg.output_folder`) instead of using `base`.

## Data and outputs

**No data are included in this repository**, to protect participant privacy. This
covers EEG recordings, PsychoPy logs, participant tables and results. The `.gitignore`
file excludes them.

Every script writes a log, the configuration it used, and README files that describe
every output table.

## How to cite

*To be completed.*

## License

This code is released under the GNU General Public License v3.0. See [LICENSE](LICENSE).
