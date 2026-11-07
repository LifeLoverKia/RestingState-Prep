% =========================================================================
% qc_dataset_overview_v2.m
% =========================================================================
%
% PURPOSE
%   READ-ONLY initial quality check (QC) of the resting-state EEG dataset
%   (BioSemi ActiveTwo, 128 scalp EEG + 8 external channels, EEGLAB .set
%   files that already contain the hand-corrected triggers 11 = EO and
%   12 = EC).
%
%   The script does NOT modify, filter, re-reference or save any EEG file.
%   All filtering / re-referencing is done on TEMPORARY copies in memory,
%   only to compute descriptive numbers.
%
% WHAT IS NEW IN v2 (compared with qc_initial_dataset_overview.m)
%   1. Electrode bridges are detected with the electrical-distance (ED)
%      method of Alschuler et al. (2014, Clin Neurophysiol 125:484-490,
%      doi:10.1016/j.clinph.2013.08.024), applied to short segments of the
%      continuous data and to ALL channel pairs. Bridged channels are
%      reported SEPARATELY and are no longer counted as "suspect".
%      The script also estimates how many data dimensions are lost because
%      of the bridges (important later for ICA).
%   2. Ocular channels: high kurtosis / high amplitude in channels that are
%      strongly correlated with the VEOG or HEOG is labelled "ocular" and is
%      NOT counted as a bad-channel reason (blinks are expected there).
%   3. Saturated channels (DC offset close to the BioSemi +/-262 mV input
%      range) are detected explicitly.
%   4. Bridged channels are NOT excluded from the posterior alpha check
%      (in v1 this left S17 and S24 without posterior channels).
%   5. Latency of every event other than 11/12 (e.g. the leftover '13') is
%      reported. These events are IGNORED by all computations.
%   6. Trigger times that look hand-placed (round values to 0.01 s) and
%      recordings ending exactly 1.000 s after the last trigger are counted
%      (needed to document the event file when converting to BIDS).
%   7. Figures are saved with a white background (MATLAB dark theme off).
%   v2.1: works on BIDS derivative folders (files in sub-XX/eeg subfolders,
%   subject = 'sub-XX') and detects narrow-band spectral peaks (24-47 Hz,
%   e.g. equipment interference) that could contaminate the beta band.
%
% HOW TO RUN IT ON THE PREPROCESSED DATA (comparison with the initial QC)
%   cfg.input_folder  = '...\BIDS\derivatives\preproc-spectral';
%   cfg.file_pattern  = 'sub-*\eeg\*_desc-clean_eeg.set';
%   cfg.output_folder = '...\QC_after_preproc';
%   Note: in the cleaned data the bridged channels are exact copies of their
%   representative, so they will (correctly) still be detected as bridges.
%
% OUTPUTS (written to cfg.output_folder)
%   QC_v2_report.txt      -> human-readable report, subject by subject,
%                            plus a GROUP SUMMARY at the end
%   QC_v2_summary.csv     -> one row per subject (open in Excel)
%   QC_v2_channels.csv    -> one row per channel per subject
%   QC_v2_epochs.csv      -> one row per trigger (alpha, blinks...)
%   QC_v2_bridges.csv     -> one row per bridged channel pair
%   QC_v2_config.mat      -> the exact settings used (traceability)
%   figures\<subject>_QC_v2.png -> PSD of all channels + posterior
%                                  EO vs EC spectra
%
% REQUIREMENTS
%   - MATLAB (written for R2026a)
%   - Signal Processing Toolbox (pwelch, hann, fir1, filtfilt, findpeaks)
%   - EEGLAB (pop_loadset, pop_select, pop_eegfiltnew = firfilt plugin,
%     readlocs). The script starts EEGLAB without GUI if needed.
%   - natural_sort_files.m on the MATLAB path (your helper function)
%
% IMPORTANT
%   All channel flags are HEURISTIC and INFORMATIONAL. They guide the
%   preprocessing decisions; they are NOT final bad-channel decisions.
%   The bridge detection is a re-implementation of the published method,
%   not the official eBridge code. For a confirmation on single subjects
%   the official eBridge function can be used:
%   https://psychophysiology.cpmc.columbia.edu/software/eBridge/index.html
% =========================================================================

clear; clc;

%% ========================================================================
%  0. CONFIGURATION  (edit only this section)
%  ========================================================================
cfg = struct();

% ------------------------------ Paths ------------------------------------
% Folder with the .set files that contain the corrected triggers
cfg.input_folder  = 'C:\Users\chiar\Documents\MATLAB\FONDECYT\DATOS\Resting State\BIDS\derivatives\preproc-spectral';
cfg.file_pattern  = 'sub-*\eeg\*_desc-clean_eeg.set';
% Where the QC outputs will be written (created if it does not exist).
% A NEW folder, so the v1 results are not overwritten.
cfg.output_folder = 'C:\Users\chiar\Documents\MATLAB\FONDECYT\DATOS\Resting State\QC_after_preproc';
% Channel coordinate file (used ONLY if the .set files lack coordinates)
cfg.coord_file    = 'C:\Users\chiar\Documents\MATLAB\FONDECYT\code\128BIOSEMI_EXG.xyz';
% Save one PNG figure per subject (true/false)
cfg.save_figures  = true;

% ------------------------- Expected recording ----------------------------
cfg.expected_srate = 512;   % Hz (files were resampled from 2048 Hz)
cfg.n_eeg          = 128;   % scalp channels
cfg.n_ext          = 8;     % external channels
% Roles of the 8 external channels, in the documented order EXT1..EXT8
% (in the files they are still called EXG1..EXG8)
cfg.ext_labels = {'HEOG_L','HEOG_R','VEOG_U','VEOG_D', ...
                  'MAST_L','MAST_R','CARD_L','CARD_R'};

% ----------------------------- Paradigm ----------------------------------
cfg.event_eo = '11';        % trigger code: Eyes Open
cfg.event_ec = '12';        % trigger code: Eyes Closed
% Events that are known leftovers and must be ignored (reported only)
cfg.ignore_events = {'13'};
% The two counterbalanced sequences (trigger codes in temporal order)
cfg.seq_C1 = {'11','12','12','11','12','11','11','12'};  % EO EC EC EO EC EO EO EC
cfg.seq_C2 = {'12','11','11','12','11','12','12','11'};  % EC EO EO EC EO EC EC EO
% Recommended epoch window (trigger = END of the condition)
cfg.epoch_win_recommended = [-58.9 0.1];   % s relative to trigger
% Plausible range for the interval between two consecutive triggers
cfg.iti_expected_range_s  = [55 75];       % s
% "Core" window used for the alpha / blink sanity check. It avoids the
% first seconds of each condition (eyes opening/closing transition) and
% the last 2 s before the audio instruction.
cfg.core_win   = [-55 -2];   % s relative to trigger
cfg.min_core_s = 10;         % minimum usable seconds to compute an epoch

% ------------------- Temporary filters (QC only) -------------------------
cfg.qc_hp_hz      = 1;       % high-pass for channel / spectral metrics
cfg.qc_lp_hz      = 40;      % low-pass for time-domain metrics
cfg.welch_win_s   = 4;       % Welch window length (s) -> 0.25 Hz resolution
cfg.welch_overlap = 0.5;     % 50 % overlap
cfg.line_freq     = 50;      % Hz (Chile)
cfg.bands = struct('delta',    [1 4], ...
                   'theta',    [4 8], ...
                   'alpha',    [8 13], ...
                   'beta',     [13 30], ...
                   'lowgamma', [30 45]);

% ------------------- Heuristic channel thresholds ------------------------
cfg.flat_sd_uV     = 0.5;    % SD below this (after 1 Hz HP) = flat
cfg.max_flat_run_s = 1;      % identical consecutive samples longer than this
cfg.saturation_uV  = 250000; % |DC offset| above this = saturated channel
                             % (BioSemi ActiveTwo input range is +/-262 mV)
cfg.offset_thr_uV  = 40000;  % |DC offset| above this = large offset
cfg.robust_z_thr   = 5;      % robust z threshold (median / MAD) across channels
cfg.neigh_r_thr    = 0.4;    % median corr. with neighbours below this = suspect
cfg.n_neighbours   = 6;      % nearest neighbours (3-D distance)
cfg.amp_thr_uV     = 100;    % amplitude threshold for 1-s windows (1-40 Hz data)
cfg.pct_win_thr    = 20;     % % of 1-s windows above amp_thr -> suspect channel
cfg.line_db_min    = 10;     % line noise flagged only if also > 10 dB
% Narrow-band peaks (non-physiological, e.g. equipment interference)
cfg.narrow_peak_range_hz = [24 47];  % searched range (50 Hz excluded)
cfg.narrow_peak_db       = 3;        % minimum height above background
cfg.narrow_peak_width_hz = 1.5;      % maximum width (physiological beta is broader)
cfg.ocular_r_thr   = 0.4;    % |corr| with VEOG or HEOG above this = ocular channel

% ------------------- Bridge detection (Alschuler et al., 2014) -----------
cfg.bridge_seg_s        = 2;     % length of the segments ("epochs") in s
cfg.bridge_bin          = 0.25;  % histogram bin of the scaled ED
cfg.bridge_interp_step  = 0.05;  % interpolation step of the histogram
cfg.bridge_peak_max     = 5;     % near-zero peak must be at scaled ED <= 5
cfg.bridge_cutoff_max   = 10;    % local minimum after it must be <= 10
cfg.bridge_min_seg_frac = 0.5;   % pair bridged if ED <= cutoff in >= 50 % of segments
% Two safeguards added for a 128-channel cap (neighbouring electrodes are
% close, so the normal ED distribution starts near zero):
cfg.bridge_valley_max   = 0.5;   % valley must be <= 50 % of the near-zero peak height
cfg.bridge_rel_ed_max   = 0.05;  % AND median ED / mean channel variance < 0.05
                                 % (i.e. correlation of about 0.975 or higher)

% ----------------------- Trigger sanity check ----------------------------
% Leave empty to select posterior channels automatically from coordinates
% (most posterior 15 % of the cap, excluding the lowest row). Or write the
% labels, e.g. {'A15','A21','A23','A28'}.
cfg.posterior_labels  = {};
cfg.posterior_frac    = 0.15;
cfg.iaf_search_hz     = [7 14];   % search range for the individual alpha peak
cfg.berger_ratio_ok   = 1.2;      % EC/EO alpha ratio considered "clear"
cfg.blink_band_hz     = [1 10];   % VEOG band-pass for blink detection
cfg.blink_min_uV      = 50;       % minimum blink amplitude
cfg.blink_mad_k       = 5;        % threshold = median + k * robust SD
cfg.blink_min_dist_s  = 0.3;      % minimum distance between blinks
cfg.blink_ratio_ok    = 1.2;      % EO/EC blink-rate ratio considered "clear"
cfg.ecg_band_hz       = [5 30];   % ECG band-pass for R-peak detection
cfg.hr_plausible_bpm  = [40 120];

% ------------------------ Group summary ----------------------------------
% A channel is listed in the group summary if it is suspect / bridged in
% at least this fraction of subjects
cfg.group_chan_frac = 0.20;

%% ========================================================================
%  1. SETUP AND CHECKS
%  ========================================================================

% --- EEGLAB: start it without GUI if its functions are not on the path ---
if ~exist('pop_eegfiltnew', 'file')
    if exist('eeglab', 'file')
        eeglab('nogui');   % adds EEGLAB plugins (firfilt, etc.) to the path
    else
        error('EEGLAB not found on the MATLAB path. Add the EEGLAB folder first.');
    end
end
if ~exist('pop_loadset', 'file') || ~exist('pop_eegfiltnew', 'file')
    error('EEGLAB functions (pop_loadset / pop_eegfiltnew) not available.');
end

% --- Signal Processing Toolbox ---
needed = {'pwelch', 'hann', 'fir1', 'filtfilt', 'findpeaks'};
for k = 1:numel(needed)
    if ~exist(needed{k}, 'file')
        error('Function %s not found: the Signal Processing Toolbox is required.', needed{k});
    end
end

% --- Your natural sorting helper ---
if ~exist('natural_sort_files', 'file')
    error('natural_sort_files.m is not on the MATLAB path.');
end

% --- Folders ---
if ~isfolder(cfg.input_folder)
    error('Input folder not found:\n  %s', cfg.input_folder);
end
if ~isfolder(cfg.output_folder), mkdir(cfg.output_folder); end
fig_folder = fullfile(cfg.output_folder, 'figures');
if cfg.save_figures && ~isfolder(fig_folder), mkdir(fig_folder); end

% --- List the files in natural order (S1, S2, ..., S10) ---
files = dir(fullfile(cfg.input_folder, cfg.file_pattern));
files = natural_sort_files(files);
nS = numel(files);
if nS == 0
    error('No files matching %s in:\n  %s', cfg.file_pattern, cfg.input_folder);
end
fprintf('Found %d file(s) to check.\n', nS);

% --- Output file names ---
report_txt  = fullfile(cfg.output_folder, 'QC_v2_report.txt');
summary_csv = fullfile(cfg.output_folder, 'QC_v2_summary.csv');
chan_csv    = fullfile(cfg.output_folder, 'QC_v2_channels.csv');
epoch_csv   = fullfile(cfg.output_folder, 'QC_v2_epochs.csv');
bridge_csv  = fullfile(cfg.output_folder, 'QC_v2_bridges.csv');

% Save the configuration for traceability
save(fullfile(cfg.output_folder, 'QC_v2_config.mat'), 'cfg');

% --- Open the report ('w' = overwrite, so re-running never duplicates) ---
fid = fopen(report_txt, 'w');
if fid < 0
    error('Cannot write the report file:\n  %s', report_txt);
end
write_report_header(fid, cfg, nS);

%% ========================================================================
%  2. MAIN LOOP OVER SUBJECTS
%  ========================================================================
template   = init_summary_row(cfg);         % empty row with all fields
rows       = repmat(template, nS, 1);       % one row per subject
chan_all   = cell(nS, 1);                   % per-channel tables
epoch_all  = cell(nS, 1);                   % per-trigger tables
bridge_all = cell(nS, 1);                   % per-bridged-pair tables

for iS = 1:nS

    fname   = files(iS).name;
    % Subject ID = file name without '_withEvents' and without '.set'
    subj_id = regexprep(fname, '(_withEvents)?\.set$', '', 'ignorecase');
    if contains(subj_id, '_task')                  % BIDS names: keep 'sub-XX'
        subj_id = extractBefore(subj_id, '_task');
    end

    fprintf('\n[%d/%d] %s\n', iS, nS, fname);

    row         = template;
    row.subject = subj_id;
    row.file    = fname;
    flags       = {};
    chanT       = table();
    epT         = table();
    brT         = table();

    try
        % Load the dataset (read-only: we never call pop_saveset)
        EEG = pop_loadset('filename', fname, 'filepath', files(iS).folder);

        % Run all the QC computations for this subject
        [row, flags, chanT, epT, brT] = qc_one_subject(EEG, row, cfg, fig_folder);
        row.status = "OK";

    catch ME
        % If anything fails, keep going with the next subject
        row.status = "ERROR";
        msg = ME.message;
        if ~isempty(ME.stack)
            msg = sprintf('%s (in %s, line %d)', msg, ME.stack(1).name, ME.stack(1).line);
        end
        row.error_message = msg;
        fprintf(2, '  ERROR: %s\n', msg);
    end

    % Store the flags inside the row and normalise the field types
    row.n_flags    = numel(flags);
    row.flags_text = strjoin(flags, ' | ');
    row            = enforce_fields(row, template);

    rows(iS)       = row;
    chan_all{iS}   = chanT;
    epoch_all{iS}  = epT;
    bridge_all{iS} = brT;

    % Write this subject's block in the TXT report
    write_subject_block(fid, row, flags, chanT, epT, cfg);

    clear EEG
end

%% ========================================================================
%  3. TABLES, GROUP SUMMARY AND SAVING
%  ========================================================================
summaryT = struct2table(rows);

chanAll   = stack_tables(chan_all);
epochAll  = stack_tables(epoch_all);
bridgeAll = stack_tables(bridge_all);

write_group_summary(fid, summaryT, chanAll, cfg);
fclose(fid);

writetable(summaryT, summary_csv);
if ~isempty(chanAll),   writetable(chanAll,   chan_csv);   end
if ~isempty(epochAll),  writetable(epochAll,  epoch_csv);  end
if ~isempty(bridgeAll), writetable(bridgeAll, bridge_csv); end

fprintf('\n%s\n', repmat('=', 1, 70));
fprintf('QC v2 finished. Outputs in:\n  %s\n', cfg.output_folder);
fprintf('%s\n', repmat('=', 1, 70));


%% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================

% -------------------------------------------------------------------------
function [row, flags, chanT, epT, brT] = qc_one_subject(EEG, row, cfg, fig_folder)
% QC_ONE_SUBJECT  Compute all QC metrics for one continuous dataset.

flags = {};

% ================= 1. BASIC FILE PROPERTIES ==============================
if EEG.trials > 1
    error('Dataset is epoched (%d trials); this QC expects continuous data.', EEG.trials);
end

srate = EEG.srate;
X     = double(EEG.data);        % channels x samples, in microvolts
[nCh, N] = size(X);
dur   = N / srate;               % duration in seconds

row.srate_hz   = srate;
row.n_channels = nCh;
row.n_samples  = N;
row.duration_s = dur;
row.data_class = class(EEG.data);
row.ref_field  = ref2str(EEG.ref);
row.has_ica    = ~isempty(EEG.icaweights);

% NaN / Inf values would break all computations: count them and set to 0
bad_val = ~isfinite(X);
row.n_nan_inf = nnz(bad_val);
if row.n_nan_inf > 0
    flags{end+1} = sprintf('%d NaN/Inf samples found (set to 0 for this QC only).', row.n_nan_inf);
    X(bad_val) = 0;
end
clear bad_val

if srate ~= cfg.expected_srate
    flags{end+1} = sprintf('Sampling rate is %g Hz (expected %g Hz).', srate, cfg.expected_srate);
end
if nCh ~= cfg.n_eeg + cfg.n_ext
    flags{end+1} = sprintf('%d channels found (expected %d).', nCh, cfg.n_eeg + cfg.n_ext);
end
if row.has_ica
    flags{end+1} = 'File already contains ICA weights (unexpected for raw data).';
end

% Remove ICA fields from the temporary copy (they are not needed here)
EEG.icaweights = []; EEG.icasphere = []; EEG.icawinv = [];
EEG.icaact = []; EEG.icachansind = [];

% ================= 2. CHANNELS AND COORDINATES ===========================
if isempty(EEG.chanlocs)
    labels = arrayfun(@(k) sprintf('Ch%d', k), 1:nCh, 'UniformOutput', false);
else
    labels = {EEG.chanlocs.labels};
end
labels = labels(:)';
types  = get_chan_types(EEG.chanlocs, nCh);

[eeg_idx, ext_idx, role, mapping] = identify_channels(labels, types, nCh, cfg);
nE = numel(eeg_idx);
eeg_labels = labels(eeg_idx);

row.n_eeg       = nE;
row.n_ext       = numel(ext_idx);
row.ext_mapping = mapping;
if nE < 10
    error('Only %d EEG channels identified: cannot run the QC.', nE);
end

% Coordinates (needed for neighbours and posterior channel selection)
xyz = get_xyz(EEG.chanlocs, eeg_idx);
row.n_eeg_valid_coords = sum(all(isfinite(xyz), 2));
if row.n_eeg_valid_coords < nE
    flags{end+1} = sprintf('Only %d/%d EEG channels have coordinates in the file.', ...
                           row.n_eeg_valid_coords, nE);
    if isfile(cfg.coord_file) && exist('readlocs', 'file')
        locs = readlocs(cfg.coord_file, 'filetype', 'xyz');
        if numel(locs) >= max(eeg_idx)
            xyz = get_xyz(locs, eeg_idx);
            flags{end+1} = 'Coordinates for this QC taken from the .xyz file (assumes the same channel order).';
        end
    end
end

% ================= 3. EVENTS / TRIGGERS ==================================
[ev_types, ev_lat] = get_events(EEG);
row.n_events_total = numel(ev_types);
row.n_boundary     = sum(strcmpi(ev_types, 'boundary'));

is_trig  = ismember(ev_types, {cfg.event_eo, cfg.event_ec});
is_other = ~is_trig & ~strcmpi(ev_types, 'boundary');

% Report type and time of every other event (they are ignored afterwards)
if any(is_other)
    oth_txt = arrayfun(@(i) sprintf('%s@%.2fs', ev_types{i}, (ev_lat(i) - 1) / srate), ...
                       find(is_other), 'UniformOutput', false);
    row.other_events_detail = strjoin(oth_txt, ' ');
    unknown = setdiff(unique(ev_types(is_other)), cfg.ignore_events);
    if ~isempty(unknown)
        flags{end+1} = sprintf('Unexpected event types present: %s', strjoin(unknown, ' '));
    end
end
if row.n_boundary > 0
    flags{end+1} = sprintf('%d boundary event(s): the recording has discontinuities.', row.n_boundary);
end

% Keep the 11/12 triggers, sorted by time
trig_types = ev_types(is_trig);
trig_lat   = ev_lat(is_trig);
[trig_lat, order] = sort(trig_lat(:)');
trig_types = trig_types(order);
trig_types = trig_types(:)';
trig_t     = (trig_lat - 1) / srate;        % seconds from recording start
nT         = numel(trig_t);

row.n_trig_EO        = sum(strcmp(trig_types, cfg.event_eo));
row.n_trig_EC        = sum(strcmp(trig_types, cfg.event_ec));
row.trigger_sequence = strjoin(trig_types, '-');
row.sequence_match   = match_sequence(trig_types, cfg.seq_C1, cfg.seq_C2);
row.n_trig_out_of_bounds = sum(trig_lat < 1 | trig_lat > N | ~isfinite(trig_lat));

if nT == 0
    flags{end+1} = 'NO 11/12 triggers found.';
else
    if row.n_trig_EO ~= 4 || row.n_trig_EC ~= 4
        flags{end+1} = sprintf('Expected 4 EO + 4 EC triggers, found %d EO + %d EC.', ...
                               row.n_trig_EO, row.n_trig_EC);
    end
    if row.sequence_match == "NO MATCH"
        flags{end+1} = sprintf('Trigger sequence (%s) matches neither C1 nor C2.', row.trigger_sequence);
    end
    if row.n_trig_out_of_bounds > 0
        flags{end+1} = sprintf('%d trigger(s) outside the data range.', row.n_trig_out_of_bounds);
    end

    row.first_trigger_s       = trig_t(1);
    row.last_trigger_to_end_s = dur - trig_t(end);

    % Hints of hand-placed triggers: times that are round to 0.01 s
    % (a real trigger falls on any sample, i.e. multiples of 1/512 s)
    is_round = abs(trig_t * 100 - round(trig_t * 100)) < 1e-6;
    row.n_trig_round_times   = sum(is_round);
    row.round_trigger_index  = strjoin(arrayfun(@num2str, find(is_round), 'UniformOutput', false), ' ');
    row.ends_1s_after_last   = abs(row.last_trigger_to_end_s - 1) < 1e-3;

    % Truncated start: not enough data before the first trigger for a full
    % backward epoch
    row.truncated_start         = (trig_t(1) + cfg.epoch_win_recommended(1)) < 0;
    row.first_epoch_available_s = min(trig_t(1), -cfg.epoch_win_recommended(1));
    if row.truncated_start
        flags{end+1} = sprintf(['TRUNCATED START: only %.1f s of data before the first ' ...
                                'trigger (first epoch shorter than %.1f s).'], ...
                               trig_t(1), -cfg.epoch_win_recommended(1));
    end
    % The last epoch needs at least +0.1 s of data after the last trigger
    if row.last_trigger_to_end_s < cfg.epoch_win_recommended(2)
        flags{end+1} = sprintf('Only %.2f s of data after the last trigger (epoch window needs %.1f s).', ...
                               row.last_trigger_to_end_s, cfg.epoch_win_recommended(2));
    end

    if nT > 1
        iti = diff(trig_t);
        row.iti_median_s         = median(iti);
        row.iti_min_s            = min(iti);
        row.iti_max_s            = max(iti);
        row.n_duplicate_triggers = sum(iti < 1);
        out_rng = iti < cfg.iti_expected_range_s(1) | iti > cfg.iti_expected_range_s(2);
        row.n_iti_out_of_range = sum(out_rng);

        if row.n_duplicate_triggers > 0
            flags{end+1} = sprintf('%d trigger pair(s) less than 1 s apart (duplicates?).', ...
                                   row.n_duplicate_triggers);
        end
        if row.n_iti_out_of_range > 0
            flags{end+1} = sprintf('%d inter-trigger interval(s) outside [%g %g] s: %s', ...
                row.n_iti_out_of_range, cfg.iti_expected_range_s, mat2str(round(iti(out_rng), 1)));
        end
    end
end

% ================= 4. RAW SIGNAL (CMS reference) =========================
% DC offset = median of each channel (BioSemi data are DC-coupled)
offset = median(X, 2);
X = X - offset;                                    % remove offset (in memory)

raw_ptp = max(X, [], 2) - min(X, [], 2);           % peak-to-peak, mostly drift

% Linear drift in microvolts per minute (least-squares slope)
t_min = (0:N-1) / srate / 60;
t_min = t_min - mean(t_min);
drift_slope = (X * t_min') / sum(t_min.^2);

% Longest run of identical consecutive samples (disconnection / saturation)
flat_run_s = zeros(nCh, 1);
for c = 1:nCh
    L = longest_run(diff(X(c, :)) == 0);
    if L > 0
        flat_run_s(c) = (L + 1) / srate;
    end
end

% Saturated channels (offset at the limit of the amplifier input range)
saturated = abs(offset) > cfg.saturation_uV;
row.n_saturated_channels = sum(saturated);
row.saturated_channels   = strjoin(labels(saturated), ' ');
if any(saturated)
    flags{end+1} = sprintf('SATURATED channel(s) (|offset| > %.0f mV): %s', ...
                           cfg.saturation_uV / 1000, row.saturated_channels);
end

row.eeg_offset_median_abs_uV = median(abs(offset(eeg_idx)));
row.eeg_offset_max_abs_uV    = max(abs(offset(eeg_idx(~saturated(eeg_idx)))));
if ~isempty(ext_idx)
    row.ext_offset_max_abs_uV = max(abs(offset(ext_idx)));
end
row.eeg_raw_ptp_median_uV           = median(raw_ptp(eeg_idx));
row.eeg_drift_median_abs_uV_per_min = median(abs(drift_slope(eeg_idx)));

% ================= 5. TEMPORARY 1 Hz HIGH-PASS (EEG channels) ============
EEGf      = pop_select(EEG, 'channel', eeg_idx);
EEGf.data = X(eeg_idx, :);                         % offset-removed raw data
EEGf      = pop_eegfiltnew(EEGf, cfg.qc_hp_hz, []); % high-pass (FIR, zero-phase)
Xf        = double(EEGf.data);

% SD in the original CMS reference: used to find flat/extreme channels
sd_cms = std(Xf, 0, 2);
flat   = sd_cms < cfg.flat_sd_uV | flat_run_s(eeg_idx) > cfg.max_flat_run_s | saturated(eeg_idx);
z_cms  = robust_z(log(sd_cms + eps));

% Channels allowed into the TEMPORARY average reference: not flat and not
% extreme. This keeps a few very bad channels from contaminating it.
ref_ok = ~flat & abs(z_cms) <= cfg.robust_z_thr;
if sum(ref_ok) < 0.5 * nE
    flags{end+1} = sprintf('Only %d EEG channels passed the reference check; using all non-flat channels.', sum(ref_ok));
    ref_ok = ~flat;
end
row.n_ref_channels = sum(ref_ok);

Xf = Xf - mean(Xf(ref_ok, :), 1);                  % temporary average reference

% ---- Slow-drift index on RAW (unfiltered) average-referenced data -------
raw_ref    = mean(X(eeg_idx(ref_ok), :), 1);       % raw average of good channels
nwin_long  = round(8 * srate);                     % 8-s windows (0.125 Hz res.)
[Praw, fraw] = pwelch((X(eeg_idx, :) - raw_ref)', hann(nwin_long), ...
                      round(nwin_long / 2), nwin_long, srate);
Praw       = Praw';                                % channels x frequencies
slow_ratio = bandpow(Praw, fraw, [0.25 1]) ./ bandpow(Praw, fraw, [1 40]);
clear Praw fraw

% ---- Power spectrum of each channel (1 Hz HP, average reference) --------
nwin  = round(cfg.welch_win_s * srate);
nover = round(nwin * cfg.welch_overlap);
[P, f] = pwelch(Xf', hann(nwin), nover, nwin, srate);
P = P';                                            % channels x frequencies

% Line noise: peak at 50 Hz relative to neighbouring frequencies (dB)
i_line  = abs(f - cfg.line_freq) <= 0.5;
i_neigh = (f >= cfg.line_freq - 5 & f <= cfg.line_freq - 2) | ...
          (f >= cfg.line_freq + 2 & f <= cfg.line_freq + 5);
line_db = 10 * log10(max(P(:, i_line), [], 2) ./ median(P(:, i_neigh), 2));

% High-frequency ratio (muscle / EMG proxy), excluding line harmonics
excl = false(size(f));
for h = cfg.line_freq * (1:floor(max(f) / cfg.line_freq))
    excl = excl | abs(f - h) <= 2;
end
hf_ratio = sum(P(:, f >= 30 & f <= 100 & ~excl), 2) ./ ...
           sum(P(:, f >= 1  & f <= 100 & ~excl), 2);

% ================= 6. TEMPORARY 40 Hz LOW-PASS -> time-domain metrics ====
EEGf.data = Xf;
clear Xf
EEGf = pop_eegfiltnew(EEGf, [], cfg.qc_lp_hz);     % low-pass (FIR, zero-phase)
Xb   = double(EEGf.data);                          % 1-40 Hz, average reference
clear EEGf

% Kurtosis (manual formula, no Statistics Toolbox needed) and covariance
C0   = Xb - mean(Xb, 2);
s2   = mean(C0.^2, 2);
kurt = mean(C0.^4, 2) ./ (s2.^2 + eps) - 3;        % excess kurtosis (Gaussian = 0)
Cv   = (C0 * C0') / (N - 1);                       % channel covariance matrix
clear C0
sd_bp = sqrt(diag(Cv));                            % SD of 1-40 Hz data
R     = Cv ./ (sd_bp * sd_bp' + eps);              % correlation matrix
ED    = diag(Cv) + diag(Cv)' - 2 * Cv;             % "electrical distance" = var(xi - xj)

% Nearest neighbours from 3-D coordinates
has_xyz = all(isfinite(xyz), 2);
D = pairwise_dist(xyz);
D(1:nE+1:end) = Inf;
[~, ord] = sort(D, 2);                             % NaN distances go last
k_nb = min(cfg.n_neighbours, nE - 1);
nb   = ord(:, 1:k_nb);

neigh_r = NaN(nE, 1);
min_ed  = NaN(nE, 1);
for c = 1:nE
    if ~has_xyz(c), continue; end
    nbc = nb(c, has_xyz(nb(c, :)));
    if isempty(nbc), continue; end
    neigh_r(c) = median(R(c, nbc));
    min_ed(c)  = min(ED(c, nbc));
end

% ---- Electrode bridges (Alschuler et al., 2014) -------------------------
% The ED of a pair does not depend on the reference (it is the variance of
% the difference between the two channels), so the average-referenced data
% can be used. Saturated / flat channels are excluded from the test.
br = detect_bridges_ed(Xb, srate, ~flat, cfg);
is_bridged = br.bridged_ch;
row.bridge_cutoff_scaled  = br.cutoff;
row.bridge_ed_scale_uV2   = br.median_ed_uV2;
row.n_bridged_pairs       = size(br.pairs, 1);
row.n_bridged_channels    = sum(is_bridged);
row.n_bridge_groups       = br.n_groups;
row.bridge_rank_loss      = br.rank_loss;
row.n_near_identical_pairs = br.n_rel_only_pairs;
row.bridge_note           = br.note;
row.bridged_channels      = strjoin(eeg_labels(is_bridged), ' ');
if row.n_bridged_pairs > 0
    flags{end+1} = sprintf(['%d bridged pair(s): %d channels in %d group(s); ' ...
                            'about %d data dimension(s) lost.'], row.n_bridged_pairs, ...
                           row.n_bridged_channels, row.n_bridge_groups, row.bridge_rank_loss);
end

% Percentage of 1-s windows exceeding the amplitude threshold
L_win = round(srate);
nW    = floor(N / L_win);
A     = reshape(abs(Xb(:, 1:nW * L_win)), nE, L_win, nW);
maxw  = reshape(max(A, [], 2), nE, nW);            % channels x windows
clear A
over   = maxw > cfg.amp_thr_uV;
pct_ch = 100 * mean(over, 2);

% ---- Correlation of each EEG channel with the bipolar EOG ---------------
% Bipolar derivations remove the CMS common-mode signal
heog = get_bipolar(X, role.HEOG_L, role.HEOG_R);
veog = get_bipolar(X, role.VEOG_U, role.VEOG_D);
ecg  = get_bipolar(X, role.CARD_L, role.CARD_R);

r_ocular = zeros(nE, 1);
for eog = {veog, heog}
    if isempty(eog{1}), continue; end
    e = fir_bandpass(eog{1}, srate, [cfg.qc_hp_hz cfg.qc_lp_hz]);
    e = e - mean(e);
    r = (Xb - mean(Xb, 2)) * e' ./ (sd_bp * std(e) * (N - 1) + eps);
    r_ocular = max(r_ocular, abs(r));
end
ocular = r_ocular > cfg.ocular_r_thr;

% Robust z-scores across channels (median / MAD)
z_sd   = robust_z(log(sd_bp + eps));
z_k    = robust_z(kurt);
z_line = robust_z(line_db);
z_hf   = robust_z(log(hf_ratio + eps));

% ---- Heuristic suspect-channel reasons (INFORMATIONAL) ------------------
% Reasons that ARE counted as "suspect":
thr = cfg.robust_z_thr;
reasons = repmat({''}, nE, 1);
reasons = add_reason(reasons, saturated(eeg_idx),                       'saturated');
reasons = add_reason(reasons, sd_cms < cfg.flat_sd_uV,                  'flat');
reasons = add_reason(reasons, flat_run_s(eeg_idx) > cfg.max_flat_run_s, 'flat-run');
reasons = add_reason(reasons, z_sd > thr & ~ocular,                     'high-SD');
reasons = add_reason(reasons, z_sd < -thr & ~is_bridged,                'low-SD');
reasons = add_reason(reasons, z_k > thr & ~ocular,                      'high-kurtosis');
reasons = add_reason(reasons, neigh_r < cfg.neigh_r_thr,                'low-neighbour-corr');
reasons = add_reason(reasons, z_line > thr & line_db > cfg.line_db_min, 'line-noise');
reasons = add_reason(reasons, z_hf > thr,                               'high-freq/EMG');
reasons = add_reason(reasons, abs(offset(eeg_idx)) > cfg.offset_thr_uV & ~saturated(eeg_idx), 'large-offset');
reasons = add_reason(reasons, pct_ch > cfg.pct_win_thr & ~ocular,       sprintf('often>%duV', cfg.amp_thr_uV));
suspect = ~cellfun(@isempty, reasons);

% Notes that are NOT counted as "suspect" (expected / reported separately)
notes = repmat({''}, nE, 1);
notes = add_reason(notes, ocular,                                  'ocular');
notes = add_reason(notes, ocular & (z_sd > thr | z_k > thr | pct_ch > cfg.pct_win_thr), ...
                   'ocular-artefact(expected)');
notes = add_reason(notes, is_bridged,                              'bridged');

row.n_suspect_channels = sum(suspect);
row.suspect_channels   = strjoin(eeg_labels(suspect), ' ');
row.n_ocular_channels  = sum(ocular);
row.ocular_channels    = strjoin(eeg_labels(ocular), ' ');

% ---- Subject-level summaries over the non-suspect channels --------------
good = ~suspect;
if sum(good) < 10, good = ~flat; end
row.eeg_sd_median_uV      = median(sd_bp(good));
row.eeg_line_db_median    = median(line_db(good));
row.eeg_hf_ratio_median   = median(hf_ratio(good));
row.eeg_slow_ratio_median = median(slow_ratio(good));
row.pct_win_any_good_gt_thr   = 100 * mean(any(over(good, :), 1));
row.pct_win_10pct_good_gt_thr = 100 * mean(mean(over(good, :), 1) > 0.10);

band_names = fieldnames(cfg.bands);
for b = 1:numel(band_names)
    bp = bandpow(P(good, :), f, cfg.bands.(band_names{b}));
    row.(['pow_' band_names{b} '_uV2']) = median(bp);
end
% Narrow-band peaks (e.g. the ~30 Hz peak of some subjects): peaks of the
% median spectrum that rise above the smoothed 1/f background by at least
% cfg.narrow_peak_db and are narrower than cfg.narrow_peak_width_hz
[row.narrow_peaks_hz, row.narrow_peak_max_db] = narrow_peaks(P(good, :), f, cfg);
if strlength(row.narrow_peaks_hz) > 0
    flags{end+1} = sprintf('Narrow-band spectral peak(s) at %s Hz (max +%.1f dB above background).', ...
                           row.narrow_peaks_hz, row.narrow_peak_max_db);
end
row.rel_alpha = median(bandpow(P(good, :), f, cfg.bands.alpha) ./ ...
                       bandpow(P(good, :), f, [1 45]));

if row.n_suspect_channels > 0.10 * nE
    flags{end+1} = sprintf('%d suspect EEG channels (> 10%% of the cap).', row.n_suspect_channels);
end

% Posterior vs anterior alpha over the whole recording: posterior alpha is
% normally larger. A ratio < 1 may indicate a coordinate orientation
% problem or poor data. Bridged channels are ALLOWED here.
post = choose_region(eeg_labels, xyz, suspect, cfg, 'posterior');
ant  = choose_region(eeg_labels, xyz, suspect, cfg, 'anterior');
row.posterior_channels = strjoin(eeg_labels(post), ' ');
if ~isempty(post) && ~isempty(ant)
    row.alpha_post_ant_ratio = median(bandpow(P(post, :), f, cfg.bands.alpha)) / ...
                               median(bandpow(P(ant, :),  f, cfg.bands.alpha));
    if row.alpha_post_ant_ratio < 1
        flags{end+1} = sprintf(['Alpha is not larger at posterior than anterior channels ' ...
            '(ratio %.2f): check data quality.'], row.alpha_post_ant_ratio);
    end
end
if numel(post) < 3
    flags{end+1} = 'Fewer than 3 usable posterior channels: alpha check unreliable.';
end

% ================= 7. EXTERNAL CHANNELS ==================================
blink_locs = [];
if ~isempty(heog)
    row.heog_sd_uV = std(fir_bandpass(heog, srate, [1 40]));
end
if ~isempty(veog)
    row.veog_sd_uV = std(fir_bandpass(veog, srate, [1 40]));
    blink_locs = detect_blinks(fir_bandpass(veog, srate, cfg.blink_band_hz), srate, cfg);
    row.blink_rate_total_pm = numel(blink_locs) / dur * 60;
else
    flags{end+1} = 'VEOG channels not available: blink check skipped.';
end
if ~isempty(ecg)
    [row.heart_rate_bpm, row.rr_cv] = estimate_hr(fir_bandpass(ecg, srate, cfg.ecg_band_hz), srate);
    if ~(row.heart_rate_bpm >= cfg.hr_plausible_bpm(1) && row.heart_rate_bpm <= cfg.hr_plausible_bpm(2)) ...
            || row.rr_cv > 0.25
        flags{end+1} = sprintf(['Heart-rate estimate unreliable (%.1f bpm, RR CV %.2f): ' ...
                                'check the CARD_L/CARD_R channels.'], row.heart_rate_bpm, row.rr_cv);
    end
end

% Mastoids, referenced to the same temporary average as the EEG
mL = []; mR = [];
if isfinite(role.MAST_L)
    mL = fir_bandpass(X(role.MAST_L, :) - raw_ref, srate, [cfg.qc_hp_hz cfg.qc_lp_hz]);
    row.mastL_sd_uV = std(mL);
end
if isfinite(role.MAST_R)
    mR = fir_bandpass(X(role.MAST_R, :) - raw_ref, srate, [cfg.qc_hp_hz cfg.qc_lp_hz]);
    row.mastR_sd_uV = std(mR);
end
if ~isempty(mL) && ~isempty(mR)
    cc = corrcoef(mL, mR);
    row.mast_LR_corr = cc(1, 2);
end
if row.mastL_sd_uV > 2 * row.eeg_sd_median_uV || row.mastR_sd_uV > 2 * row.eeg_sd_median_uV
    flags{end+1} = sprintf('Noisy mastoid(s): SD L = %.1f, R = %.1f uV (EEG median %.1f uV).', ...
                           row.mastL_sd_uV, row.mastR_sd_uV, row.eeg_sd_median_uV);
end
clear mL mR heog veog ecg

% ================= 8. TRIGGER SANITY CHECK (Berger effect + blinks) ======
isEO     = strcmp(trig_types(:), cfg.event_eo);
ep_t0    = NaN(nT, 1);
ep_t1    = NaN(nT, 1);
ep_avail = zeros(nT, 1);
ep_alpha = NaN(nT, 1);
ep_blink = NaN(nT, 1);
sumP = struct('EO', 0, 'EC', 0);
nP   = struct('EO', 0, 'EC', 0);
fk   = [];
has_veog = isfinite(role.VEOG_U) && isfinite(role.VEOG_D);

for k = 1:nT
    % Core window relative to the trigger, clipped to the recording
    t0 = max(0,   trig_t(k) + cfg.core_win(1));
    t1 = min(dur, trig_t(k) + cfg.core_win(2));
    ep_t0(k) = t0;
    ep_t1(k) = t1;
    ep_avail(k) = max(0, t1 - t0);
    if ep_avail(k) < cfg.min_core_s, continue; end

    s0 = floor(t0 * srate) + 1;
    s1 = floor(t1 * srate);

    % Posterior alpha power in this window
    if ~isempty(post)
        [Pk, fk] = pwelch(Xb(post, s0:s1)', hann(nwin), nover, nwin, srate);
        pm = mean(Pk, 2);                          % mean PSD across posterior channels
        ep_alpha(k) = log10(bandpow(pm', fk, cfg.bands.alpha));
        if isEO(k), cond = 'EO'; else, cond = 'EC'; end
        sumP.(cond) = sumP.(cond) + pm;
        nP.(cond)   = nP.(cond) + 1;
    end

    % Blink rate in this window (blinks per minute)
    if has_veog
        ep_blink(k) = sum(blink_locs >= s0 & blink_locs <= s1) / ep_avail(k) * 60;
    end
end
clear Xb

valid = isfinite(ep_alpha);
row.n_valid_epochs = sum(valid);
row.alpha_EO_log10 = mean(ep_alpha(valid &  isEO));
row.alpha_EC_log10 = mean(ep_alpha(valid & ~isEO));
row.berger_ratio   = 10^(row.alpha_EC_log10 - row.alpha_EO_log10);   % EC / EO

% Epoch-by-epoch agreement: are the high-alpha epochs the EC epochs?
vv = find(valid);
if numel(vv) >= 4 && any(isEO(vv)) && any(~isEO(vv))
    hi = ep_alpha(vv) > median(ep_alpha(vv));
    row.alpha_label_agreement = mean(hi == ~isEO(vv));
end

if ~isfinite(row.berger_ratio)
    row.alpha_verdict = "NOT COMPUTED";
elseif row.berger_ratio >= cfg.berger_ratio_ok && row.alpha_label_agreement >= 0.75
    row.alpha_verdict = "OK (EC > EO)";
elseif row.berger_ratio <= 1 / cfg.berger_ratio_ok && row.alpha_label_agreement <= 0.25
    row.alpha_verdict = "LIKELY INVERTED (EO > EC)";
else
    row.alpha_verdict = "WEAK/AMBIGUOUS";
end

% Individual alpha peak frequency (eyes closed, posterior channels)
if nP.EC > 0 && ~isempty(fk)
    pEC  = sumP.EC / nP.EC;
    isrc = find(fk >= cfg.iaf_search_hz(1) & fk <= cfg.iaf_search_hz(2));
    [~, im] = max(pEC(isrc));
    if im > 1 && im < numel(isrc)                  % a true peak, not an edge
        row.IAF_EC_hz = fk(isrc(im));
    end
end

% Blink rate EO vs EC
row.blink_rate_EO_pm = mean(ep_blink(isfinite(ep_blink) &  isEO));
row.blink_rate_EC_pm = mean(ep_blink(isfinite(ep_blink) & ~isEO));
if ~isfinite(row.blink_rate_EO_pm) || ~isfinite(row.blink_rate_EC_pm)
    row.blink_verdict = "NOT COMPUTED";
elseif row.blink_rate_EO_pm > cfg.blink_ratio_ok * row.blink_rate_EC_pm
    row.blink_verdict = "OK (EO > EC)";
elseif row.blink_rate_EC_pm > cfg.blink_ratio_ok * row.blink_rate_EO_pm
    row.blink_verdict = "LIKELY INVERTED (EC > EO)";
else
    row.blink_verdict = "WEAK/AMBIGUOUS";
end

% Overall verdict: inversion is suggested only when BOTH indicators agree;
% the alpha indicator is the primary one (blinks can be confused with eye
% movements or eyelid flutter with closed eyes).
a_ok  = row.alpha_verdict == "OK (EC > EO)";
a_inv = startsWith(row.alpha_verdict, "LIKELY INVERTED");
b_inv = startsWith(row.blink_verdict, "LIKELY INVERTED");
if a_inv && b_inv
    row.trigger_check_summary = "BOTH alpha and blinks suggest 11/12 labels are INVERTED";
elseif a_inv
    row.trigger_check_summary = "Alpha suggests inversion: inspect";
elseif a_ok && b_inv
    row.trigger_check_summary = "Alpha OK; blink pattern unusual (probably not inversion): inspect VEOG";
elseif a_ok
    row.trigger_check_summary = "Consistent";
else
    row.trigger_check_summary = "Weak/absent Berger effect: inspect";
end
if row.trigger_check_summary ~= "Consistent"
    flags{end+1} = sprintf('TRIGGER CHECK: %s (Berger ratio %.2f, agreement %.2f, blinks EO/EC %.1f/%.1f per min).', ...
        row.trigger_check_summary, row.berger_ratio, row.alpha_label_agreement, ...
        row.blink_rate_EO_pm, row.blink_rate_EC_pm);
end

% ================= 9. OUTPUT TABLES ======================================
% Bridge partners of each channel (text)
partners = repmat({''}, nE, 1);
for i = 1:size(br.pairs, 1)
    a = br.pairs(i, 1); b = br.pairs(i, 2);
    partners{a} = strtrim([partners{a} ' ' eeg_labels{b}]);
    partners{b} = strtrim([partners{b} ' ' eeg_labels{a}]);
end

subjE = repmat(string(row.subject), nE, 1);
chanT = table(subjE, eeg_idx(:), string(eeg_labels(:)), repmat("EEG", nE, 1), ...
    offset(eeg_idx), saturated(eeg_idx), raw_ptp(eeg_idx), drift_slope(eeg_idx), flat_run_s(eeg_idx), ...
    sd_cms, sd_bp, z_sd, kurt, z_k, neigh_r, min_ed, line_db, hf_ratio, slow_ratio, ...
    pct_ch, r_ocular, ocular, ref_ok, is_bridged, string(partners), ...
    suspect, string(reasons), string(notes), ...
    'VariableNames', chan_var_names());

nX = numel(ext_idx);
if nX > 0
    nanX = NaN(nX, 1);
    fX   = false(nX, 1);
    eX   = repmat("", nX, 1);
    extT = table(repmat(string(row.subject), nX, 1), ext_idx(:), string(labels(ext_idx))', ...
        repmat("EXT", nX, 1), offset(ext_idx), saturated(ext_idx), raw_ptp(ext_idx), ...
        drift_slope(ext_idx), flat_run_s(ext_idx), nanX, nanX, nanX, nanX, nanX, nanX, nanX, ...
        nanX, nanX, nanX, nanX, nanX, fX, fX, fX, eX, fX, eX, eX, ...
        'VariableNames', chan_var_names());
    chanT = [chanT; extT];
end

cond_lbl = repmat("EC", nT, 1);
cond_lbl(isEO) = "EO";
epT = table(repmat(string(row.subject), nT, 1), (1:nT)', string(trig_types(:)), cond_lbl, ...
    trig_t(:), ep_t0, ep_t1, ep_avail, ep_alpha, ep_blink, ...
    'VariableNames', {'subject','trigger_index','trigger_code','condition', ...
                      'trigger_time_s','core_start_s','core_end_s','core_available_s', ...
                      'posterior_alpha_log10_uV2','blinks_per_min'});

nB = size(br.pairs, 1);
if nB > 0
    brT = table(repmat(string(row.subject), nB, 1), ...
        string(eeg_labels(br.pairs(:, 1)))', string(eeg_labels(br.pairs(:, 2)))', ...
        br.pair_frac(:), br.pair_median_ed_uV2(:), pair_distance(xyz, br.pairs), ...
        'VariableNames', {'subject','chan_1','chan_2','frac_segments_below_cutoff', ...
                          'median_ED_uV2','distance_coord_units'});
else
    brT = table();
end

% ================= 10. FIGURE ============================================
if cfg.save_figures
    fig = [];
    try
        fig = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 1300 480]);
        % Force the light theme (MATLAB R2025a+ uses the desktop theme)
        try
            theme(fig, 'light');
        catch
        end

        % Left panel: spectra of all EEG channels
        ax1 = subplot(1, 2, 1); hold(ax1, 'on');
        fs  = f >= 1 & f <= 100;
        PdB = 10 * log10(P(:, fs) + eps);
        grey = ~suspect & ~is_bridged;
        if any(grey),       plot(ax1, f(fs), PdB(grey, :)',       'Color', [0.75 0.75 0.75]); end
        if any(is_bridged & ~suspect)
            plot(ax1, f(fs), PdB(is_bridged & ~suspect, :)', 'Color', [0.45 0.65 0.95]);
        end
        if any(suspect),    plot(ax1, f(fs), PdB(suspect, :)',    'Color', [0.85 0.2 0.2]); end
        plot(ax1, f(fs), median(PdB, 1), 'k', 'LineWidth', 2);
        xlabel(ax1, 'Frequency (Hz)');
        ylabel(ax1, 'Power (dB re 1 \muV^2/Hz)');
        title(ax1, sprintf('%s: grey ok, blue bridged, red suspect, black median', row.subject), ...
              'Interpreter', 'none');
        grid(ax1, 'on');

        % Right panel: posterior EO vs EC spectra
        ax2 = subplot(1, 2, 2); hold(ax2, 'on');
        if ~isempty(fk) && (nP.EO > 0 || nP.EC > 0)
            fs2 = fk >= 1 & fk <= 30;
            if nP.EO > 0
                plot(ax2, fk(fs2), 10*log10(sumP.EO(fs2) / nP.EO), 'b', 'LineWidth', 1.8, ...
                     'DisplayName', sprintf('EO (%s)', cfg.event_eo));
            end
            if nP.EC > 0
                plot(ax2, fk(fs2), 10*log10(sumP.EC(fs2) / nP.EC), 'r', 'LineWidth', 1.8, ...
                     'DisplayName', sprintf('EC (%s)', cfg.event_ec));
            end
            yl = ylim(ax2);
            hp = patch(ax2, [cfg.bands.alpha fliplr(cfg.bands.alpha)], [yl(1) yl(1) yl(2) yl(2)], ...
                       [1 0.9 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
            uistack(hp, 'bottom');
            if isfinite(row.IAF_EC_hz)
                xline(ax2, row.IAF_EC_hz, '--k', sprintf('IAF %.2f Hz', row.IAF_EC_hz), ...
                      'HandleVisibility', 'off');
            end
            legend(ax2, 'Location', 'northeast');
        end
        xlabel(ax2, 'Frequency (Hz)');
        ylabel(ax2, 'Power (dB re 1 \muV^2/Hz)');
        title(ax2, sprintf('Posterior (%d ch) | EC/EO = %.2f | %s', numel(post), ...
              row.berger_ratio, row.alpha_verdict), 'Interpreter', 'none');
        grid(ax2, 'on');

        exportgraphics(fig, fullfile(fig_folder, [char(row.subject) '_QC_v2.png']), ...
                       'Resolution', 110, 'BackgroundColor', 'white');
        close(fig);
    catch MEf
        flags{end+1} = ['Figure not saved: ' MEf.message];
        if ~isempty(fig) && isvalid(fig), close(fig); end
    end
end

end % qc_one_subject

% -------------------------------------------------------------------------
function br = detect_bridges_ed(Xb, srate, usable, cfg)
% DETECT_BRIDGES_ED  Electrode bridges from electrical-distance (ED)
% distributions, following Alschuler et al. (2014):
%   1. Cut the data in segments; for each segment and each channel pair
%      compute ED = variance of the difference between the two channels.
%   2. Scale all EDs by 100 / median(ED).
%   3. Histogram of the scaled EDs (bin 0.25), interpolated to 0.05.
%   4. Bridging is present only if there is a local peak near zero
%      (scaled ED <= 5) followed by a local minimum at scaled ED <= 10;
%      that minimum is the cut-off.
%   5. A pair is bridged if its ED is below the cut-off in at least 50 %
%      of the segments.
%   Safeguards added here (not in the original paper), because in a
%   128-channel cap the ED of normal neighbouring electrodes is also small:
%   - the valley after the near-zero peak must be clearly lower than the
%     peak (cfg.bridge_valley_max), i.e. two separate populations;
%   - the pair must also have a median relative ED (ED divided by the mean
%     variance of the two channels) below cfg.bridge_rel_ed_max.
%   Pairs meeting only the relative-ED criterion are counted separately
%   ("near-identical pairs") for information.
%   Unusable (flat / saturated) channels are excluded.

nE  = size(Xb, 1);
use = find(usable(:)');
br  = struct('pairs', zeros(0, 2), 'pair_frac', [], 'pair_median_ed_uV2', [], ...
             'bridged_ch', false(nE, 1), 'cutoff', NaN, 'median_ed_uV2', NaN, ...
             'n_groups', 0, 'rank_loss', 0, 'n_rel_only_pairs', 0, 'note', "");

L   = round(cfg.bridge_seg_s * srate);
nSg = floor(size(Xb, 2) / L);
nU  = numel(use);
if nSg < 10 || nU < 2, return; end

iu = find(triu(true(nU), 1));                      % index of each pair
EDs = zeros(numel(iu), nSg, 'single');
REL = zeros(numel(iu), nSg, 'single');
for k = 1:nSg
    seg = Xb(use, (k-1)*L+1 : k*L);
    seg = seg - mean(seg, 2);
    C   = (seg * seg') / (L - 1);                  % covariance in this segment
    v   = diag(C);
    E   = v + v' - 2 * C;                          % ED of all pairs
    Rel = E ./ ((v + v') / 2 + eps);               % ED relative to channel variance
    EDs(:, k) = E(iu);
    REL(:, k) = Rel(iu);
end
rel_ok = median(REL, 2) < cfg.bridge_rel_ed_max;   % near-identical pairs
clear REL

med_ed = median(EDs(:));
br.median_ed_uV2 = double(med_ed);
S = EDs * (100 / med_ed);                          % scaled ED

% Histogram of the scaled EDs between 0 and 100
edges   = 0:cfg.bridge_bin:100;
counts  = histcounts(S(:), edges);
centers = edges(1:end-1) + cfg.bridge_bin / 2;
fine    = 0:cfg.bridge_interp_step:100;
dens    = interp1(centers, counts, fine, 'linear', 0);

% Near-zero local peak followed by a local minimum
% Pairs that are near-identical by the relative criterion alone (info)
br.n_rel_only_pairs = sum(rel_ok);

is_max = islocalmax(dens) & fine > 0 & fine <= cfg.bridge_peak_max;
if ~any(is_max)
    br.note = "no near-zero peak in the ED distribution";
    return
end
cand  = find(is_max);
[~, j] = max(dens(cand));
i_pk  = cand(j);
is_min = islocalmin(dens, 'FlatSelection', 'first');  % first point of a flat valley
i_mn  = find(is_min & (1:numel(fine)) > i_pk, 1, 'first');
if isempty(i_mn) || fine(i_mn) > cfg.bridge_cutoff_max
    br.note = "near-zero peak not followed by a minimum <= cutoff_max";
    return
end
if dens(i_mn) > cfg.bridge_valley_max * dens(i_pk)
    br.note = "near-zero peak not separated from the main distribution";
    return
end
br.cutoff = fine(i_mn);
br.note   = "near-zero peak found";

% Pairs below the cut-off in at least 50 % of the segments
frac   = mean(S <= br.cutoff, 2);
is_br  = frac >= cfg.bridge_min_seg_frac & rel_ok;
[r, c] = ind2sub([nU nU], iu(is_br));
br.pairs              = [use(r(:))' use(c(:))'];
br.pair_frac          = double(frac(is_br));
br.pair_median_ed_uV2 = double(median(EDs(is_br, :), 2));

if isempty(br.pairs), return; end
br.bridged_ch(br.pairs(:)) = true;

% Groups of mutually bridged channels (connected components of a graph).
% A group of k bridged channels carries about 1 independent signal, so
% about (k - 1) data dimensions are lost per group.
G     = graph(br.pairs(:, 1), br.pairs(:, 2), ones(size(br.pairs, 1), 1), nE);
comp  = conncomp(G);
grp   = unique(comp(br.bridged_ch));
br.n_groups  = numel(grp);
br.rank_loss = sum(br.bridged_ch) - br.n_groups;
end

% -------------------------------------------------------------------------
function d = pair_distance(xyz, pairs)
% Euclidean distance between the two electrodes of each pair
if isempty(pairs)
    d = zeros(0, 1);
else
    d = sqrt(sum((xyz(pairs(:, 1), :) - xyz(pairs(:, 2), :)).^2, 2));
end
end

% -------------------------------------------------------------------------
function [txt, mx] = narrow_peaks(P, f, cfg)
% NARROW_PEAKS  Narrow peaks of the median spectrum above the 1/f background.
txt = ""; mx = NaN;
m  = 10 * log10(median(P, 1) + eps);
df = f(2) - f(1);
bg = movmedian(m, round(4 / df));                   % smooth 4-Hz background
r  = m - bg;
sel = f >= cfg.narrow_peak_range_hz(1) & f <= cfg.narrow_peak_range_hz(2) & abs(f - cfg.line_freq) > 2;
r(~sel) = 0;
[pk, loc, w] = findpeaks(r, 'MinPeakHeight', cfg.narrow_peak_db);
keep = w * df <= cfg.narrow_peak_width_hz;
if any(keep)
    txt = string(strjoin(arrayfun(@(x) sprintf('%.2f', x), f(loc(keep)), 'UniformOutput', false), ' '));
    mx  = max(pk(keep));
end
end

% -------------------------------------------------------------------------
function names = chan_var_names()
% Column names of the per-channel table (same for EEG and EXT rows)
names = {'subject','chan_index','label','type', ...
         'offset_uV','saturated','raw_ptp_uV','drift_uV_per_min','flat_run_s', ...
         'sd_hp_cms_uV','sd_bp_uV','z_sd','kurtosis','z_kurtosis', ...
         'neigh_r_median','min_neigh_ED_uV2','line_db','hf_ratio','slow_ratio', ...
         'pct_win_gt_thr','r_ocular','ocular','used_in_ref','bridged','bridge_partners', ...
         'suspect','reasons','notes'};
end

% -------------------------------------------------------------------------
function T = stack_tables(c)
% Concatenate non-empty tables of a cell array
c = c(~cellfun(@isempty, c));
if isempty(c)
    T = table();
else
    T = vertcat(c{:});
end
end

% -------------------------------------------------------------------------
function row = init_summary_row(cfg)
% INIT_SUMMARY_ROW  Empty summary row: strings = "", numbers = NaN.
s = {'subject','file','status','error_message','data_class','ref_field','ext_mapping', ...
     'other_events_detail','trigger_sequence','sequence_match','round_trigger_index', ...
     'saturated_channels','suspect_channels','ocular_channels','bridged_channels','bridge_note', ...
     'posterior_channels','alpha_verdict','blink_verdict','trigger_check_summary','flags_text', ...
     'narrow_peaks_hz'};
n = {'srate_hz','n_channels','n_eeg','n_ext','n_samples','duration_s','n_nan_inf', ...
     'n_eeg_valid_coords','n_events_total','n_boundary','n_trig_EO','n_trig_EC', ...
     'first_trigger_s','last_trigger_to_end_s','iti_median_s','iti_min_s','iti_max_s', ...
     'n_iti_out_of_range','n_duplicate_triggers','n_trig_out_of_bounds','n_trig_round_times', ...
     'first_epoch_available_s','n_saturated_channels','eeg_offset_median_abs_uV', ...
     'eeg_offset_max_abs_uV','ext_offset_max_abs_uV','eeg_raw_ptp_median_uV', ...
     'eeg_drift_median_abs_uV_per_min','eeg_slow_ratio_median','n_ref_channels', ...
     'eeg_sd_median_uV','eeg_line_db_median','eeg_hf_ratio_median','narrow_peak_max_db'};
band_names = fieldnames(cfg.bands);
for b = 1:numel(band_names)
    n{end+1} = ['pow_' band_names{b} '_uV2']; %#ok<AGROW>
end
n = [n, {'rel_alpha','pct_win_any_good_gt_thr','pct_win_10pct_good_gt_thr', ...
         'n_suspect_channels','n_ocular_channels', ...
         'bridge_cutoff_scaled','bridge_ed_scale_uV2','n_bridged_pairs','n_bridged_channels', ...
         'n_bridge_groups','bridge_rank_loss','n_near_identical_pairs','alpha_post_ant_ratio', ...
         'heog_sd_uV','veog_sd_uV','blink_rate_total_pm','heart_rate_bpm','rr_cv', ...
         'mastL_sd_uV','mastR_sd_uV','mast_LR_corr','n_valid_epochs', ...
         'alpha_EO_log10','alpha_EC_log10','berger_ratio','alpha_label_agreement', ...
         'IAF_EC_hz','blink_rate_EO_pm','blink_rate_EC_pm','n_flags'}];
l = {'has_ica','truncated_start','ends_1s_after_last'};

row = struct();
for k = 1:numel(s), row.(s{k}) = ""; end
for k = 1:numel(n), row.(n{k}) = NaN; end
for k = 1:numel(l), row.(l{k}) = false; end
end

% -------------------------------------------------------------------------
function out = enforce_fields(row, template)
% ENFORCE_FIELDS  Keep only template fields with consistent types, so all
% subjects can be stacked into one table.
out = template;
fn  = fieldnames(template);
for k = 1:numel(fn)
    if ~isfield(row, fn{k}), continue; end
    v = row.(fn{k});
    if isstring(template.(fn{k}))
        if iscell(v), v = strjoin(v, ' '); end
        v = string(v);
        if isempty(v), v = ""; end
        v = v(1);
    elseif islogical(template.(fn{k}))
        if isempty(v), v = false; end
        v = logical(v(1));
    else
        if isempty(v), v = NaN; end
        v = double(v(1));
    end
    out.(fn{k}) = v;
end
end

% -------------------------------------------------------------------------
function types = get_chan_types(chanlocs, nCh)
% Channel 'type' field as a cell array of char ('' if missing)
if isempty(chanlocs) || ~isfield(chanlocs, 'type')
    types = repmat({''}, 1, nCh);
    return
end
types = arrayfun(@(c) to_char(c.type), chanlocs, 'UniformOutput', false);
types = types(:)';
end

% -------------------------------------------------------------------------
function [eeg_idx, ext_idx, role, mapping] = identify_channels(labels, types, nCh, cfg)
% IDENTIFY_CHANNELS  Separate scalp EEG and external channels and find the
% index of each external role (HEOG_L, ..., CARD_R).
is_ext  = strcmpi(types, 'EXT');
mapping = "EXT identified by channel type";
if ~any(is_ext)
    is_ext = ismember(upper(labels), upper(cfg.ext_labels)) | ...
             ~cellfun(@isempty, regexp(upper(labels), '^EXG\d+$', 'once'));
    mapping = "EXT identified by label";
end
if ~any(is_ext) && nCh == cfg.n_eeg + cfg.n_ext
    is_ext = false(1, nCh);
    is_ext(cfg.n_eeg+1:end) = true;
    mapping = "EXT identified by position (last 8 channels)";
end
ext_idx = find(is_ext);
eeg_idx = find(~is_ext);

role = struct();
by_position = false;
for r = 1:numel(cfg.ext_labels)
    k = find(strcmpi(labels, cfg.ext_labels{r}), 1);
    if isempty(k)
        if numel(ext_idx) == numel(cfg.ext_labels)
            k = ext_idx(r);                 % documented order EXT1..EXT8
            by_position = true;
        else
            k = NaN;
        end
    end
    role.(cfg.ext_labels{r}) = k;
end
if by_position
    mapping = mapping + "; roles by order (EXG1..EXG8)";
else
    mapping = mapping + "; roles by label";
end
end

% -------------------------------------------------------------------------
function xyz = get_xyz(locs, idx)
% GET_XYZ  N x 3 matrix of coordinates (NaN where missing)
xyz = NaN(numel(idx), 3);
if isempty(locs) || ~isfield(locs, 'X'), return; end
for k = 1:numel(idx)
    c = locs(idx(k));
    xyz(k, :) = [num_or_nan(c.X), num_or_nan(c.Y), num_or_nan(c.Z)];
end
end

% -------------------------------------------------------------------------
function [types, lat] = get_events(EEG)
% GET_EVENTS  Event types as char and latencies (samples, 1-based)
if isempty(EEG.event)
    types = {};
    lat   = [];
    return
end
types = arrayfun(@(e) to_char(e.type), EEG.event, 'UniformOutput', false);
lat   = arrayfun(@(e) num_or_nan(e.latency), EEG.event);
types = types(:)';
lat   = lat(:)';
end

% -------------------------------------------------------------------------
function s = match_sequence(seq, C1, C2)
% MATCH_SEQUENCE  Compare the observed trigger sequence with C1 and C2.
if isempty(seq)
    s = "NO TRIGGERS";
elseif isequal(seq, C1)
    s = "C1";
elseif isequal(seq, C2)
    s = "C2";
else
    p1 = is_contig_sub(seq, C1);
    p2 = is_contig_sub(seq, C2);
    if p1 && p2
        s = "C1 or C2 (partial, ambiguous)";
    elseif p1
        s = "C1 (partial)";
    elseif p2
        s = "C2 (partial)";
    else
        s = "NO MATCH";
    end
end
end

% -------------------------------------------------------------------------
function tf = is_contig_sub(q, T)
% TRUE if q appears as a contiguous block inside T
tf = false;
n  = numel(q);
if n > numel(T), return; end
for s = 1:numel(T) - n + 1
    if isequal(T(s:s+n-1), q)
        tf = true;
        return
    end
end
end

% -------------------------------------------------------------------------
function L = longest_run(v)
% LONGEST_RUN  Length of the longest run of TRUE values in a logical vector
v = v(:)';
d = diff([false v false]);
starts = find(d == 1);
ends   = find(d == -1);
if isempty(starts)
    L = 0;
else
    L = max(ends - starts);
end
end

% -------------------------------------------------------------------------
function z = robust_z(v)
% ROBUST_Z  (v - median) / (1.4826 * MAD): z-score resistant to outliers
v   = v(:);
med = median(v, 'omitnan');
s   = 1.4826 * median(abs(v - med), 'omitnan');
z   = (v - med) / (s + eps);
end

% -------------------------------------------------------------------------
function s = robust_sd(x)
x = x(:);
s = 1.4826 * median(abs(x - median(x)));
end

% -------------------------------------------------------------------------
function q = qtile(v, p)
% QTILE  Simple quantile (avoids the Statistics Toolbox)
v = sort(v(isfinite(v)));
if isempty(v), q = NaN; return; end
idx = max(1, min(numel(v), round(p * (numel(v) - 1)) + 1));
q = v(idx);
end

% -------------------------------------------------------------------------
function bp = bandpow(P, f, band)
% BANDPOW  Integrated power in a band (rows of P = channels)
idx = f >= band(1) & f <= band(2);
bp  = trapz(f(idx), P(:, idx), 2);
end

% -------------------------------------------------------------------------
function D = pairwise_dist(xyz)
% PAIRWISE_DIST  Euclidean distance between all pairs of electrodes
sq = sum(xyz.^2, 2);
D  = sqrt(max(0, sq + sq' - 2 * (xyz * xyz')));
end

% -------------------------------------------------------------------------
function reasons = add_reason(reasons, mask, txt)
% ADD_REASON  Append a text label to the channels where mask is true
for i = find(mask(:))'
    if isempty(reasons{i})
        reasons{i} = txt;
    else
        reasons{i} = [reasons{i} ',' txt];
    end
end
end

% -------------------------------------------------------------------------
function idx = choose_region(eeg_labels, xyz, suspect, cfg, which)
% CHOOSE_REGION  Posterior or anterior channel group.
%   EEGLAB convention: +X points to the nose, +Z to the vertex.
%   The lowest 10 % of channels (neck / face) are excluded.
%   Only "suspect" channels are removed (bridged channels are kept).
if strcmp(which, 'posterior') && ~isempty(cfg.posterior_labels)
    idx = find(ismember(upper(eeg_labels), upper(cfg.posterior_labels)));
else
    cand = find(all(isfinite(xyz), 2));
    if isempty(cand), idx = []; return; end
    zthr = qtile(xyz(cand, 3), 0.10);
    if strcmp(which, 'posterior')
        xthr = qtile(xyz(cand, 1), cfg.posterior_frac);
        idx  = cand(xyz(cand, 1) <= xthr & xyz(cand, 3) > zthr);
    else
        xthr = qtile(xyz(cand, 1), 1 - cfg.posterior_frac);
        idx  = cand(xyz(cand, 1) >= xthr & xyz(cand, 3) > zthr);
    end
end
idx = idx(:)';
idx = idx(~suspect(idx));
end

% -------------------------------------------------------------------------
function y = get_bipolar(X, a, b)
% GET_BIPOLAR  Channel a minus channel b (empty if either is missing)
if isfinite(a) && isfinite(b)
    y = X(a, :) - X(b, :);
else
    y = [];
end
end

% -------------------------------------------------------------------------
function y = fir_bandpass(x, srate, band)
% FIR_BANDPASS  Zero-phase FIR band-pass (Hamming window) for one signal.
%   Filter order follows the usual rule: order = 3.3 * fs / transition.
x     = double(x(:)');
trans = min(band(1), 2);                      % transition bandwidth (Hz)
order = 2 * ceil(1.65 * srate / trans);       % even order
b     = fir1(order, band / (srate / 2), 'bandpass');
if numel(x) <= 3 * order
    y = x - mean(x);
    return
end
y = filtfilt(b, 1, x - mean(x));
end

% -------------------------------------------------------------------------
function locs = detect_blinks(x, srate, cfg)
% DETECT_BLINKS  Approximate blink detection on band-passed bipolar VEOG.
x = x(:)';
% Make blinks positive (polarity depends on the electrode montage)
if abs(qtile(x, 0.005)) > abs(qtile(x, 0.995))
    x = -x;
end
thr = max(cfg.blink_min_uV, median(x) + cfg.blink_mad_k * robust_sd(x));
[~, locs] = findpeaks(x, 'MinPeakHeight', thr, ...
                         'MinPeakDistance', round(cfg.blink_min_dist_s * srate));
end

% -------------------------------------------------------------------------
function [hr, rr_cv] = estimate_hr(x, srate)
% ESTIMATE_HR  Heart rate from band-passed bipolar ECG (R-peak detection)
hr = NaN; rr_cv = NaN;
x = x(:)';
if abs(qtile(x, 0.001)) > abs(qtile(x, 0.999))   % make R peaks positive
    x = -x;
end
thr = 0.4 * qtile(x, 0.999);
if ~(thr > 0), return; end
[~, locs] = findpeaks(x, 'MinPeakHeight', thr, 'MinPeakDistance', round(0.33 * srate));
rr = diff(locs) / srate;
rr = rr(rr > 0.33 & rr < 2);                       % plausible RR intervals only
if numel(rr) < 10, return; end
hr    = 60 / median(rr);
rr_cv = std(rr) / mean(rr);
end

% -------------------------------------------------------------------------
function s = to_char(v)
% TO_CHAR  Convert numbers / strings / chars to a trimmed char row
if isnumeric(v) || islogical(v)
    if isempty(v), s = ''; else, s = num2str(v); end
elseif isstring(v)
    s = char(v);
elseif ischar(v)
    s = strtrim(v);
else
    s = '';
end
end

% -------------------------------------------------------------------------
function v = num_or_nan(x)
if isnumeric(x) && ~isempty(x), v = double(x(1)); else, v = NaN; end
end

% -------------------------------------------------------------------------
function s = ref2str(r)
if ischar(r) || isstring(r)
    s = char(r);
elseif isnumeric(r)
    s = mat2str(r);
else
    s = 'unknown';
end
end

% -------------------------------------------------------------------------
function s = describe(x)
% DESCRIBE  "median [min, max] (n)" of a numeric vector
x = x(isfinite(x));
if isempty(x)
    s = 'n/a';
else
    s = sprintf('median %.2f  [min %.2f, max %.2f]  (n=%d)', median(x), min(x), max(x), numel(x));
end
end

% -------------------------------------------------------------------------
function write_report_header(fid, cfg, nS)
p = @(varargin) fprintf(fid, varargin{:});
p('%s\n', repmat('=', 1, 78));
p('  QUALITY CHECK v2 - RESTING-STATE EEG (BioSemi 128 + 8 EXT)\n');
p('  Generated : %s\n', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));
p('  Input     : %s\n', cfg.input_folder);
p('  Files     : %d\n', nS);
p('%s\n\n', repmat('=', 1, 78));
p('HOW TO READ THIS REPORT\n');
p('  - Input files are NOT modified. All filtering and re-referencing below is\n');
p('    done on temporary copies, only to compute descriptive numbers.\n');
p('  - Events %s are known leftovers: listed with their time, then IGNORED.\n', strjoin(cfg.ignore_events, ', '));
p('  - Offsets / drifts: raw data in the BioSemi CMS reference.\n');
p('  - Channel and spectral metrics: %g Hz high-pass (+ %g Hz low-pass for\n', cfg.qc_hp_hz, cfg.qc_lp_hz);
p('    time-domain metrics), temporary average reference computed WITHOUT the\n');
p('    flat / saturated / extreme channels.\n');
p('  - robust z = (value - median) / (1.4826 * MAD) across the EEG channels of\n');
p('    the same subject. |z| > %g is flagged.\n', cfg.robust_z_thr);
p('  - line_dB    : power at %g Hz relative to the surrounding frequencies.\n', cfg.line_freq);
p('  - hf_ratio   : power 30-100 Hz / power 1-100 Hz (muscle / EMG proxy).\n');
p('  - slow_ratio : power 0.25-1 Hz / power 1-40 Hz in raw data (drift / sweat).\n');
p('  - often>%duV : %% of 1-s windows where the channel exceeds %g uV (1-40 Hz).\n', cfg.amp_thr_uV, cfg.amp_thr_uV);
p('  - SATURATED  : |DC offset| > %.0f mV (amplifier range +/-262 mV).\n', cfg.saturation_uV / 1000);
p('  - OCULAR     : |corr| with VEOG or HEOG > %.2f. High SD / kurtosis in\n', cfg.ocular_r_thr);
p('                 these channels is expected (blinks) and NOT counted as bad.\n');
p('  - BRIDGES    : Alschuler et al. (2014) electrical-distance method, %g-s\n', cfg.bridge_seg_s);
p('                 segments, all channel pairs. Reported SEPARATELY, not as bad.\n');
p('                 "dimensions lost" = bridged channels - bridge groups\n');
p('                 (each group of mutually bridged channels carries ~1 signal).\n');
p('                 Safeguards: separated near-zero peak (valley <= %.0f%% of peak)\n', 100 * cfg.bridge_valley_max);
p('                 and median ED / channel variance < %.2f.\n', cfg.bridge_rel_ed_max);
p('  - Suspect channels are HEURISTIC flags, not final bad-channel decisions.\n');
p('  - TRIGGER CHECK: posterior alpha (%g-%g Hz) and blink rate in the window\n', cfg.bands.alpha);
p('    [%g %g] s relative to each trigger. Expected: alpha EC > EO (Berger\n', cfg.core_win);
p('    effect) and blinks EO > EC. "Agreement" = fraction of epochs in which\n');
p('    high alpha coincides with the EC label (1 = perfect, 0 = fully inverted).\n');
p('    The alpha indicator is the primary one.\n\n');
end

% -------------------------------------------------------------------------
function write_subject_block(fid, row, flags, chanT, epT, cfg)
% WRITE_SUBJECT_BLOCK  Human-readable block for one subject
p = @(varargin) fprintf(fid, varargin{:});
p('\n%s\n', repmat('=', 1, 78));
p('SUBJECT: %s    (file: %s)    status: %s\n', row.subject, row.file, row.status);
p('%s\n', repmat('=', 1, 78));

if row.status == "ERROR"
    p('  ERROR: %s\n', row.error_message);
    return
end

p('\n[1] FILE\n');
p('  Sampling rate               : %g Hz\n', row.srate_hz);
p('  Duration                    : %.1f s (%.2f min)\n', row.duration_s, row.duration_s/60);
p('  Channels (EEG / EXT)        : %d (%d / %d)\n', row.n_channels, row.n_eeg, row.n_ext);
p('  EXT identification          : %s\n', row.ext_mapping);
p('  EEG channels with coords    : %d\n', row.n_eeg_valid_coords);
p('  NaN/Inf samples             : %d\n', row.n_nan_inf);

p('\n[2] TRIGGERS\n');
p('  Events total / boundaries   : %d / %d\n', row.n_events_total, row.n_boundary);
p('  Ignored events (type@time)  : %s\n', row.other_events_detail);
p('  EO (%s) / EC (%s)           : %d / %d\n', cfg.event_eo, cfg.event_ec, row.n_trig_EO, row.n_trig_EC);
p('  Sequence                    : %s\n', row.trigger_sequence);
p('  Matches counterbalancing    : %s\n', row.sequence_match);
p('  Start -> first trigger      : %.3f s\n', row.first_trigger_s);
p('  Last trigger -> end         : %.3f s   (exactly 1.000 s: %s)\n', ...
  row.last_trigger_to_end_s, string(row.ends_1s_after_last));
p('  ITI median [min, max]       : %.2f [%.2f, %.2f] s\n', row.iti_median_s, row.iti_min_s, row.iti_max_s);
p('  Round trigger times (0.01 s): %d   (trigger #: %s)\n', row.n_trig_round_times, row.round_trigger_index);
p('  Truncated start             : %s (first epoch available: %.1f s)\n', ...
  string(row.truncated_start), row.first_epoch_available_s);

p('\n[3] RAW SIGNAL (CMS reference)\n');
p('  Saturated channels (%d)     : %s\n', row.n_saturated_channels, row.saturated_channels);
p('  EEG |offset| median / max   : %.0f / %.0f uV (max excludes saturated)\n', ...
  row.eeg_offset_median_abs_uV, row.eeg_offset_max_abs_uV);
p('  EXT |offset| max            : %.0f uV\n', row.ext_offset_max_abs_uV);
p('  EEG peak-to-peak (median)   : %.0f uV\n', row.eeg_raw_ptp_median_uV);
p('  EEG |linear drift| (median) : %.1f uV/min\n', row.eeg_drift_median_abs_uV_per_min);
p('  EEG slow_ratio (median)     : %.3f\n', row.eeg_slow_ratio_median);

p('\n[4] ELECTRODE BRIDGES (Alschuler et al., 2014)\n');
p('  Median ED of all pairs      : %.2f uV^2 (scale = 100)\n', row.bridge_ed_scale_uV2);
p('  Distribution                : %s\n', row.bridge_note);
p('  Cut-off (scaled ED)         : %.2f   (NaN = no bridges detected)\n', row.bridge_cutoff_scaled);
p('  Bridged pairs / channels    : %d / %d   in %d group(s)\n', ...
  row.n_bridged_pairs, row.n_bridged_channels, row.n_bridge_groups);
p('  Dimensions lost (approx.)   : %d\n', row.bridge_rank_loss);
p('  Near-identical pairs (rel.) : %d   (relative-ED criterion alone, info)\n', row.n_near_identical_pairs);
p('  Bridged channels            : %s\n', row.bridged_channels);

p('\n[5] EEG CHANNEL QUALITY (heuristic, informational)\n');
p('  Channels in temporary ref   : %d\n', row.n_ref_channels);
p('  SD 1-40 Hz (median)         : %.2f uV\n', row.eeg_sd_median_uV);
p('  line_dB (median)            : %.1f dB\n', row.eeg_line_db_median);
p('  hf_ratio (median)           : %.3f\n', row.eeg_hf_ratio_median);
p('  Narrow-band peaks           : %s Hz (max %.1f dB)\n', row.narrow_peaks_hz, row.narrow_peak_max_db);
p('  1-s windows > %g uV         : any good channel %.1f %% | >10%% of good channels %.1f %%\n', ...
  cfg.amp_thr_uV, row.pct_win_any_good_gt_thr, row.pct_win_10pct_good_gt_thr);
p('  Ocular channels (%d)        : %s\n', row.n_ocular_channels, row.ocular_channels);
p('  Suspect channels (%d)       : %s\n', row.n_suspect_channels, row.suspect_channels);
if ~isempty(chanT)
    S = chanT(chanT.type == "EEG" & chanT.suspect, :);
    if height(S) > 0
        p('    %-6s %8s %7s %8s %6s %7s %6s %8s  %s\n', 'label', 'SD(uV)', 'z_SD', 'kurt', ...
          'r_nb', 'lineDB', 'hf', '%>thr', 'reasons [notes]');
        for i = 1:height(S)
            p('    %-6s %8.2f %7.1f %8.1f %6.2f %7.1f %6.3f %8.1f  %s [%s]\n', S.label(i), ...
              S.sd_bp_uV(i), S.z_sd(i), S.kurtosis(i), S.neigh_r_median(i), S.line_db(i), ...
              S.hf_ratio(i), S.pct_win_gt_thr(i), S.reasons(i), S.notes(i));
        end
    end
end

p('\n[6] SPECTRAL OVERVIEW (median of non-suspect channels, 1 Hz HP, avg ref)\n');
band_names = fieldnames(cfg.bands);
for b = 1:numel(band_names)
    p('  %-10s %5.1f-%5.1f Hz     : %10.3f uV^2\n', band_names{b}, cfg.bands.(band_names{b}), ...
      row.(['pow_' band_names{b} '_uV2']));
end
p('  Relative alpha (8-13/1-45)  : %.3f\n', row.rel_alpha);
p('  Alpha posterior/anterior    : %.2f\n', row.alpha_post_ant_ratio);

p('\n[7] EXTERNAL CHANNELS\n');
p('  HEOG / VEOG SD (1-40 Hz)    : %.1f / %.1f uV\n', row.heog_sd_uV, row.veog_sd_uV);
p('  Blink rate (whole record)   : %.1f per min\n', row.blink_rate_total_pm);
p('  Heart rate (RR CV)          : %.1f bpm (%.2f)\n', row.heart_rate_bpm, row.rr_cv);
p('  Mastoid SD L / R (avg ref)  : %.1f / %.1f uV   (L-R corr %.2f)\n', ...
  row.mastL_sd_uV, row.mastR_sd_uV, row.mast_LR_corr);

p('\n[8] TRIGGER SANITY CHECK\n');
p('  Posterior channels          : %s\n', row.posterior_channels);
p('  log10 alpha EO / EC         : %.3f / %.3f\n', row.alpha_EO_log10, row.alpha_EC_log10);
p('  Berger ratio (EC/EO)        : %.2f\n', row.berger_ratio);
p('  Epoch label agreement       : %.2f\n', row.alpha_label_agreement);
p('  Alpha verdict               : %s\n', row.alpha_verdict);
p('  IAF (EC, posterior)         : %.2f Hz\n', row.IAF_EC_hz);
p('  Blinks/min EO / EC          : %.1f / %.1f   -> %s\n', ...
  row.blink_rate_EO_pm, row.blink_rate_EC_pm, row.blink_verdict);
p('  OVERALL                     : %s\n', row.trigger_check_summary);
if ~isempty(epT) && height(epT) > 0
    p('    %3s %5s %5s %10s %8s %18s %8s %10s %9s\n', '#', 'code', 'cond', 'trig(s)', ...
      'ITI(s)', 'core window (s)', 'avail', 'alpha_log', 'blink/min');
    for i = 1:height(epT)
        if i == 1, iti = NaN; else, iti = epT.trigger_time_s(i) - epT.trigger_time_s(i-1); end
        p('    %3d %5s %5s %10.3f %8.2f %8.1f-%-9.1f %8.1f %10.3f %9.1f\n', i, ...
          epT.trigger_code(i), epT.condition(i), epT.trigger_time_s(i), iti, ...
          epT.core_start_s(i), epT.core_end_s(i), epT.core_available_s(i), ...
          epT.posterior_alpha_log10_uV2(i), epT.blinks_per_min(i));
    end
end

p('\n[9] FLAGS (%d)\n', numel(flags));
if isempty(flags)
    p('  none\n');
else
    for i = 1:numel(flags)
        p('  - %s\n', flags{i});
    end
end
end

% -------------------------------------------------------------------------
function write_group_summary(fid, T, chanAll, cfg)
% WRITE_GROUP_SUMMARY  Overview across all subjects (end of the report)
p = @(varargin) fprintf(fid, varargin{:});
p('\n\n%s\n', repmat('#', 1, 78));
p('  GROUP SUMMARY\n');
p('%s\n', repmat('#', 1, 78));

ok  = T.status == "OK";
Tok = T(ok, :);
p('  Files: %d | processed OK: %d | errors: %d\n', height(T), sum(ok), sum(~ok));
if any(~ok)
    p('  Files with ERROR: %s\n', strjoin(T.subject(~ok)', ', '));
end
if height(Tok) == 0, return; end

p('\n  --- Recording ---\n');
p('  Sampling rates found       : %s\n', strjoin(string(unique(Tok.srate_hz))', ', '));
p('  Duration (s)               : %s\n', describe(Tok.duration_s));

p('\n  --- Triggers ---\n');
cats = unique(Tok.sequence_match);
for i = 1:numel(cats)
    sel = Tok.sequence_match == cats(i);
    p('  %-30s : %2d  (%s)\n', cats(i), sum(sel), strjoin(Tok.subject(sel)', ', '));
end
tr = Tok.truncated_start;
p('  Truncated start            : %d', sum(tr));
if any(tr)
    txt = arrayfun(@(i) sprintf('%s (%.1f s)', Tok.subject(i), Tok.first_epoch_available_s(i)), ...
                   find(tr)', 'UniformOutput', false);
    p('  -> %s', strjoin(txt, ', '));
end
p('\n');
p('  ITI median (s)             : %s\n', describe(Tok.iti_median_s));
e1 = Tok.ends_1s_after_last;
p('  Ends exactly 1 s after last trigger : %d  (not: %s)\n', sum(e1), strjoin(Tok.subject(~e1)', ', '));
p('  Subjects with round trigger times   : %d\n', sum(Tok.n_trig_round_times > 0));
p('  Ignored events, per subject:\n');
for i = 1:height(Tok)
    p('    %-6s %s\n', Tok.subject(i), Tok.other_events_detail(i));
end

p('\n  --- Trigger sanity check (Berger effect / blinks) ---\n');
p('  Berger ratio EC/EO         : %s\n', describe(Tok.berger_ratio));
p('  IAF (Hz)                   : %s\n', describe(Tok.IAF_EC_hz));
for v = ["alpha_verdict", "blink_verdict", "trigger_check_summary"]
    cats = unique(Tok.(v));
    p('  %s:\n', v);
    for i = 1:numel(cats)
        sel = Tok.(v) == cats(i);
        p('    %-70s : %2d  (%s)\n', cats(i), sum(sel), strjoin(Tok.subject(sel)', ', '));
    end
end

p('\n  --- Bridges ---\n');
p('  Subjects with bridges      : %d / %d\n', sum(Tok.n_bridged_pairs > 0), height(Tok));
p('  Bridged channels           : %s\n', describe(Tok.n_bridged_channels));
p('  Dimensions lost            : %s\n', describe(Tok.bridge_rank_loss));
[~, o] = sort(Tok.n_bridged_channels, 'descend');
o = o(1:min(8, numel(o)));
p('  Most bridged subjects      : %s\n', strjoin(arrayfun(@(i) sprintf('%s (%d ch)', ...
  Tok.subject(i), Tok.n_bridged_channels(i)), o', 'UniformOutput', false), ', '));

p('\n  --- Signal (across subjects) ---\n');
p('  Saturated channels         : %s\n', strjoin(arrayfun(@(i) sprintf('%s:%s', Tok.subject(i), ...
  Tok.saturated_channels(i)), find(Tok.n_saturated_channels > 0)', 'UniformOutput', false), ', '));
p('  EEG |offset| max (uV)      : %s\n', describe(Tok.eeg_offset_max_abs_uV));
p('  EEG drift (uV/min)         : %s\n', describe(Tok.eeg_drift_median_abs_uV_per_min));
p('  slow_ratio                 : %s\n', describe(Tok.eeg_slow_ratio_median));
p('  SD 1-40 Hz (uV)            : %s\n', describe(Tok.eeg_sd_median_uV));
p('  line_dB                    : %s\n', describe(Tok.eeg_line_db_median));
p('  hf_ratio                   : %s\n', describe(Tok.eeg_hf_ratio_median));
np_sel = strlength(Tok.narrow_peaks_hz) > 0;
p('  Narrow-band peaks          : %d subject(s)\n', sum(np_sel));
for i = find(np_sel)'
    p('    %-8s %s Hz (+%.1f dB)\n', Tok.subject(i), Tok.narrow_peaks_hz(i), Tok.narrow_peak_max_db(i));
end
p('  %% windows > %g uV (>10%%)  : %s\n', cfg.amp_thr_uV, describe(Tok.pct_win_10pct_good_gt_thr));
p('  Suspect EEG channels       : %s\n', describe(Tok.n_suspect_channels));
p('  Heart rate (bpm)           : %s\n', describe(Tok.heart_rate_bpm));
p('  Blinks/min (whole record)  : %s\n', describe(Tok.blink_rate_total_pm));

% Per-channel frequencies across subjects
if ~isempty(chanAll) && height(chanAll) > 0
    E = chanAll(chanAll.type == "EEG", :);
    [g, labs] = findgroups(E.label);
    nsub  = splitapply(@numel, E.suspect, g);
    nsusp = splitapply(@sum, double(E.suspect), g);
    nbr   = splitapply(@sum, double(E.bridged), g);

    frac = nsusp ./ nsub;
    [~, o] = sort(frac, 'descend');
    sel = o(frac(o) >= cfg.group_chan_frac);
    p('\n  --- EEG channels SUSPECT in >= %.0f%% of subjects ---\n', 100 * cfg.group_chan_frac);
    if isempty(sel), p('  none\n'); end
    for i = sel'
        p('  %-6s : %2d/%2d subjects  main reasons: %s\n', labs(i), nsusp(i), nsub(i), ...
          top_reasons(E.reasons(g == i & E.suspect)));
    end

    frac = nbr ./ nsub;
    [~, o] = sort(frac, 'descend');
    sel = o(frac(o) >= cfg.group_chan_frac);
    p('\n  --- EEG channels BRIDGED in >= %.0f%% of subjects ---\n', 100 * cfg.group_chan_frac);
    if isempty(sel), p('  none\n'); end
    for i = sel'
        p('  %-6s : %2d/%2d subjects\n', labs(i), nbr(i), nsub(i));
    end
end

p('\n  --- Subjects with flags ---\n');
for i = find(Tok.n_flags > 0)'
    p('  %-6s : %s\n', Tok.subject(i), Tok.flags_text(i));
end
p('\n(See each subject block above for details.)\n');
end

% -------------------------------------------------------------------------
function s = top_reasons(rs)
% TOP_REASONS  Three most frequent reasons in a list of comma-separated texts
tok = string(split(strjoin(cellstr(rs)', ','), ','));
tok = tok(tok ~= "");
if isempty(tok), s = ''; return; end
[u, ~, ic] = unique(tok);
cnt = accumarray(ic, 1);
[~, oo] = sort(cnt, 'descend');
s = strjoin(arrayfun(@(j) sprintf('%s x%d', u(j), cnt(j)), oo(1:min(3, numel(oo)))', ...
            'UniformOutput', false), ', ');
end
