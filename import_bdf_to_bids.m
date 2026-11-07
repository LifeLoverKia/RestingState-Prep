% =========================================================================
% import_bdf_to_bids.m
% =========================================================================
%
% PURPOSE
%   Build, from the ORIGINAL BioSemi .bdf files and the PsychoPy logs:
%
%   1. A BIDS raw dataset (EEG-BIDS):
%        BIDS/
%          dataset_description.json, README, participants.tsv/.json,
%          task-rest_events.json
%          sub-XX/eeg/
%            sub-XX_task-rest_eeg.bdf        (original file, unmodified copy)
%            sub-XX_task-rest_eeg.json       (recording metadata)
%            sub-XX_task-rest_channels.tsv   (channel types / status)
%            sub-XX_electrodes.tsv           (template coordinates)
%            sub-XX_coordsystem.json
%            sub-XX_task-rest_events.tsv     (CORRECTED condition onsets)
%          sourcedata/psychopy/
%            sub-XX_task-rest_psychopy.csv   (log WITHOUT the names column)
%
%   2. A BASE EEGLAB dataset for all future analyses (BIDS derivative):
%        BIDS/derivatives/eeglab-import/sub-XX/eeg/
%            sub-XX_task-rest_desc-import_eeg.set / .fdt
%      with:
%        - 128 scalp EEG + 8 external channels (renamed and typed),
%        - template coordinates,
%        - DC offset removed, resampled to 512 Hz,
%        - average reference of the 128 scalp channels (lossless, see below),
%        - corrected events (see EVENTS),
%        - NO filtering (besides the anti-aliasing filter of the
%          resampling), NO channel removal, NO ICA: those decisions belong
%          to the preprocessing step and depend on the analysis.
%
% EVENTS (the important part)
%   The hardware trigger (Status channel) has the same code for every
%   instruction, so the CONDITION of each period comes from PsychoPy.
%   For every subject:
%     a) the PsychoPy clock is aligned with the EEG clock using the
%        hardware triggers (offset that matches the most triggers; ties
%        resolved with the end of the recording);
%     b) for every instruction (= onset of a condition) the time is the
%        hardware pulse CLOSEST to the PsychoPy prediction (so, for double
%        pulses, the temporally coherent one is kept); if no pulse exists
%        (recording started late, or end of condition 8 not recorded) the
%        PsychoPy prediction is used;
%     c) events written to the .set:
%          'EO' / 'EC' : condition ONSET (instruction onset), with duration
%                        = recorded length of the condition
%          '11' / '12' : condition END (previous convention, kept for
%                        compatibility: 11 = end of EO, 12 = end of EC).
%                        If the true end was not recorded (condition 8),
%                        the marker is placed cfg.end_margin_s before the
%                        end of the file.
%      Every event carries its source ('hardware', 'psychopy',
%      'recording_start', 'recording_end') and its residual against the
%      PsychoPy prediction.
%
% REFERENCE
%   BioSemi data are recorded against CMS/DRL and must be re-referenced
%   offline. The base dataset is referenced to the average of the 128
%   scalp channels (all 136 channels are re-referenced to it). This is
%   LOSSLESS: any other reference (robust average without bad channels,
%   mastoids, REST, Laplacian/CSD) can still be computed exactly later.
%
% SUBJECT NAMES
%   S1 -> sub-01, S41 -> sub-41 (the original number is kept, so the
%   association with age, sex, etc. stays trivial; missing numbers are
%   simply absent).
%
% PRIVACY
%   The PsychoPy column "participant" (first names) is removed from the
%   copy saved in sourcedata. Nothing else from the logs is written to
%   the BIDS metadata.
%
% REQUIREMENTS
%   MATLAB R2026a, Signal Processing Toolbox, EEGLAB (eeg_emptyset,
%   pop_chanedit, pop_resample, pop_saveset, pop_loadset),
%   natural_sort_files.m. The .bdf files are read with a built-in reader
%   (BioSig plugin not needed).
%
% RUN TIME / DISK
%   About 1 min per subject. Disk: ~140 MB per .set/.fdt and one copy of
%   every .bdf (~8.5 GB in total for the BIDS raw folder).
% =========================================================================

clear; clc;

%% ========================================================================
%  0. CONFIGURATION
%  ========================================================================
cfg = struct();
base = 'C:\Users\chiar\Documents\MATLAB\FONDECYT\DATOS\Resting State';
cfg.bdf_folder  = fullfile(base, 'RAW');
cfg.psy_folder  = fullfile(base, 'EXCEL RS');
cfg.bids_root   = fullfile(base, 'BIDS');
cfg.deriv_name  = 'eeglab-import';
cfg.coord_file  = 'C:\Users\chiar\Documents\MATLAB\FONDECYT\code\128BIOSEMI_EXG.xyz';
% Optional: bridges found by qc_dataset_overview_v2.m (stored as metadata
% only, nothing is removed). Leave '' if not available.
cfg.bridges_csv = fullfile(base, 'QC_v2', 'QC_v2_bridges.csv');

cfg.overwrite   = false;   % false = skip subjects already converted
cfg.copy_bdf    = true;    % copy the original .bdf into the BIDS folder

% Connector swaps found during the QC (field = original subject code).
% The value gives, for the data recorded on inputs A, B, C, D, the cap
% bundle they actually belong to. 'ACBD' = inputs B and C were swapped
% (data recorded on B are the C electrodes and vice versa).
% S34: blinks appear on B16/B29/B28 instead of the frontopolar C16/C17/C29,
% and the spatial coherence of the data (median correlation with the 6
% nearest electrodes) rises from 0.70 (as labelled) to 0.86 with B<->C
% swapped; in the other 36 subjects the original assignment is the best.
cfg.connector_swap = struct('S34', 'ACBD');

% Signals
cfg.n_eeg        = 128;
cfg.n_ext        = 8;
cfg.target_srate = 512;    % Hz (Nyquist 256 Hz: enough for spectral,
                           % connectivity and gamma up to ~200 Hz)
cfg.ext_labels   = {'HEOG_L','HEOG_R','VEOG_U','VEOG_D','MAST_L','MAST_R','CARD_L','CARD_R'};
cfg.ext_types    = {'EOG','EOG','EOG','EOG','MAST','MAST','ECG','ECG'};          % EEGLAB
cfg.ext_bids     = {'EOG','EOG','EOG','EOG','MISC','MISC','ECG','ECG'};          % BIDS
% (type 'EOG' for the 4 EOG channels: the BIDS validator compares
%  EOGChannelCount with the number of channels of type 'EOG'; horizontal /
%  vertical is written in the description column)
cfg.ext_desc     = {'horizontal EOG, left outer canthus','horizontal EOG, right outer canthus', ...
                    'vertical EOG, above the eye','vertical EOG, below the eye', ...
                    'left mastoid','right mastoid','ECG left','ECG right'};
cfg.saturation_uV = 250000;   % |DC offset| above this = saturated (range +/-262 mV)

% Events
cfg.status_bits   = 65535;    % trigger bits of the BioSemi Status channel
cfg.match_tol_s   = 0.15;     % hardware pulse <-> PsychoPy prediction
cfg.end_margin_s  = 0.2;      % '11'/'12' placed this far from the file end if needed
cfg.n_conditions  = 8;
cfg.psychopy_version = '2022.1.3';   % from the PsychoPy logs (column psychopyVersion)
cfg.seq_C1 = {'EO','EC','EC','EO','EC','EO','EO','EC'};
cfg.seq_C2 = {'EC','EO','EO','EC','EO','EC','EC','EO'};

% Metadata for BIDS
cfg.dataset_name = 'FONDECYT resting-state EEG (eyes open / eyes closed)';
cfg.authors      = {'Chiara Saracini'};   % add co-authors here
cfg.license      = '';                    % e.g. 'CC0' or 'CC-BY-4.0' ('' = not written)
cfg.institution  = '';                    % e.g. university name ('' = not written)
cfg.power_line   = 50;

%% ========================================================================
%  1. SETUP
%  ========================================================================
if ~exist('pop_saveset', 'file')
    if exist('eeglab', 'file'), eeglab('nogui'); else, error('EEGLAB not found.'); end
end
if ~exist('resample', 'file'), error('Signal Processing Toolbox is required.'); end
if ~exist('natural_sort_files', 'file'), error('natural_sort_files.m not on the path.'); end
for f = {cfg.bdf_folder, cfg.psy_folder}
    if ~isfolder(f{1}), error('Folder not found:\n  %s', f{1}); end
end
if ~isfile(cfg.coord_file), error('Coordinate file not found:\n  %s', cfg.coord_file); end

deriv_root = fullfile(cfg.bids_root, 'derivatives', cfg.deriv_name);
make_dir(cfg.bids_root);
make_dir(fullfile(cfg.bids_root, 'sourcedata', 'psychopy'));
make_dir(deriv_root);

% Template coordinates (read once)
locs = readlocs(cfg.coord_file, 'filetype', 'xyz');
if numel(locs) < cfg.n_eeg + cfg.n_ext
    error('The coordinate file has %d rows (expected %d).', numel(locs), cfg.n_eeg + cfg.n_ext);
end

% Optional bridge table
bridges = table();
if ~isempty(cfg.bridges_csv) && isfile(cfg.bridges_csv)
    bridges = readtable(cfg.bridges_csv, 'TextType', 'string');
end

bdf_files = natural_sort_files(dir(fullfile(cfg.bdf_folder, '*.bdf')));
nS = numel(bdf_files);
fprintf('Found %d .bdf file(s).\n', nS);

log_fid = fopen(fullfile(deriv_root, 'import_log.txt'), 'w');
logp = @(varargin) fprintf(log_fid, varargin{:});
logp('import_bdf_to_bids.m  -  %s\n%s\n', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm')), repmat('=', 1, 78));

summary = {};
part    = {};

%% ========================================================================
%  2. LOOP OVER SUBJECTS
%  ========================================================================
for iS = 1:nS
    [~, orig_id] = fileparts(bdf_files(iS).name);
    num = str2double(regexp(orig_id, '\d+', 'match', 'once'));
    sub = sprintf('sub-%02d', num);
    fprintf('\n[%d/%d] %s -> %s\n', iS, nS, orig_id, sub);
    logp('\n%s\n%s  (%s)\n%s\n', repmat('-', 1, 78), sub, orig_id, repmat('-', 1, 78));

    eeg_dir   = fullfile(cfg.bids_root, sub, 'eeg');
    deriv_dir = fullfile(deriv_root, sub, 'eeg');
    stem      = [sub '_task-rest'];
    set_name  = [stem '_desc-import_eeg.set'];
    if ~cfg.overwrite && isfile(fullfile(deriv_dir, set_name))
        fprintf('  already converted: skipped (cfg.overwrite = false)\n');
        logp('  skipped (already converted)\n');
        continue
    end

    S = struct('subject', string(sub), 'original_id', string(orig_id), 'status', "OK", 'message', "");
    try
        % -------------------- 2a. Read the .bdf ------------------------
        bdf = read_bdf_full(fullfile(cfg.bdf_folder, bdf_files(iS).name), cfg);
        nKeep = cfg.n_eeg + cfg.n_ext;
        fs0   = bdf.srate;
        dur   = bdf.duration;
        S.original_srate = fs0;
        S.duration_s     = dur;

        % -------------------- 2b. PsychoPy log and events --------------
        psy = read_psychopy(fullfile(cfg.psy_folder, [orig_id '.csv']));
        if numel(psy.t) ~= cfg.n_conditions
            error('PsychoPy log has %d instructions (expected %d).', numel(psy.t), cfg.n_conditions);
        end
        ev = build_events(bdf.chg_t, bdf.grp_t, psy, dur, cfg);
        S.sequence          = string(ev.sequence);
        S.psy_offset_s      = ev.offset;
        S.n_hw_matched      = ev.n_matched;
        r_hw = abs(ev.resid_ms(ev.src_on == "hardware"));
        if isempty(r_hw), S.max_abs_resid_ms = NaN; else, S.max_abs_resid_ms = max(r_hw); end
        S.cond1_true_onset_s = ev.onset(1);
        S.cond8_true_end_s   = ev.offset_t(end);
        S.event_sources     = string(strjoin(cellstr(ev.src_on), ' '));

        % -------------------- 2c. Build the EEGLAB dataset -------------
        X = bdf.data;                                   % uV, CMS reference (136 channels)
        swap = '';
        if isfield(cfg.connector_swap, orig_id)
            swap = cfg.connector_swap.(orig_id);
            X = apply_connector_swap(X, swap);
            fprintf('  connector swap corrected: inputs ABCD -> cap bundles %s\n', swap);
        end
        bdf.swap = swap;
        bdf.data = [];                                  % free memory
        offset = median(X, 2);                          % DC offset of each channel
        saturated = abs(offset) > cfg.saturation_uV;
        X = X - offset;                                 % remove DC (before resampling)

        EEG = eeg_emptyset;
        EEG.setname  = [stem '_desc-import'];
        EEG.data     = single(X);
        clear X
        EEG.srate    = fs0;
        EEG.nbchan   = nKeep;
        EEG.pnts     = size(EEG.data, 2);
        EEG.trials   = 1;
        EEG.xmin     = 0;
        EEG.xmax     = (EEG.pnts - 1) / fs0;
        EEG.chanlocs = struct('labels', bdf.labels(1:nKeep));
        EEG = eeg_checkset(EEG);

        % Resample (EEGLAB pop_resample: anti-aliasing FIR + resampling)
        if fs0 ~= cfg.target_srate
            EEG = pop_resample(EEG, cfg.target_srate);
        end
        fs = EEG.srate;

        % Average reference of the 128 scalp channels, applied to all 136
        % channels (lossless). Then remove the remaining mean.
        EEG.data = EEG.data - mean(EEG.data(1:cfg.n_eeg, :), 1);
        EEG.data = EEG.data - mean(EEG.data, 2);
        EEG.ref  = 'average of the 128 scalp channels (all 136 channels re-referenced)';

        % Channel labels, types and template coordinates
        EEG = set_channels(EEG, locs, cfg);
        coord_warn = check_orientation(EEG);

        % Events
        EEG.event = make_eeglab_events(ev, fs, EEG.pnts, cfg);
        EEG = eeg_checkset(EEG, 'eventconsistency');

        % Documentation inside the dataset
        EEG.comments = sprintf(['Imported by import_bdf_to_bids.m on %s from %s (%g Hz). ' ...
            'DC removed, resampled to %g Hz, average reference (128 scalp). ' ...
            'Events from PsychoPy + hardware triggers.'], ...
            char(datetime('now', 'Format', 'yyyy-MM-dd')), bdf_files(iS).name, fs0, fs);
        EEG.etc.import.original_file  = bdf_files(iS).name;
        EEG.etc.import.original_id    = orig_id;
        EEG.etc.import.original_srate = fs0;
        EEG.etc.import.hardware_lowpass = bdf.prefilter;
        EEG.etc.import.dc_offset_uV   = offset;
        EEG.etc.import.saturated_channels = {EEG.chanlocs(saturated).labels};
        EEG.etc.import.psychopy_offset_s  = ev.offset;
        EEG.etc.import.event_table    = ev.table;
        EEG.etc.import.coordinate_check = coord_warn;
        EEG.etc.import.connector_swap = swap;
        if ~isempty(bridges) && isempty(swap)           % old QC labels are wrong if swapped
            b = bridges(bridges.subject == string(orig_id), :);
            EEG.etc.qc.bridged_pairs = [cellstr(b.chan_1) cellstr(b.chan_2)];
            EEG.etc.qc.bridges_source = cfg.bridges_csv;
        end

        % -------------------- 2d. Save the base .set -------------------
        make_dir(deriv_dir);
        pop_saveset(EEG, 'filename', set_name, 'filepath', deriv_dir, 'savemode', 'twofiles');

        % Verify that coordinates survived saving (known EEGLAB pitfall)
        chk = pop_loadset('filename', set_name, 'filepath', deriv_dir, 'loadmode', 'info');
        dxyz = max(abs([[chk.chanlocs(1:cfg.n_eeg).X] - [locs(1:cfg.n_eeg).X], ...
                        [chk.chanlocs(1:cfg.n_eeg).Y] - [locs(1:cfg.n_eeg).Y], ...
                        [chk.chanlocs(1:cfg.n_eeg).Z] - [locs(1:cfg.n_eeg).Z]]));
        S.coord_max_diff = dxyz;
        if dxyz > 1e-6
            coord_warn = [coord_warn ' | coordinates changed after saving!']; %#ok<AGROW>
        end
        S.coord_check = string(coord_warn);

        % -------------------- 2e. BIDS raw files -----------------------
        make_dir(eeg_dir);
        if cfg.copy_bdf
            copyfile(fullfile(cfg.bdf_folder, bdf_files(iS).name), fullfile(eeg_dir, [stem '_eeg.bdf']));
        end
        write_eeg_json(fullfile(eeg_dir, [stem '_eeg.json']), bdf, ev, cfg);
        % channels.tsv describes the .bdf file, whose labels follow the
        % amplifier inputs: map the saturation flags back to input order
        sat_in = saturated;
        for k = 1:4 * ~isempty(swap)
            sat_in((k - 1) * 32 + (1:32)) = saturated((double(swap(k)) - 65) * 32 + (1:32));
        end
        write_channels_tsv(fullfile(eeg_dir, [stem '_channels.tsv']), bdf, sat_in, cfg);
        write_electrodes(eeg_dir, stem, locs, bdf.labels, cfg);
        write_events_tsv(fullfile(eeg_dir, [stem '_events.tsv']), ev);
        copy_psychopy_anonymized(fullfile(cfg.psy_folder, [orig_id '.csv']), ...
            fullfile(cfg.bids_root, 'sourcedata', 'psychopy', [stem '_psychopy.csv']));

        % -------------------- 2f. Log ----------------------------------
        S.saturated = string(strjoin({EEG.chanlocs(saturated).labels}, ' '));
        logp('BDF %g Hz, %.1f s -> %g Hz | sequence %s | PsychoPy offset %.3f s | %d/9 matched\n', ...
             fs0, dur, fs, ev.sequence, ev.offset, ev.n_matched);
        logp('saturated: %s | coordinates: %s\n', S.saturated, coord_warn);
        logp('  %-4s %-4s %10s %10s %-16s %10s %10s %-14s %9s\n', 'cond', 'type', 'onset', 'end', ...
             'onset source', 'resid ms', 'rec. dur', 'end source', 'trunc.');
        for k = 1:cfg.n_conditions
            logp('  %-4d %-4s %10.3f %10.3f %-16s %10.1f %10.2f %-14s %9s\n', k, ev.cond{k}, ...
                 ev.onset(k), ev.offset_t(k), ev.src_on(k), ev.resid_ms(k), ev.rec_dur(k), ...
                 ev.src_off(k), string(ev.truncated(k)));
        end

        P = struct('participant_id', string(sub), 'original_id', string(orig_id), ...
                   'age', "n/a", 'sex', "n/a", 'handedness', "n/a", ...
                   'sequence', string(ev.sequence), 'original_srate_hz', string(fs0), ...
                   'cond1_recorded_s', string(round(ev.rec_dur(1), 2)), ...
                   'cond8_recorded_s', string(round(ev.rec_dur(end), 2)));
        part{end+1} = P; %#ok<SAGROW>

    catch ME
        S.status  = "ERROR";
        S.message = string(ME.message);
        if ~isempty(ME.stack)
            S.message = S.message + sprintf(' (%s, line %d)', ME.stack(1).name, ME.stack(1).line);
        end
        fprintf(2, '  ERROR: %s\n', S.message);
        logp('  ERROR: %s\n', S.message);
    end
    summary{end+1} = S; %#ok<SAGROW>
    clear EEG bdf chk
end

%% ========================================================================
%  3. DATASET-LEVEL FILES
%  ========================================================================
write_dataset_files(cfg, deriv_root, part);

% Summary table (fields can differ between OK and ERROR rows)
if ~isempty(summary)
    SumT = cell_structs_to_table(summary);
    writetable(SumT, fullfile(deriv_root, 'import_summary.csv'));
end
fclose(log_fid);

fprintf('\n%s\nDone. BIDS dataset in:\n  %s\nBase .set files in:\n  %s\n', ...
        repmat('=', 1, 70), cfg.bids_root, deriv_root);
fprintf(['Check the dataset with the BIDS validator ' ...
         '(https://bids-standard.github.io/bids-validator/).\n']);


%% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================

% -------------------------------------------------------------------------
function bdf = read_bdf_full(fn, cfg)
% READ_BDF_FULL  Read a BioSemi .bdf file: header, the first 136 channels
% (in microvolts, single precision) and the Status channel.
fid = fopen(fn, 'r', 'ieee-le');
if fid < 0, error('Cannot open %s', fn); end
cleanup = onCleanup(@() fclose(fid));

h = fread(fid, 256, '*char')';
bdf.start_date = strtrim(h(169:176));
header_bytes   = str2double(h(185:192));
n_rec          = str2double(h(237:244));
rec_dur        = str2double(h(245:252));
ns             = str2double(h(253:256));

rd = @(w) strtrim(cellstr(reshape(fread(fid, w * ns, '*char'), w, ns)'));
bdf.labels = rd(16)';
rd(80);                                       % transducer
bdf.units = rd(8)';
pmin = str2double(rd(8));
pmax = str2double(rd(8));
dmin = str2double(rd(8));
dmax = str2double(rd(8));
pref = rd(80);
nsamp = str2double(rd(8));
if any(nsamp ~= nsamp(1)), error('Channels with different sampling rates are not supported.'); end
spr = nsamp(1);

info = dir(fn);
n_rec_file = floor((info.bytes - header_bytes) / (ns * spr * 3));
if ~(n_rec > 0), n_rec = n_rec_file; end
n_rec = min(n_rec, n_rec_file);

bdf.ns        = ns;
bdf.srate     = spr / rec_dur;
bdf.duration  = n_rec * rec_dur;
bdf.prefilter = pref{1};
gain = (pmax - pmin) ./ (dmax - dmin);

nKeep = cfg.n_eeg + cfg.n_ext;
if ns < nKeep + 1, error('Only %d signals in the file.', ns); end
bdf.data = zeros(nKeep, spr * n_rec, 'single');
status   = zeros(spr * n_rec, 1);
fseek(fid, header_bytes, 'bof');
for r = 1:n_rec
    rec = reshape(fread(fid, ns * spr, 'bit24=>double'), spr, ns);
    idx = (r - 1) * spr + (1:spr);
    bdf.data(:, idx) = single(((rec(:, 1:nKeep) - dmin(1:nKeep)') .* gain(1:nKeep)' + pmin(1:nKeep)')');
    status(idx) = rec(:, ns);
end

% Status channel: every change to a non-baseline value is a candidate
st = double(bitand(uint32(mod(status, 2^24)), uint32(cfg.status_bits)));
baseline = mode(st);
chg = find(st(2:end) ~= st(1:end-1) & st(2:end) ~= baseline) + 1;
bdf.chg_t = (chg - 1) / bdf.srate;                  % all pulse onsets
keep = [true; diff(bdf.chg_t) > 0.5];
bdf.grp_t = bdf.chg_t(keep);                        % first pulse of each group
end

% -------------------------------------------------------------------------
function psy = read_psychopy(fn)
% READ_PSYCHOPY  Instruction onsets, conditions and audio files.
if ~isfile(fn), error('PsychoPy file not found: %s', fn); end
opts = detectImportOptions(fn, 'Delimiter', ',', 'VariableNamingRule', 'preserve');
opts = setvartype(opts, opts.VariableNames, 'char');
T = readtable(fn, opts);
instr = strtrim(T.('Instruction_trials'));
rows  = ~cellfun(@isempty, instr);
psy.t     = str2double(T.('sound_trial.started')(rows));
psy.instr = instr(rows)';
psy.cond  = cell(1, numel(psy.instr));
for i = 1:numel(psy.instr)
    if contains(psy.instr{i}, 'close', 'IgnoreCase', true)   % close / leave_closed
        psy.cond{i} = 'EC';
    else                                                      % open / leave_open
        psy.cond{i} = 'EO';
    end
end
end

% -------------------------------------------------------------------------
function ev = build_events(chg_t, grp_t, psy, dur, cfg)
% BUILD_EVENTS  Condition onsets / ends in EEG time.
%   Reference times (PsychoPy clock): 8 instruction onsets + end of
%   condition 8 (= last instruction + median condition length).
L   = median(diff(psy.t));
ref = [psy.t(:); psy.t(end) + L];
tol = cfg.match_tol_s;

% (a) coarse clock offset: most matches; ties -> condition 8 ends closest
%     to the end of the recording (the recording stops after condition 8)
cand = grp_t(:) - ref(:)';
cand = cand(:);
best_n = -1; best_gap = Inf; off = NaN;
for i = 1:numel(cand)
    n = 0;
    for j = 1:numel(ref)
        if any(abs(grp_t - (ref(j) + cand(i))) <= tol), n = n + 1; end
    end
    gap = abs(ref(end) + cand(i) - dur);
    if n > best_n || (n == best_n && gap < best_gap)
        best_n = n; best_gap = gap; off = cand(i);
    end
end
if best_n < 3
    error('Only %d hardware triggers match the PsychoPy log: alignment unreliable.', best_n);
end

% (b) refine twice: for every reference time take the pulse closest to the
%     prediction (handles double pulses), then update the offset
for it = 1:2
    t_ev = NaN(numel(ref), 1);
    for j = 1:numel(ref)
        [d, k] = min(abs(chg_t - (ref(j) + off)));
        if ~isempty(d) && d <= tol, t_ev(j) = chg_t(k); end
    end
    m = isfinite(t_ev);
    off = off + median(t_ev(m) - (ref(m) + off));
end
pred = ref + off;

% (c) final times and sources
t = pred;
src = repmat("psychopy", numel(ref), 1);
t(m) = t_ev(m);
src(m) = "hardware";
resid_ms = 1000 * (t - pred);

nC = cfg.n_conditions;
ev.offset    = off;
ev.n_matched = sum(m);
ev.cond      = psy.cond;
ev.instr     = psy.instr;
ev.onset     = t(1:nC);                 % true onset (can be < 0)
ev.offset_t  = t(2:nC+1);               % true end   (can be > duration)
ev.src_on    = src(1:nC);
ev.src_off   = src(2:nC+1);
ev.resid_ms  = resid_ms(1:nC);
ev.onset_rec = max(ev.onset, 0);
ev.end_rec   = min(ev.offset_t, dur);
ev.src_on(ev.onset < 0)          = "recording_start";
ev.src_off(ev.offset_t > dur)    = "recording_end";
ev.rec_dur   = ev.end_rec - ev.onset_rec;
ev.truncated = ev.onset < 0 | ev.offset_t > dur;
ev.dur_file  = dur;
if isequal(psy.cond, cfg.seq_C1)
    ev.sequence = 'C1';
elseif isequal(psy.cond, cfg.seq_C2)
    ev.sequence = 'C2';
else
    ev.sequence = 'other';
end
ev.table = table((1:nC)', string(ev.cond(:)), string(ev.instr(:)), ev.onset, ev.offset_t, ...
    ev.onset_rec, ev.end_rec, ev.rec_dur, ev.src_on, ev.src_off, ev.resid_ms, ev.truncated, ...
    'VariableNames', {'condition_index','condition','instruction','true_onset_s','true_end_s', ...
    'onset_s','end_s','recorded_duration_s','onset_source','end_source','residual_ms','truncated'});
end

% -------------------------------------------------------------------------
function E = make_eeglab_events(ev, fs, pnts, cfg)
% MAKE_EEGLAB_EVENTS  'EO'/'EC' at onset (with duration) and '11'/'12' at end
E = struct('type', {}, 'latency', {}, 'duration', {}, 'condition_index', {}, ...
           'source', {}, 'residual_ms', {}, 'truncated', {});
for k = 1:cfg.n_conditions
    % Onset event
    E(end+1).type = ev.cond{k}; %#ok<AGROW>
    E(end).latency  = ev.onset_rec(k) * fs + 1;
    E(end).duration = ev.rec_dur(k) * fs;
    E(end).condition_index = k;
    E(end).source   = char(ev.src_on(k));
    E(end).residual_ms = ev.resid_ms(k);
    E(end).truncated = double(ev.onset(k) < 0);

    % End marker (previous convention: 11 = end of EO, 12 = end of EC)
    t_end = ev.offset_t(k);
    src   = char(ev.src_off(k));
    if t_end > (pnts - 1) / fs - cfg.end_margin_s
        t_end = (pnts - 1) / fs - cfg.end_margin_s;
        src   = 'recording_end';
    end
    if strcmp(ev.cond{k}, 'EO'), code = '11'; else, code = '12'; end
    E(end+1).type = code; %#ok<AGROW>
    E(end).latency  = t_end * fs + 1;
    E(end).duration = 0;
    E(end).condition_index = k;
    E(end).source   = src;
    E(end).residual_ms = NaN;
    E(end).truncated = double(ev.offset_t(k) > ev.dur_file);
end
[~, o] = sort([E.latency]);
E = E(o);
end

% -------------------------------------------------------------------------
function EEG = set_channels(EEG, locs, cfg)
% SET_CHANNELS  Labels, types and template coordinates for the 136 channels
nE = cfg.n_eeg;
lab_bdf = {EEG.chanlocs.labels};
lab_xyz = {locs(1:nE).labels};
if ~isequal(upper(lab_bdf(1:nE)), upper(lab_xyz))
    warning('Labels of the .xyz file differ from the .bdf labels: the .bdf order is used.');
end
new = locs(1:nE + cfg.n_ext);
for c = 1:nE
    new(c).labels = lab_bdf{c};
    new(c).type   = 'EEG';
end
for c = 1:cfg.n_ext
    new(nE + c).labels = cfg.ext_labels{c};
    new(nE + c).type   = cfg.ext_types{c};
end
EEG.chanlocs = new;
EEG.chaninfo.nosedir = '+X';
EEG = pop_chanedit(EEG, 'convert', {'cart2all'});   % fill theta/radius/sph from X/Y/Z
end

% -------------------------------------------------------------------------
function msg = check_orientation(EEG)
% CHECK_ORIENTATION  Plausibility of the coordinates (EEGLAB: +X = nose,
% +Z = vertex). BioSemi 128: A1 ~ Cz (top), A19 ~ Pz (back), C21 ~ Fz (front).
L = {EEG.chanlocs.labels};
g = @(lab) EEG.chanlocs(strcmp(L, lab));
msg = 'ok';
try
    a1 = g('A1'); a19 = g('A19'); c21 = g('C21');
    Z = sort([EEG.chanlocs(1:128).Z]);
    if a1.Z < Z(round(0.9 * numel(Z))), msg = 'A1 is not near the vertex'; end
    if ~(a19.X < 0 && c21.X > 0), msg = 'front/back orientation unexpected (A19 should be posterior, C21 anterior)'; end
catch
    msg = 'check not possible';
end
end

% -------------------------------------------------------------------------
function Y = apply_connector_swap(X, swap)
% APPLY_CONNECTOR_SWAP  Move the four 32-channel blocks of the scalp data to
% their true cap positions. swap(k) = cap bundle of the data recorded on
% input k (k = A, B, C, D). External channels (129-136) are not changed.
Y = X;
for k = 1:4
    src = (k - 1) * 32 + (1:32);
    dst = (double(swap(k)) - double('A')) * 32 + (1:32);
    Y(dst, :) = X(src, :);
end
end

% -------------------------------------------------------------------------
function write_eeg_json(fn, bdf, ev, cfg)
J = struct();
J.TaskName = 'rest';
J.TaskDescription = ['Resting state: 8 periods of about 63 s alternating eyes open (EO) and ' ...
    'eyes closed (EC), in one of two counterbalanced sequences (C1: EO-EC-EC-EO-EC-EO-EO-EC; ' ...
    'C2: EC-EO-EO-EC-EO-EC-EC-EO). Each period starts with an audio instruction.'];
J.Instructions = 'Audio instructions to open / close the eyes or keep them open / closed.';
J.Manufacturer = 'BioSemi';
J.ManufacturersModelName = 'ActiveTwo';
J.CapManufacturer = 'BioSemi';
J.CapManufacturersModelName = 'headcap 128 channels, ABCD layout';
J.SamplingFrequency = bdf.srate;
J.EEGChannelCount = cfg.n_eeg;
J.EOGChannelCount = 4;
J.ECGChannelCount = 2;
J.EMGChannelCount = 0;
J.MISCChannelCount = 2 + max(0, bdf.ns - 1 - cfg.n_eeg - cfg.n_ext);
J.TriggerChannelCount = 1;
J.EEGReference = 'CMS/DRL (BioSemi active electrodes); re-referencing required offline';
J.EEGGround = 'CMS/DRL';
J.PowerLineFrequency = cfg.power_line;
J.SoftwareFilters = 'n/a';
J.HardwareFilters = hardware_filters(bdf.prefilter);
J.EEGPlacementScheme = 'BioSemi 128-channel ABCD radial layout';
if ~isempty(cfg.institution), J.InstitutionName = cfg.institution; end
J.RecordingType = 'continuous';
J.RecordingDuration = bdf.duration;
J.EpochLength = 0;
if bdf.srate ~= 2048
    J.SamplingFrequencyNote = sprintf(['Recorded at 2048 Hz and decimated to %g Hz with the ' ...
        'BioSemi ActiveTools decimate tool before archiving.'], bdf.srate);
end
if ~isempty(bdf.swap)
    J.ConnectorSwapNote = sprintf(['Electrode connectors were swapped during the recording: the ' ...
        'data of inputs A, B, C, D belong to cap bundles %s. The .bdf file is unmodified (its ' ...
        'labels follow the inputs); the swap is corrected in derivatives/eeglab-import.'], bdf.swap);
end
if any(ev.truncated)
    J.RecordingNote = 'Condition 1 and/or 8 only partially recorded: see events.tsv (column truncated).';
end
write_json(fn, J);
end

% -------------------------------------------------------------------------
function H = hardware_filters(prefilter)
% HARDWARE_FILTERS  BIDS object for the BioSemi hardware filter, built from
% the header field (e.g. 'HP: DC; LP: 208 Hz'). BIDS requires an object of
% objects: {"filter name": {"parameter": value}}.
lp = str2double(regexp(prefilter, 'LP:\s*([\d.]+)', 'tokens', 'once'));
F = struct();
F.Description = 'BioSemi ActiveTwo 5th-order sinc anti-aliasing (decimation) filter; DC coupled (no high-pass)';
if isfinite(lp), F.LowpassCutoffHz = lp; end
F.HeaderPrefilterField = prefilter;
H = struct('BioSemi_decimation_filter', F);
end

% -------------------------------------------------------------------------
function write_channels_tsv(fn, bdf, saturated, cfg)
% One row per signal of the .bdf file (BIDS requires all channels)
nE = cfg.n_eeg; nX = cfg.n_ext;
rows = cell(bdf.ns, 6);
for c = 1:bdf.ns
    name = bdf.labels{c};
    units = bdf.units{c};
    status = 'good'; sdesc = 'n/a'; desc = 'n/a';
    if c <= nE
        type = 'EEG';
        if saturated(c)
            status = 'bad'; sdesc = 'saturated (DC offset at the amplifier limit)';
        end
    elseif c <= nE + nX
        type = cfg.ext_bids{c - nE};
        desc = sprintf('%s (%s)', cfg.ext_labels{c - nE}, cfg.ext_desc{c - nE});
        if saturated(c)
            status = 'bad'; sdesc = 'saturated (DC offset at the amplifier limit)';
        end
    elseif strcmpi(name, 'Status')
        type = 'TRIG'; desc = 'BioSemi Status channel (trigger bits 0-15)';
    else
        type = 'MISC'; status = 'n/a'; desc = 'BioSemi auxiliary sensor input, not used in this study';
    end
    rows(c, :) = {name, type, units, desc, status, sdesc};
end
write_tsv(fn, {'name','type','units','description','status','status_description'}, rows);
end

% -------------------------------------------------------------------------
function write_electrodes(eeg_dir, stem, locs, labels, cfg)
nK = cfg.n_eeg + cfg.n_ext;
names = [labels(1:cfg.n_eeg), cfg.ext_labels];
rows = cell(nK, 4);
for c = 1:nK
    rows(c, :) = {names{c}, locs(c).X, locs(c).Y, locs(c).Z};
end
sub_stem = extractBefore(stem, '_task');           % electrodes take no task entity
write_tsv(fullfile(eeg_dir, [sub_stem '_electrodes.tsv']), {'name','x','y','z'}, rows);

r = median(sqrt([locs(1:cfg.n_eeg).X].^2 + [locs(1:cfg.n_eeg).Y].^2 + [locs(1:cfg.n_eeg).Z].^2));
if r > 50 && r < 150, u = 'mm'; elseif r > 5 && r < 15, u = 'cm'; elseif r > 0.05 && r < 0.15, u = 'm'; else, u = 'n/a'; end
J = struct();
J.EEGCoordinateSystem = 'Other';
J.EEGCoordinateUnits  = u;
J.EEGCoordinateSystemDescription = ['Template positions of the BioSemi 128-channel ABCD cap ' ...
    '(file 128BIOSEMI_EXG.xyz), not digitized. EEGLAB orientation: +X towards the nose, ' ...
    '+Y towards the left ear, +Z towards the vertex.'];
write_json(fullfile(eeg_dir, [sub_stem '_coordsystem.json']), J);
end

% -------------------------------------------------------------------------
function write_events_tsv(fn, ev)
n = numel(ev.onset);
rows = cell(n, 12);
for k = 1:n
    if strcmp(ev.cond{k}, 'EO'), tt = 'eyes_open'; v = 11; else, tt = 'eyes_closed'; v = 12; end
    rows(k, :) = {ev.onset_rec(k), ev.rec_dur(k), tt, v, k, ev.instr{k}, ...
                  ev.onset(k), ev.offset_t(k), char(ev.src_on(k)), char(ev.src_off(k)), ...
                  ev.resid_ms(k), lower(string(ev.truncated(k)))};
end
write_tsv(fn, {'onset','duration','trial_type','value','condition_index','instruction', ...
               'true_onset','true_end','onset_source','end_source','residual_ms','truncated'}, rows);
end

% -------------------------------------------------------------------------
function copy_psychopy_anonymized(src, dst)
% Copy the PsychoPy log without the 'participant' column (first names)
opts = detectImportOptions(src, 'Delimiter', ',', 'VariableNamingRule', 'preserve');
opts = setvartype(opts, opts.VariableNames, 'char');
T = readtable(src, opts);
drop = strcmpi(T.Properties.VariableNames, 'participant') | ...
       startsWith(T.Properties.VariableNames, 'Var');      % empty trailing column
T(:, drop) = [];
writetable(T, dst, 'Delimiter', ',');
end

% -------------------------------------------------------------------------
function write_dataset_files(cfg, deriv_root, part)
% dataset_description.json (raw and derivative), README, participants
D = struct();
D.Name = cfg.dataset_name;
D.BIDSVersion = '1.9.0';
D.DatasetType = 'raw';
D.Authors = cfg.authors;
if ~isempty(cfg.license), D.License = cfg.license; end
D.GeneratedBy = {struct('Name', 'import_bdf_to_bids.m', 'Description', ...
    'Conversion of the original BioSemi .bdf files and PsychoPy logs to BIDS')};
write_json(fullfile(cfg.bids_root, 'dataset_description.json'), D);

DD = struct();
DD.Name = [cfg.dataset_name ' - EEGLAB base import'];
DD.BIDSVersion = '1.9.0';
DD.DatasetType = 'derivative';
DD.GeneratedBy = {struct('Name', 'import_bdf_to_bids.m', 'Description', ...
    ['DC removal, resampling to ' num2str(cfg.target_srate) ' Hz, average reference of 128 scalp ' ...
     'channels, template coordinates, events from PsychoPy and hardware triggers. No filtering, ' ...
     'no channel removal, no ICA.'])};
write_json(fullfile(deriv_root, 'dataset_description.json'), DD);

readme = fullfile(cfg.bids_root, 'README');
if ~isfile(readme)
    f = fopen(readme, 'w');
    fprintf(f, ['%s\n\nResting-state EEG, BioSemi ActiveTwo, 128 scalp channels + 8 external.\n' ...
        'Eight periods of ~63 s (EO/EC), two counterbalanced sequences (see participants.tsv).\n\n' ...
        'Events: onsets of each period = onset of the audio instruction. The hardware trigger has\n' ...
        'the same code for all instructions; the condition comes from the PsychoPy log\n' ...
        '(sourcedata/psychopy, names removed). Times are the hardware pulse closest to the\n' ...
        'PsychoPy prediction, or the PsychoPy prediction when no pulse was recorded.\n' ...
        'Original subject codes S<n> are mapped to sub-<nn>.\n'], cfg.dataset_name);
    fclose(f);
end

T = fullfile(cfg.bids_root, 'task-rest_events.json');
J = struct();
% (onset and duration are standard BIDS columns: not redefined here.
%  onset = recorded onset, 0 if the period started before the recording;
%  duration = recorded duration of the period)
J.StimulusPresentation = struct('SoftwareName', 'PsychoPy', 'SoftwareVersion', cfg.psychopy_version);
J.trial_type = struct('Description', 'Condition', 'Levels', struct('eyes_open', 'eyes open', 'eyes_closed', 'eyes closed'));
J.value = struct('Description', 'Code used in previous analyses (11 = EO, 12 = EC)');
J.condition_index = struct('Description', 'Position of the period in the sequence (1-8)');
J.instruction = struct('Description', 'Audio instruction played at the onset (PsychoPy)');
J.true_onset = struct('Description', 'Onset in EEG time, negative if before the recording start (s)');
J.true_end = struct('Description', 'End in EEG time (next instruction / predicted end), may exceed the recording (s)');
J.onset_source = struct('Description', 'hardware = Status pulse closest to the PsychoPy prediction; psychopy = predicted from the log; recording_start = period started before the recording');
J.end_source = struct('Description', 'Same as onset_source, for the end of the period; recording_end = end not recorded');
J.residual_ms = struct('Description', 'Hardware time minus PsychoPy prediction (ms)');
J.truncated = struct('Description', 'true if the period is only partially recorded');
write_json(T, J);

% participants.tsv. If the file already exists (e.g. age / sex already
% filled in by hand), only the columns written by this script are updated
% and every other column is preserved.
if isempty(part), return; end
P = cell_structs_to_table(part);                 % all columns are strings
pf = fullfile(cfg.bids_root, 'participants.tsv');
user_cols = {'age', 'sex', 'handedness'};        % never overwritten
if isfile(pf)
    o = detectImportOptions(pf, 'FileType', 'text', 'Delimiter', '\t');
    o = setvartype(o, o.VariableNames, 'string');
    old = readtable(pf, o);
    for v = P.Properties.VariableNames             % add missing columns
        if ~ismember(v{1}, old.Properties.VariableNames)
            old.(v{1}) = repmat("n/a", height(old), 1);
        end
    end
    for i = 1:height(P)
        j = find(old.participant_id == P.participant_id(i), 1);
        if isempty(j)                               % new participant: append
            newrow = old(1, :);
            for v = old.Properties.VariableNames, newrow.(v{1}) = "n/a"; end
            for v = P.Properties.VariableNames, newrow.(v{1}) = P.(v{1})(i); end
            old = [old; newrow]; %#ok<AGROW>
        else                                        % existing: update script columns
            for v = setdiff(P.Properties.VariableNames, user_cols)
                old.(v{1})(j) = P.(v{1})(i);
            end
        end
    end
    P = old;
end
P = sortrows(P, 'participant_id');
rows = table2cell(P);
write_tsv(pf, P.Properties.VariableNames, rows);
J = struct();
J.original_id = struct('Description', 'Subject code used during the study');
J.age = struct('Description', 'Age', 'Units', 'years');
J.sex = struct('Description', 'Sex', 'Levels', struct('M', 'male', 'F', 'female'));
J.handedness = struct('Description', 'Handedness', 'Levels', struct('R', 'right', 'L', 'left', 'A', 'ambidextrous'));
J.sequence = struct('Description', 'Counterbalanced sequence', 'Levels', ...
    struct('C1', 'EO-EC-EC-EO-EC-EO-EO-EC', 'C2', 'EC-EO-EO-EC-EO-EC-EC-EO'));
J.original_srate_hz = struct('Description', 'Sampling rate of the archived .bdf file', 'Units', 'Hz');
J.cond1_recorded_s = struct('Description', 'Recorded duration of period 1', 'Units', 's');
J.cond8_recorded_s = struct('Description', 'Recorded duration of period 8', 'Units', 's');
write_json(fullfile(cfg.bids_root, 'participants.json'), J);
end

% -------------------------------------------------------------------------
function write_tsv(fn, header, rows)
% WRITE_TSV  Tab-separated file; empty / NaN / missing values -> 'n/a'
f = fopen(fn, 'w');
fprintf(f, '%s\n', strjoin(header, '\t'));
for i = 1:size(rows, 1)
    cells = cell(1, size(rows, 2));
    for j = 1:size(rows, 2)
        v = rows{i, j};
        if isstring(v), v = char(v); end
        if isempty(v) || (isnumeric(v) && any(isnan(v))) || (ischar(v) && strcmp(v, '<missing>'))
            cells{j} = 'n/a';
        elseif isnumeric(v) || islogical(v)
            cells{j} = num2str(double(v), '%.6g');
        else
            cells{j} = v;
        end
    end
    fprintf(f, '%s\n', strjoin(cells, '\t'));
end
fclose(f);
end

% -------------------------------------------------------------------------
function write_json(fn, S)
f = fopen(fn, 'w');
fprintf(f, '%s\n', jsonencode(S, 'PrettyPrint', true));
fclose(f);
end

% -------------------------------------------------------------------------
function make_dir(d)
if ~isfolder(d), mkdir(d); end
end

% -------------------------------------------------------------------------
function T = cell_structs_to_table(L)
% Stack structs that may have different fields (missing ones -> "n/a"/NaN)
fields = {};
for i = 1:numel(L), fields = union(fields, fieldnames(L{i}), 'stable'); end
proto = struct();
for f = fields(:)'
    for i = 1:numel(L)
        if isfield(L{i}, f{1}), proto.(f{1}) = L{i}.(f{1}); break; end
    end
end
for i = 1:numel(L)
    for f = fields(:)'
        if ~isfield(L{i}, f{1})
            p = proto.(f{1});
            if isnumeric(p), L{i}.(f{1}) = NaN; else, L{i}.(f{1}) = "n/a"; end
        end
    end
    L{i} = orderfields(L{i}, fields);
end
T = struct2table([L{:}]', 'AsArray', true);
end
