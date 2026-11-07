% =========================================================================
% preprocess_rs_eeg.m
% =========================================================================
%
% PURPOSE
%   Automated preprocessing of the resting-state EEG (eyes open / eyes
%   closed) for SPECTRAL analyses (theta, alpha + IAF, beta, theta/beta and
%   alpha/beta ratios). Input: the base datasets written by
%   import_bdf_to_bids.m (BIDS/derivatives/eeglab-import).
%
%   Designed for THIS dataset, where electrode bridges are present in every
%   subject (median 49 channels): no subject is discarded because of
%   bridges; bridged channels are handled explicitly (see step 4).
%
% PIPELINE (one subject at a time; parameters in the CONFIGURATION section)
%   1.  Load the base dataset (512 Hz, average reference, corrected events).
%   2.  Two copies of the data:
%         - ICA copy      : high-pass 1 Hz   (better ICA decomposition;
%                           Winkler et al. 2015; Klug & Gramann 2021)
%         - analysis copy : high-pass 0.5 Hz (keeps 1-4 Hz for future use)
%       Zero-phase FIR filters (pop_eegfiltnew; Widmann et al. 2015).
%       No low-pass and no notch: the bands of interest are <= 30 Hz and the
%       50 Hz line noise of this dataset is small; ICLabel needs the
%       broadband signal to recognise muscle activity.
%   3.  Bad channels (scalp EEG only), on the ICA copy, with clean_rawdata
%       (Kothe & Makeig; parameters of Pernet et al. 2021 and DISCOVER-EEG,
%       Gil Avila et al. 2023): flat > 5 s, line-noise z > 4, correlation
%       with the prediction from neighbours < 0.8. Saturated channels (DC
%       offset at the amplifier limit, from the import) are always bad.
%       Ocular protection: channels flagged only because of eye activity
%       (|r| with VEOG/HEOG > 0.4) are kept, so that ICA can find and remove
%       the eye components (otherwise frontal channels are interpolated and
%       blinks survive, as happened in the first run for sub-34).
%   4.  Electrode bridges (Alschuler et al. 2014 electrical-distance method,
%       same implementation as qc_dataset_overview_v2.m). Channels bridged
%       together carry (almost) the same signal, so each GROUP of mutually
%       bridged channels is ONE independent signal. For ICA, only one
%       channel per group (the one closest to the group centre) is kept;
%       the others are temporarily removed. This removes the rank
%       deficiency that would otherwise create "ghost" components
%       (Kim et al. 2023). After cleaning, the removed members of a group
%       receive the cleaned signal of their representative (they recorded
%       that same signal). Bridged channels are NOT interpolated from
%       neighbours that are often bridged themselves.
%   5.  ICA training data: ICA copy, good + representative channels,
%       average reference, bursts removed with ASR (burst criterion 20,
%       segments removed, not reconstructed; Chang et al. 2020). ICA:
%       extended Infomax (runica; Picard optional), with the number of
%       components = data rank (n channels - 1, average reference).
%   6.  ICLabel (Pion-Tonachini et al. 2019): components classified as Eye
%       or Muscle with probability >= 0.80 are removed (Pernet et al. 2021;
%       DISCOVER-EEG). Heart / line / channel-noise removal can be switched
%       on in the configuration.
%   7.  The ICA weights are applied to the analysis copy (same channels,
%       same reference) and the artefactual components are subtracted.
%   8.  Full montage rebuilt: bridged-group members receive the cleaned
%       representative signal; bad channels are interpolated with spherical
%       splines (Perrin et al. 1989) from the cleaned data.
%   9.  Average reference of the 128 scalp channels.
%   10. Epochs for spectral analysis: inside every EO / EC period, from
%       cfg.skip_after_instruction_s after the audio instruction to the
%       end of the period, 2-s epochs with 50 % overlap (DISCOVER-EEG).
%       Epochs with any scalp channel beyond +/- cfg.epoch_reject_uV in the
%       1-45 Hz band are rejected.
%
% OUTPUTS (BIDS/derivatives/preproc-spectral/sub-XX/eeg/)
%   sub-XX_task-rest_desc-clean_eeg.set     continuous cleaned data
%   sub-XX_task-rest_desc-epochsEO_eeg.set  clean 2-s epochs, eyes open
%   sub-XX_task-rest_desc-epochsEC_eeg.set  clean 2-s epochs, eyes closed
%   sub-XX_task-rest_desc-qc.png            QC figure (see MANUAL CHECK)
%   sub-XX_task-rest_desc-ica.mat           ICA weights, ICLabel results
%   sub-XX_task-rest_desc-provenance.mat    everything done to the subject
%                                           (channels, bridges, ICs, epochs,
%                                           pre/post metrics) for tables and
%                                           teaching figures
% and in the derivative root:
%   preproc_summary.csv  one row per subject, written while running
%   preproc_log.txt      details
%   table_subjects.csv, table_channels.csv, table_ics.csv, table_bridges.csv,
%   table_epochs.csv, table_prepost.csv  descriptive tables for papers and
%                        theses (rebuilt at the end from the provenance files)
%   tables_README.txt    meaning of every column
%
% MANUAL CHECK (about 1 minute per subject)
%   Open preproc_summary.csv and look at the column "flags". For every
%   flagged subject, open its QC figure (see the instructions at the end of
%   this header). Unflagged subjects need no manual work.
%
% REQUIREMENTS
%   MATLAB R2026a, Signal Processing Toolbox, EEGLAB with the plugins
%   clean_rawdata, ICLabel and firfilt (all installed by default in EEGLAB
%   2023+), optionally picard (faster ICA: File > Manage EEGLAB extensions).
%   natural_sort_files.m on the path.
%
% RUN TIME
%   ~15-30 min per subject with runica (all 37: run it overnight),
%   ~3-6 min with Picard.
%   Subjects already processed are skipped (cfg.overwrite = false), so the
%   script can be stopped and restarted.
%
% -------------------------------------------------------------------------
% HOW TO READ THE QC FIGURE (only for flagged subjects)
%   Top row   : channel map. Red = bad (interpolated), blue = bridged,
%               black dot = representative of a bridged group.
%   Middle    : removed components (scalp map + ICLabel class). Eye
%               components look frontal (blinks: symmetric frontal; lateral
%               movements: left-right opposite poles). Muscle components
%               are small, at the edge (temporal / neck) with a spectrum
%               rising above 20 Hz. If a removed component looks like a
%               broad posterior map with an alpha peak, write it down.
%   Bottom    : mean spectrum (all scalp channels) before / after cleaning
%               for EO and EC. After cleaning the curves should be smoother
%               above 20 Hz (less muscle) and below 4 Hz (less eye), and
%               the alpha peak must still be there (larger in EC).
% =========================================================================

clear; clc;

%% ========================================================================
%  0. CONFIGURATION
%  ========================================================================
cfg = struct();
base = 'C:\Users\chiar\Documents\MATLAB\FONDECYT\DATOS\Resting State';
cfg.in_root   = fullfile(base, 'BIDS', 'derivatives', 'eeglab-import');
cfg.out_root  = fullfile(base, 'BIDS', 'derivatives', 'preproc-spectral');
cfg.overwrite = false;

cfg.n_eeg = 128;                    % scalp channels are 1..128

% Filters
cfg.hp_ica_hz      = 1.0;           % ICA copy
cfg.hp_analysis_hz = 0.5;           % analysis copy

% Bad channels (clean_rawdata)
cfg.flatline_s     = 5;
cfg.chan_corr      = 0.8;
cfg.line_noise_z   = 4;
cfg.max_bad_frac   = 0.15;          % more than 15 % bad channels -> flag
% Ocular protection: blinks / eye movements can make frontal channels look
% "unpredictable" to clean_rawdata. A channel flagged ONLY by clean_rawdata
% (not saturated) whose signal correlates with the bipolar VEOG or HEOG
% above this value is kept: its ocular activity is removed by ICA instead.
cfg.protect_ocular = true;
cfg.ocular_r_thr   = 0.4;

% Bridges (Alschuler et al. 2014; same settings as qc_dataset_overview_v2.m)
cfg.bridge_seg_s        = 2;
cfg.bridge_bin          = 0.25;
cfg.bridge_interp_step  = 0.05;
cfg.bridge_peak_max     = 5;
cfg.bridge_cutoff_max   = 10;
cfg.bridge_min_seg_frac = 0.5;
cfg.bridge_valley_max   = 0.5;
cfg.bridge_rel_ed_max   = 0.05;
cfg.bridge_lp_hz        = 40;       % ED computed on 1-40 Hz data

% ICA
% 'runica' = extended Infomax (used by Pernet et al. 2021 and DISCOVER-EEG;
% ~15-30 min per subject). 'picard' = same model, much faster, but needs
% the picard plugin (File > Manage EEGLAB extensions).
cfg.ica_algorithm  = 'runica';
cfg.asr_burst      = 20;            % ASR burst criterion (training data only)
cfg.ica_seed       = 42;            % reproducible ICA
cfg.min_ica_chans  = 30;            % fewer independent channels -> flag

% ICLabel: class -> probability threshold (Inf = never removed)
cfg.iclabel_thr = struct('Muscle', 0.80, 'Eye', 0.80, 'Heart', Inf, ...
                         'LineNoise', Inf, 'ChannelNoise', Inf);
cfg.max_ic_removed_frac = 0.30;     % more than 30 % of ICs removed -> flag

% Epochs
cfg.skip_after_instruction_s = 4;   % discard the first 4 s of each period
cfg.epoch_len_s     = 2;
cfg.epoch_step_s    = 1;            % 50 % overlap
cfg.epoch_reject_uV = 100;
% The amplitude criterion is evaluated on a 1-45 Hz filtered copy, so that
% slow drifts (sweat; e.g. sub-32, sub-39) below the bands of interest do
% not reject epochs. The saved epochs keep the 0.5 Hz high-pass data.
cfg.reject_band_hz  = [1 45];
cfg.min_epochs      = 30;           % fewer clean epochs per condition -> flag

cfg.save_figure = true;

%% ========================================================================
%  1. SETUP
%  ========================================================================
if ~exist('pop_clean_rawdata', 'file') || ~exist('iclabel', 'file')
    if exist('eeglab', 'file'), eeglab('nogui'); else, error('EEGLAB not found.'); end
end
assert(exist('pop_clean_rawdata', 'file') == 2, 'clean_rawdata plugin missing.');
assert(exist('iclabel', 'file') == 2, 'ICLabel plugin missing.');
assert(exist('pop_eegfiltnew', 'file') == 2, 'firfilt plugin missing.');
if strcmpi(cfg.ica_algorithm, 'picard') && exist('picard', 'file') ~= 2
    warning('Picard not installed: using runica.');
    cfg.ica_algorithm = 'runica';
end
fprintf('ICA algorithm: %s\n', cfg.ica_algorithm);

if ~isfolder(cfg.in_root), error('Input folder not found:\n  %s', cfg.in_root); end
if ~isfolder(cfg.out_root), mkdir(cfg.out_root); end
save(fullfile(cfg.out_root, 'preproc_config.mat'), 'cfg');

files = dir(fullfile(cfg.in_root, 'sub-*', 'eeg', '*_desc-import_eeg.set'));
files = natural_sort_files(files);
nS = numel(files);
fprintf('Found %d subject(s).\n', nS);

log_fid = fopen(fullfile(cfg.out_root, 'preproc_log.txt'), 'a');
fprintf(log_fid, '\n%s\npreprocess_rs_eeg.m  -  %s\n%s\n', repmat('=', 1, 78), ...
        char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm')), repmat('=', 1, 78));

summary_file = fullfile(cfg.out_root, 'preproc_summary.csv');
t_all = tic;                                        % total run time

%% ========================================================================
%  2. LOOP OVER SUBJECTS
%  ========================================================================
for iS = 1:nS
    sub  = extractBefore(files(iS).name, '_task');
    stem = [sub '_task-rest'];
    out_dir = fullfile(cfg.out_root, sub, 'eeg');
    clean_name = [stem '_desc-clean_eeg.set'];
    fprintf('\n[%d/%d] %s\n', iS, nS, sub);

    if ~cfg.overwrite && isfile(fullfile(out_dir, clean_name))
        fprintf('  already processed: skipped\n');
        continue
    end
    if ~isfolder(out_dir), mkdir(out_dir); end
    fprintf(log_fid, '\n--- %s ---\n', sub);

    t_sub = tic;                                    % run time of this subject
    R = struct('subject', string(sub), 'status', "OK", 'message', "", 'flags', "");
    flags = strings(0);
    try
        %% 1. Load
        EEG0 = pop_loadset('filename', files(iS).name, 'filepath', files(iS).folder);
        scalp = 1:cfg.n_eeg;
        labels = {EEG0.chanlocs.labels};

        %% 2. Filtered copies (all 136 channels; EXT are filtered but not cleaned)
        EEGi = pop_eegfiltnew(EEG0, cfg.hp_ica_hz, []);       % ICA copy
        EEGa = pop_eegfiltnew(EEG0, cfg.hp_analysis_hz, []);  % analysis copy

        %% 3. Bad channels (scalp only)
        sat = {};
        if isfield(EEG0.etc, 'import') && isfield(EEG0.etc.import, 'saturated_channels')
            sat = EEG0.etc.import.saturated_channels;
        end
        tmp = pop_select(EEGi, 'channel', scalp);
        tmp = pop_clean_rawdata(tmp, 'FlatlineCriterion', cfg.flatline_s, ...
            'ChannelCriterion', cfg.chan_corr, 'LineNoiseCriterion', cfg.line_noise_z, ...
            'Highpass', 'off', 'BurstCriterion', 'off', 'WindowCriterion', 'off', ...
            'BurstRejection', 'off', 'Distance', 'Euclidian');
        kept = {tmp.chanlocs.labels};
        clear tmp
        bad = reshape(setdiff(labels(scalp), kept, 'stable'), 1, []);
        bad_cleanrawdata = bad;                         % before ocular protection
        r_oc = ocular_correlation(EEGi, scalp);         % saved for the tables

        % Ocular protection (see cfg.protect_ocular)
        rescued = {};
        if cfg.protect_ocular
            cand = ismember(labels(scalp), bad) & ~ismember(labels(scalp), sat);
            resc = cand & (r_oc(:)' > cfg.ocular_r_thr);
            rescued = labels(scalp(resc));
            bad = setdiff(bad, rescued, 'stable');
        end
        R.n_rescued_ocular = numel(rescued);
        R.rescued_channels = string(strjoin(rescued, ' '));

        bad = unique([bad, reshape(intersect(sat, labels(scalp)), 1, [])], 'stable');
        bad_idx = find(ismember(labels, bad));
        R.n_bad = numel(bad_idx);
        R.bad_channels = string(strjoin(bad, ' '));
        if R.n_bad > cfg.max_bad_frac * cfg.n_eeg
            flags(end+1) = sprintf('%d bad channels', R.n_bad);
        end

        %% 4. Bridges among the good scalp channels
        good = setdiff(scalp, bad_idx);
        tmp  = pop_select(EEGi, 'channel', good);
        tmp  = pop_eegfiltnew(tmp, [], cfg.bridge_lp_hz);
        br   = detect_bridges_ed(double(tmp.data), tmp.srate, true(numel(good), 1), cfg);
        clear tmp
        pairs = good(br.pairs);                         % indices in the full montage
        if isvector(pairs) && numel(pairs) == 2, pairs = pairs(:)'; end
        [groups, rep] = bridge_groups(pairs, EEG0.chanlocs);
        members = setdiff([groups{:}], rep);            % removed for ICA
        R.n_bridged = numel([groups{:}]);
        R.n_bridge_groups = numel(groups);
        R.largest_group = max([0, cellfun(@numel, groups)]);
        R.bridge_note = string(br.note);

        %% 5. ICA training data
        ica_ch = setdiff(scalp, [bad_idx, members]);    % independent channels
        R.n_ica_channels = numel(ica_ch);
        if numel(ica_ch) < cfg.min_ica_chans
            flags(end+1) = sprintf('only %d independent channels for ICA', numel(ica_ch));
        end
        % Average reference of the independent channels, applied to ALL
        % channels of both copies (keeps every channel in the same reference)
        EEGi.data = EEGi.data - mean(EEGi.data(ica_ch, :), 1);
        EEGa.data = EEGa.data - mean(EEGa.data(ica_ch, :), 1);
        T = pop_select(EEGi, 'channel', ica_ch);
        T.ref  = 'average';
        T = pop_clean_rawdata(T, 'FlatlineCriterion', 'off', 'ChannelCriterion', 'off', ...
            'LineNoiseCriterion', 'off', 'Highpass', 'off', 'BurstCriterion', cfg.asr_burst, ...
            'WindowCriterion', 'off', 'BurstRejection', 'on', 'Distance', 'Euclidian');
        R.ica_training_s = T.pnts / T.srate;

        rk = data_rank(T.data);
        R.ica_rank = rk;
        rng(cfg.ica_seed);
        if strcmpi(cfg.ica_algorithm, 'picard')
            T = pop_runica(T, 'icatype', 'picard', 'pca', rk);
        else
            T = pop_runica(T, 'icatype', 'runica', 'extended', 1, 'pca', rk);
        end

        % Segments removed by ASR from the ICA training data (for the tables)
        asr_removed_s = zeros(0, 2);
        if isfield(T.etc, 'clean_sample_mask')
            msk = T.etc.clean_sample_mask(:)';
            d = diff([true msk true]);
            asr_removed_s = [find(d == -1)' (find(d == 1) - 1)'] / T.srate;
        end
        R.asr_removed_s = sum(diff(asr_removed_s, 1, 2));

        %% 6. ICLabel
        T = iclabel(T);
        pvaf = ic_pvaf(T);                              % % variance of each IC
        cls  = T.etc.ic_classification.ICLabel.classes;
        prob = T.etc.ic_classification.ICLabel.classifications;
        rej  = false(size(prob, 1), 1);
        thr_names = fieldnames(cfg.iclabel_thr);
        for k = 1:numel(thr_names)
            c = find(strcmpi(strrep(cls, ' ', ''), thr_names{k}));
            if ~isempty(c)
                rej = rej | prob(:, c) >= cfg.iclabel_thr.(thr_names{k});
            end
        end
        rej_idx = find(rej)';
        [~, top] = max(prob, [], 2);
        R.n_ics = size(prob, 1);
        R.n_ics_removed = numel(rej_idx);
        R.n_eye_removed = sum(rej & top == find(strcmp(cls, 'Eye')));
        R.n_muscle_removed = sum(rej & top == find(strcmp(cls, 'Muscle')));
        R.n_brain_ics = sum(top == find(strcmp(cls, 'Brain')));
        if R.n_ics_removed > cfg.max_ic_removed_frac * R.n_ics
            flags(end+1) = sprintf('%d of %d ICs removed', R.n_ics_removed, R.n_ics);
        end

        %% 7. Apply ICA to the analysis copy
        EEGa_before = EEGa;                             % for the QC figure
        A = pop_select(EEGa, 'channel', ica_ch);        % same channels and reference as T
        A.icaweights  = T.icaweights;
        A.icasphere   = T.icasphere;
        A.icawinv     = T.icawinv;
        A.icachansind = 1:A.nbchan;
        A = eeg_checkset(A, 'ica');
        var_before = mean(var(double(A.data), 0, 2));
        if ~isempty(rej_idx)
            A = pop_subcomp(A, rej_idx, 0);
        end
        R.pct_variance_removed = 100 * (1 - mean(var(double(A.data), 0, 2)) / var_before);

        %% 8. Rebuild the full montage
        C = EEGa;                                       % 136 channels, analysis copy
        C.data(ica_ch, :) = A.data;
        for g = 1:numel(groups)                         % bridged members <- representative
            r = rep(g);
            for m = setdiff(groups{g}, r)
                C.data(m, :) = C.data(r, :);
            end
        end
        C.icaweights = []; C.icasphere = []; C.icawinv = []; C.icaact = []; C.icachansind = [];
        % Interpolation from the SCALP channels only (the external channels
        % have positions in the file but must not be used as sources)
        if ~isempty(bad_idx)
            Sc = pop_select(C, 'channel', scalp);
            Sc = pop_interp(Sc, bad_idx, 'spherical');
            C.data(scalp, :) = Sc.data;
            clear Sc
        end

        %% 9. Average reference of the 128 scalp channels (all channels)
        C.data = C.data - mean(C.data(scalp, :), 1);
        C.ref  = 'average of the 128 scalp channels';
        C = eeg_checkset(C);

        % Documentation
        C.etc.preproc.script        = 'preprocess_rs_eeg.m';
        C.etc.preproc.date          = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm'));
        C.etc.preproc.cfg           = cfg;
        C.etc.preproc.bad_channels  = bad;
        C.etc.preproc.bridge_groups = cellfun(@(g) labels(g), groups, 'UniformOutput', false);
        C.etc.preproc.bridge_representatives = labels(rep);
        C.etc.preproc.ica_channels  = labels(ica_ch);
        C.etc.preproc.ica_rank      = rk;
        C.etc.preproc.ic_removed    = rej_idx;
        C.etc.preproc.effective_rank = numel(ica_ch) - 1 - numel(rej_idx);
        C.etc.preproc.note = ['Bridged-group members carry the cleaned signal of their representative; ' ...
            'bad channels are interpolated. Effective rank = independent channels - 1 - removed ICs.'];
        C.comments = sprintf('%s\nPreprocessed by preprocess_rs_eeg.m on %s', C.comments, C.etc.preproc.date);
        R.effective_rank = C.etc.preproc.effective_rank;

        %% Per-channel description and pre / post metrics (tables)
        chanT = channel_table(labels, scalp, sat, bad_cleanrawdata, rescued, bad_idx, ...
                              r_oc, groups, rep, EEGa_before, C, cfg);
        prepost = prepost_table(EEGa_before, C, cfg);
        R.n_interpolated     = numel(bad_idx);
        R.n_bad_saturated    = numel(intersect(sat, labels(scalp)));
        R.n_bad_cleanrawdata = numel(setdiff(bad_cleanrawdata, rescued));
        R.n_bridged_copied   = numel(members);
        R.n_valid_connectivity = sum(chanT.valid_for_connectivity);
        R.pvaf_removed       = sum(pvaf(rej_idx));
        R.median_sd_pre_uV   = median(chanT.sd_pre_uV);
        R.median_sd_post_uV  = median(chanT.sd_post_uV);

        C.setname = [stem '_desc-clean'];
        pop_saveset(C, 'filename', clean_name, 'filepath', out_dir, 'savemode', 'twofiles');

        ica = struct('weights', T.icaweights, 'sphere', T.icasphere, 'winv', T.icawinv, ...
                     'chan_labels', {labels(ica_ch)}, 'iclabel_classes', {cls}, ...
                     'iclabel_prob', prob, 'removed', rej_idx, 'rank', rk);
        save(fullfile(out_dir, [stem '_desc-ica.mat']), '-struct', 'ica');

        %% 10. Clean epochs per condition
        n_ep = struct('EO', 0, 'EC', 0, 'EO_total', 0, 'EC_total', 0);
        epT = table();
        Cf = pop_eegfiltnew(C, cfg.reject_band_hz(1), cfg.reject_band_hz(2));   % rejection copy
        for cond = {'EO', 'EC'}
            [E, n_tot, seg_t, seg_ci] = make_epochs(C, cond{1}, cfg);
            if isempty(E)
                flags(end+1) = sprintf('no %s epochs', cond{1}); %#ok<SAGROW>
                continue
            end
            Ef = make_epochs(Cf, cond{1}, cfg);            % same segments, filtered
            mx = squeeze(max(max(abs(Ef.data(scalp, :, :)), [], 1), [], 2));
            clear Ef
            if numel(mx) == numel(seg_t)
                epT = [epT; table(repmat(string(cond{1}), numel(seg_t), 1), seg_ci(:), ...
                    seg_t(:), mx(:), mx(:) > cfg.epoch_reject_uV, ...
                    'VariableNames', {'condition','condition_index','onset_s', ...
                    'max_abs_uV_1_45Hz','rejected'})]; %#ok<AGROW>
            end
            badep = find(mx > cfg.epoch_reject_uV);
            if ~isempty(badep) && numel(badep) < E.trials
                E = pop_rejepoch(E, badep, 0);
            elseif numel(badep) == E.trials
                flags(end+1) = sprintf('all %s epochs exceed %g uV', cond{1}, cfg.epoch_reject_uV); %#ok<SAGROW>
            end
            n_ep.(cond{1}) = E.trials;
            n_ep.([cond{1} '_total']) = n_tot;
            if E.trials < cfg.min_epochs
                flags(end+1) = sprintf('%d clean %s epochs', E.trials, cond{1}); %#ok<SAGROW>
            end
            E.setname = [stem '_desc-epochs' cond{1}];
            pop_saveset(E, 'filename', [E.setname '_eeg.set'], 'filepath', out_dir, 'savemode', 'twofiles');
        end
        R.epochs_EO = n_ep.EO; R.epochs_EO_total = n_ep.EO_total;
        R.epochs_EC = n_ep.EC; R.epochs_EC_total = n_ep.EC_total;
        R.pct_epochs_rejected = 100 * (1 - (R.epochs_EO + R.epochs_EC) / ...
                                max(1, R.epochs_EO_total + R.epochs_EC_total));

        % ---------- Provenance file: everything needed for tables / teaching
        icT = table((1:size(prob, 1))', string(cls(top))', max(prob, [], 2), pvaf(:), ...
                    ismember((1:size(prob, 1))', rej_idx(:)), ...
                    'VariableNames', {'ic','iclabel_class','iclabel_prob','pvaf','removed'});
        for c = 1:numel(cls)
            icT.(['p_' strrep(cls{c}, ' ', '')]) = prob(:, c);
        end
        brT = table();
        if ~isempty(br.pairs)
            brT = table(string(labels(good(br.pairs(:, 1))))', string(labels(good(br.pairs(:, 2))))', ...
                br.pair_med_ed(:), br.pair_frac(:), ...
                'VariableNames', {'chan_1','chan_2','median_ED_uV2','frac_segments_bridged'});
        end
        P = struct();
        P.subject = sub;  P.summary = R;  P.cfg = cfg;
        P.channels = chanT;  P.ics = icT;  P.bridges = brT;  P.epochs = epT;
        P.prepost = prepost;
        P.bridge_hist = struct('scaled_ED', br.hist_x, 'density', br.hist_y, ...
                               'cutoff', br.cutoff, 'median_ED_uV2', br.scale, 'note', br.note);
        P.bridge_groups = cellfun(@(g) labels(g), groups, 'UniformOutput', false);
        P.bridge_representatives = labels(rep);
        P.asr_removed_s = asr_removed_s;
        P.ica_channels = labels(ica_ch);
        save(fullfile(out_dir, [stem '_desc-provenance.mat']), '-struct', 'P');

        %% QC figure
        if cfg.save_figure
            try
                qc_figure(EEGa_before, C, T, rej_idx, cls, prob, bad_idx, groups, rep, cfg, ...
                          fullfile(out_dir, [stem '_desc-qc.png']), sub);
            catch MEf
                flags(end+1) = "figure not saved: " + MEf.message; %#ok<SAGROW>
            end
        end

        fprintf(log_fid, ['bad (%d): %s\nbridges: %d channels in %d groups (largest %d); ' ...
            'ICA on %d channels, rank %d, %.0f s; ICs removed %d/%d (eye %d, muscle %d); ' ...
            'variance removed %.1f%%; epochs EO %d/%d, EC %d/%d\n'], R.n_bad, R.bad_channels, ...
            R.n_bridged, R.n_bridge_groups, R.largest_group, R.n_ica_channels, rk, R.ica_training_s, ...
            R.n_ics_removed, R.n_ics, R.n_eye_removed, R.n_muscle_removed, R.pct_variance_removed, ...
            R.epochs_EO, R.epochs_EO_total, R.epochs_EC, R.epochs_EC_total);

    catch ME
        R.status  = "ERROR";
        R.message = string(ME.message);
        if ~isempty(ME.stack)
            R.message = R.message + sprintf(' (%s, line %d)', ME.stack(1).name, ME.stack(1).line);
        end
        fprintf(2, '  ERROR: %s\n', R.message);
        fprintf(log_fid, 'ERROR: %s\n', R.message);
    end
    R.flags = strjoin(flags, ' | ');
    if strlength(R.flags) > 0
        fprintf(log_fid, 'FLAGS: %s\n', R.flags);
        fprintf('  FLAGS: %s\n', R.flags);
    end
    R.minutes = toc(t_sub) / 60;
    fprintf(log_fid, 'time: %.1f min\n', R.minutes);
    append_summary(summary_file, R);
    % Store the final summary (with flags and run time) in the provenance file
    prov = fullfile(out_dir, [stem '_desc-provenance.mat']);
    if R.status == "OK" && isfile(prov)
        upd = struct('summary', R);
        save(prov, '-struct', 'upd', '-append');
    end
    clear EEG0 EEGi EEGa EEGa_before T A C Cf E
end

% Group-level descriptive tables (rebuilt from all provenance files, so
% they are always complete, also after an interrupted / resumed run)
build_group_tables(cfg.out_root);

fprintf(log_fid, '\nTOTAL RUN TIME: %.2f h (finished %s)\n', toc(t_all) / 3600, ...
        char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm')));
fclose(log_fid);
fprintf('\nDone in %.2f h. Outputs in:\n  %s\n', toc(t_all) / 3600, cfg.out_root);


%% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================

% -------------------------------------------------------------------------
function pvaf = ic_pvaf(T)
% IC_PVAF  Percent of the variance of the ICA training data accounted for
% by each component (back-projection of one component at a time).
X = double(T.data);
X = X - mean(X, 2);
act = (T.icaweights * T.icasphere) * X;
v0 = sum(var(X, 0, 2));
pvaf = zeros(size(T.icawinv, 2), 1);
for k = 1:numel(pvaf)
    proj = T.icawinv(:, k) * act(k, :);
    pvaf(k) = 100 * (1 - sum(var(X - proj, 0, 2)) / v0);
end
end

% -------------------------------------------------------------------------
function chanT = channel_table(labels, scalp, sat, bad_cr, rescued, bad_idx, ...
                               r_oc, groups, rep, B, C, cfg)
% CHANNEL_TABLE  One row per scalp channel: what was found and what was done.
%   handling:
%     good                    measured signal, used as is
%     bridged_representative  measured signal, one per bridged group (used in ICA)
%     bridged_copy            bridged with another channel: receives the
%                             cleaned signal of its representative
%     interpolated_saturated  saturated electrode, spherical-spline interpolation
%     interpolated_bad        bad channel (clean_rawdata), interpolation
%   valid_for_connectivity: true only for independent measured signals
%   (good + bridged_representative). All channels are usable for regional
%   spectral power and topographies.
n   = numel(scalp);
lab = string(labels(scalp))';
is_sat  = ismember(lab, string(sat));
is_cr   = ismember(lab, string(bad_cr));
is_resc = ismember(lab, string(rescued));
is_bad  = ismember(scalp(:), bad_idx(:));
group_id = zeros(n, 1);
is_rep   = false(n, 1);
copied_from = strings(n, 1);
for g = 1:numel(groups)
    group_id(groups{g}) = g;
    is_rep(rep(g)) = true;
    for m = setdiff(groups{g}, rep(g))
        copied_from(m) = labels{rep(g)};
    end
end
handling = repmat("good", n, 1);
handling(group_id > 0 &  is_rep) = "bridged_representative";
handling(group_id > 0 & ~is_rep) = "bridged_copy";
handling(is_bad & ~is_sat) = "interpolated_bad";
handling(is_bad &  is_sat) = "interpolated_saturated";
bad_reason = strings(n, 1);
bad_reason(is_cr)  = "clean_rawdata";
bad_reason(is_sat) = "saturated";
valid_conn = ismember(handling, ["good", "bridged_representative"]);

% Standard deviation (1-45 Hz, average reference) before / after cleaning
Bs = pop_select(B, 'channel', scalp);
Bs = pop_eegfiltnew(Bs, 1, 45);
X  = double(Bs.data);  X = X - mean(X, 1);
sd_pre = std(X, 0, 2);
clear Bs X
Cs = pop_select(C, 'channel', scalp);
Cs = pop_eegfiltnew(Cs, 1, 45);
sd_post = std(double(Cs.data), 0, 2);
clear Cs

chanT = table(lab, is_bad, bad_reason, is_resc, r_oc(:), group_id, is_rep, copied_from, ...
    handling, valid_conn, sd_pre, sd_post, ...
    'VariableNames', {'channel','interpolated','bad_reason','ocular_rescued','r_ocular', ...
    'bridge_group','bridge_representative','copied_from','handling', ...
    'valid_for_connectivity','sd_pre_uV','sd_post_uV'});
end

% -------------------------------------------------------------------------
function T = prepost_table(B, C, cfg)
% PREPOST_TABLE  Band power of the mean scalp spectrum before / after the
% cleaning (ICA + interpolation), for EO and EC.
bands = {'delta', [1 4]; 'theta', [4 8]; 'alpha', [8 13]; 'beta', [13 30]; 'gamma', [30 45]};
T = table();
for cond = {'EO', 'EC'}
    [f, p0] = cond_spectrum(B, cond{1}, cfg);
    [~, p1] = cond_spectrum(C, cond{1}, cfg);
    if isempty(f), continue; end
    for b = 1:size(bands, 1)
        sel = f >= bands{b, 2}(1) & f <= bands{b, 2}(2);
        pw0 = trapz(f(sel), p0(sel));
        pw1 = trapz(f(sel), p1(sel));
        T = [T; table(string(cond{1}), string(bands{b, 1}), bands{b, 2}(1), bands{b, 2}(2), ...
             pw0, pw1, 100 * (pw1 - pw0) / pw0, ...
             'VariableNames', {'condition','band','f_low_hz','f_high_hz', ...
             'power_pre_uV2','power_post_uV2','change_pct'})]; %#ok<AGROW>
    end
end
end

% -------------------------------------------------------------------------
function build_group_tables(root)
% BUILD_GROUP_TABLES  Descriptive tables for papers / theses, rebuilt from
% the provenance file of every subject.
files = dir(fullfile(root, 'sub-*', 'eeg', '*_desc-provenance.mat'));
if isempty(files), return; end
S = {}; Ch = {}; Ic = {}; Br = {}; Ep = {}; Pp = {};
for i = 1:numel(files)
    P = load(fullfile(files(i).folder, files(i).name));
    sub = string(P.subject);
    S{end+1} = P.summary; %#ok<AGROW>
    Ch{end+1} = add_subject(P.channels, sub); %#ok<AGROW>
    Ic{end+1} = add_subject(P.ics, sub); %#ok<AGROW>
    Br{end+1} = add_subject(P.bridges, sub); %#ok<AGROW>
    Ep{end+1} = add_subject(P.epochs, sub); %#ok<AGROW>
    Pp{end+1} = add_subject(P.prepost, sub); %#ok<AGROW>
end
writetable(stack_structs(S), fullfile(root, 'table_subjects.csv'));
write_stack(Ch, fullfile(root, 'table_channels.csv'));
write_stack(Ic, fullfile(root, 'table_ics.csv'));
write_stack(Br, fullfile(root, 'table_bridges.csv'));
write_stack(Ep, fullfile(root, 'table_epochs.csv'));
write_stack(Pp, fullfile(root, 'table_prepost.csv'));
write_dictionary(fullfile(root, 'tables_README.txt'));
end

function T = add_subject(T, sub)
if isempty(T) || height(T) == 0, return; end
T.subject = repmat(sub, height(T), 1);
T = movevars(T, 'subject', 'Before', 1);
end

function write_stack(L, fn)
L = L(~cellfun(@(t) isempty(t) || height(t) == 0, L));
if ~isempty(L), writetable(vertcat(L{:}), fn); end
end

function T = stack_structs(L)
% Stack structs with possibly different fields (missing -> NaN or "")
fields = {};
for i = 1:numel(L), fields = union(fields, fieldnames(L{i}), 'stable'); end
for i = 1:numel(L)
    for f = fields(:)'
        if ~isfield(L{i}, f{1})
            proto = [];
            for j = 1:numel(L)
                if isfield(L{j}, f{1}), proto = L{j}.(f{1}); break; end
            end
            if isnumeric(proto) || islogical(proto), L{i}.(f{1}) = NaN; else, L{i}.(f{1}) = ""; end
        end
    end
    L{i} = orderfields(L{i}, fields);
end
T = struct2table([L{:}]', 'AsArray', true);
end

function write_dictionary(fn)
% Plain-text description of every table and column
txt = {
'DESCRIPTIVE TABLES - preprocess_rs_eeg.m'
'=========================================='
''
'table_subjects.csv  (one row per subject)'
'  n_bad / n_interpolated   scalp channels interpolated (saturated + clean_rawdata)'
'  n_bad_saturated          electrodes at the amplifier limit (no signal)'
'  n_bad_cleanrawdata       flat > 5 s, line noise z > 4 or low correlation with neighbours'
'  n_rescued_ocular         channels flagged only because of eye activity, kept for ICA'
'  n_bridged                channels involved in electrolyte bridges'
'  n_bridge_groups          groups of mutually bridged channels (1 signal per group)'
'  largest_group            channels in the largest bridged group'
'  n_bridged_copied         bridged channels that received their representative''s signal'
'  n_ica_channels           independent channels used for ICA'
'  n_valid_connectivity     channels with an independent measured signal'
'  ica_rank / n_ics         ICA dimension (independent channels - 1, average reference)'
'  n_ics_removed            components removed (n_eye_removed, n_muscle_removed)'
'  pvaf_removed             % variance of the ICA training data explained by removed ICs'
'  asr_removed_s            seconds removed by ASR from the ICA training data only'
'  epochs_EO / EC (_total)  clean 2-s epochs (out of the total available)'
'  pct_epochs_rejected      % epochs beyond +/-100 uV in the 1-45 Hz band'
'  median_sd_pre/post_uV    median channel SD (1-45 Hz) before / after cleaning'
'  flags                    reasons to look at the QC figure'
''
'table_channels.csv  (one row per subject x scalp channel)'
'  interpolated, bad_reason (saturated / clean_rawdata), ocular_rescued, r_ocular'
'  bridge_group (0 = not bridged), bridge_representative, copied_from'
'  handling: good | bridged_representative | bridged_copy | interpolated_bad | interpolated_saturated'
'  valid_for_connectivity: true = independent measured signal (good or representative)'
'  sd_pre_uV / sd_post_uV: SD 1-45 Hz before / after cleaning'
''
'table_ics.csv  (one row per subject x independent component)'
'  iclabel_class / iclabel_prob: most probable ICLabel class; p_<class>: all probabilities'
'  pvaf: % of variance of the training data; removed: true if subtracted'
''
'table_bridges.csv  (one row per bridged pair)'
'  median_ED_uV2: electrical distance (variance of the difference); frac_segments_bridged'
''
'table_epochs.csv  (one row per 2-s epoch)'
'  condition, condition_index (1-8), onset_s, max_abs_uV_1_45Hz, rejected'
''
'table_prepost.csv  (one row per subject x condition x band)'
'  power of the mean scalp spectrum before / after cleaning and % change'
};
fid = fopen(fn, 'w');
fprintf(fid, '%s\n', txt{:});
fclose(fid);
end

% -------------------------------------------------------------------------
function r = ocular_correlation(EEG, scalp)
% OCULAR_CORRELATION  Maximum |correlation| of every scalp channel with the
% bipolar VEOG (VEOG_U - VEOG_D) and HEOG (HEOG_L - HEOG_R).
L = {EEG.chanlocs.labels};
idx = @(lab) find(strcmpi(L, lab), 1);
r = zeros(numel(scalp), 1);
pairs = {{'VEOG_U', 'VEOG_D'}, {'HEOG_L', 'HEOG_R'}};
X = double(EEG.data(scalp, :));
X = X - mean(X, 2);
sx = sqrt(sum(X.^2, 2));
for p = 1:numel(pairs)
    a = idx(pairs{p}{1}); b = idx(pairs{p}{2});
    if isempty(a) || isempty(b), continue; end
    e = double(EEG.data(a, :) - EEG.data(b, :));
    e = e - mean(e);
    r = max(r, abs((X * e') ./ (sx * norm(e) + eps)));
end
end

% -------------------------------------------------------------------------
function [E, n_tot, seg_t, seg_ci] = make_epochs(C, cond, cfg)
% MAKE_EPOCHS  2-s epochs (50 % overlap) inside every period of a condition,
% starting cfg.skip_after_instruction_s after the instruction.
fs  = C.srate;
ev  = C.event;
typ = {ev.type};
on  = find(strcmp(typ, cond));
seg = struct('type', {}, 'latency', {}, 'duration', {});
seg_t = []; seg_ci = [];
for k = on
    if isfield(ev, 'condition_index'), ci = ev(k).condition_index; else, ci = NaN; end
    t0 = (ev(k).latency - 1) / fs + cfg.skip_after_instruction_s;
    t1 = (ev(k).latency - 1 + ev(k).duration) / fs;           % end of the period
    for t = t0 : cfg.epoch_step_s : (t1 - cfg.epoch_len_s)
        seg(end+1).type = [cond '_seg']; %#ok<AGROW>
        seg(end).latency  = t * fs + 1;
        seg(end).duration = 0;
        seg_t(end+1)  = t; %#ok<AGROW>
        seg_ci(end+1) = ci; %#ok<AGROW>
    end
end
n_tot = numel(seg);
if n_tot == 0, E = []; return; end
% Add the segment markers (only the needed fields) and epoch
C2 = C;
C2.event = struct('type', {seg.type}, 'latency', {seg.latency}, 'duration', {seg.duration});
C2.urevent = [];
C2 = eeg_checkset(C2, 'eventconsistency');
E = pop_epoch(C2, {[cond '_seg']}, [0 cfg.epoch_len_s - 1 / fs], 'epochinfo', 'yes');
E = eeg_checkset(E);
end

% -------------------------------------------------------------------------
function [groups, rep] = bridge_groups(pairs, chanlocs)
% BRIDGE_GROUPS  Connected groups of bridged channels and, for each group,
% the representative = member closest to the group centre.
groups = {}; rep = [];
if isempty(pairs), return; end
n = numel(chanlocs);
G = graph(pairs(:, 1), pairs(:, 2), ones(size(pairs, 1), 1), n);
comp = conncomp(G);
inv_ch = unique(pairs(:));
cc = unique(comp(inv_ch));
for c = cc(:)'
    m = find(comp == c);
    xyz = [[chanlocs(m).X]' [chanlocs(m).Y]' [chanlocs(m).Z]'];
    [~, i] = min(sum((xyz - mean(xyz, 1)).^2, 2));
    groups{end+1} = m; %#ok<AGROW>
    rep(end+1) = m(i); %#ok<AGROW>
end
end

% -------------------------------------------------------------------------
function r = data_rank(X)
% DATA_RANK  Numerical rank of channels x samples data (eigenvalues of the
% covariance above 1e-7 of the largest), never more than n - 1 (average
% reference).
X = double(X);
X = X - mean(X, 2);
ev = eig((X * X') / (size(X, 2) - 1));
r = sum(ev > max(ev) * 1e-7);
r = min(r, size(X, 1) - 1);
end

% -------------------------------------------------------------------------
function br = detect_bridges_ed(Xb, srate, usable, cfg)
% DETECT_BRIDGES_ED  Electrical-distance bridge detection (Alschuler et al.
% 2014) with two safeguards for 128-channel caps (see qc_dataset_overview_v2.m).
nE  = size(Xb, 1);
use = find(usable(:)');
br  = struct('pairs', zeros(0, 2), 'note', "no bridges", 'pair_med_ed', [], 'pair_frac', [], ...
             'hist_x', [], 'hist_y', [], 'cutoff', NaN, 'scale', NaN);
L   = round(cfg.bridge_seg_s * srate);
nSg = floor(size(Xb, 2) / L);
nU  = numel(use);
if nSg < 10 || nU < 2, return; end
iu  = find(triu(true(nU), 1));
EDs = zeros(numel(iu), nSg, 'single');
REL = zeros(numel(iu), nSg, 'single');
for k = 1:nSg
    seg = Xb(use, (k-1)*L+1 : k*L);
    seg = seg - mean(seg, 2);
    C   = (seg * seg') / (L - 1);
    v   = diag(C);
    E   = v + v' - 2 * C;
    Rl  = E ./ ((v + v') / 2 + eps);
    EDs(:, k) = E(iu);
    REL(:, k) = Rl(iu);
end
rel_ok = median(REL, 2) < cfg.bridge_rel_ed_max;
clear REL
br.scale = double(median(EDs(:)));
S = EDs * (100 / br.scale);
edges   = 0:cfg.bridge_bin:100;
counts  = histcounts(S(:), edges);
centers = edges(1:end-1) + cfg.bridge_bin / 2;
fine    = 0:cfg.bridge_interp_step:100;
dens    = interp1(centers, counts, fine, 'linear', 0);
br.hist_x = fine(fine <= 30);  br.hist_y = dens(fine <= 30);     % for the teaching figures
is_max  = islocalmax(dens) & fine > 0 & fine <= cfg.bridge_peak_max;
if ~any(is_max), br.note = "no near-zero ED peak"; return; end
cand = find(is_max);
[~, j] = max(dens(cand));
i_pk = cand(j);
i_mn = find(islocalmin(dens, 'FlatSelection', 'first') & (1:numel(fine)) > i_pk, 1, 'first');
if isempty(i_mn) || fine(i_mn) > cfg.bridge_cutoff_max || dens(i_mn) > cfg.bridge_valley_max * dens(i_pk)
    br.note = "near-zero ED peak not separated";
    return
end
frac  = mean(S <= fine(i_mn), 2);
is_br = frac >= cfg.bridge_min_seg_frac & rel_ok;
br.cutoff = fine(i_mn);
[r, c] = ind2sub([nU nU], iu(is_br));
br.pairs = [use(r(:))' use(c(:))'];
br.pair_med_ed = double(median(EDs(is_br, :), 2));
br.pair_frac   = double(frac(is_br));
br.note  = "bridges found";
if nE == 0, br.pairs = zeros(0, 2); end
end

% -------------------------------------------------------------------------
function qc_figure(EEG0, C, T, rej_idx, cls, prob, bad_idx, groups, rep, cfg, fn, sub)
% QC_FIGURE  Channel map, removed ICs, spectra before / after.
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [20 20 1500 1000]);
try, theme(fig, 'light'); catch, end
scalp = 1:cfg.n_eeg;
loc = C.chanlocs(scalp);

% (1) channel map
subplot(3, 6, 1:2);
bridged = [groups{:}];
cols = repmat([0.6 0.6 0.6], cfg.n_eeg, 1);
cols(intersect(bridged, scalp), :) = repmat([0.3 0.5 1], numel(intersect(bridged, scalp)), 1);
cols(bad_idx, :) = repmat([0.9 0.1 0.1], numel(bad_idx), 1);
[x, y] = topo_xy(loc);
scatter(y, x, 40, cols, 'filled'); hold on
if ~isempty(rep)
    [xr, yr] = topo_xy(C.chanlocs(rep));
    plot(yr, xr, 'k.', 'MarkerSize', 10);
end
axis equal off
title(sprintf('%s | red bad (%d), blue bridged (%d in %d groups)', sub, numel(bad_idx), ...
      numel(bridged), numel(groups)), 'Interpreter', 'none');

% (2) removed ICs (up to 10)
nshow = min(10, numel(rej_idx));
for k = 1:nshow
    subplot(3, 6, 2 + k);
    ic = rej_idx(k);
    topoplot(T.icawinv(:, ic), T.chanlocs, 'electrodes', 'off');
    [p, c] = max(prob(ic, :));
    title(sprintf('IC%d %s %.2f', ic, cls{c}, p), 'FontSize', 8);
end

% (3) spectra before / after, EO and EC
for ci = 1:2
    cond = {'EO', 'EC'};
    subplot(3, 6, 13 + 3 * (ci - 1) : 15 + 3 * (ci - 1)); hold on
    [f, p0] = cond_spectrum(EEG0, cond{ci}, cfg);
    [~, p1] = cond_spectrum(C, cond{ci}, cfg);
    if ~isempty(f)
        sel = f >= 1 & f <= 45;
        plot(f(sel), 10 * log10(p0(sel)), 'Color', [0.6 0.6 0.6], 'LineWidth', 1.5);
        plot(f(sel), 10 * log10(p1(sel)), 'k', 'LineWidth', 1.5);
        legend({'before ICA', 'after ICA'}, 'Location', 'northeast');
    end
    xlabel('Hz'); ylabel('dB (\muV^2/Hz)'); grid on
    title(sprintf('%s: mean spectrum of scalp channels', cond{ci}));
end
exportgraphics(fig, fn, 'Resolution', 90, 'BackgroundColor', 'white');
close(fig);
end

% -------------------------------------------------------------------------
function [f, p] = cond_spectrum(EEG, cond, cfg)
% Mean Welch spectrum over the scalp channels inside the periods of a
% condition (first cfg.skip_after_instruction_s s excluded).
fs = EEG.srate; f = []; p = [];
on = find(strcmp({EEG.event.type}, cond));
X = [];
for k = on
    s0 = round(EEG.event(k).latency + cfg.skip_after_instruction_s * fs);
    s1 = round(EEG.event(k).latency + EEG.event(k).duration) - 1;
    if s1 - s0 > 4 * fs
        X = [X, double(EEG.data(1:cfg.n_eeg, s0:s1))]; %#ok<AGROW>
    end
end
if isempty(X), return; end
X = X - mean(X(1:cfg.n_eeg, :), 1);      % average reference for comparison
[pp, f] = pwelch(X', hann(2 * fs), fs, 2 * fs, fs);
p = mean(pp, 2);
end

% -------------------------------------------------------------------------
function [x, y] = topo_xy(loc)
% 2-D projection for the channel map (EEGLAB polar coordinates)
th = [loc.theta] * pi / 180;
rd = [loc.radius];
x = rd .* cos(th);
y = rd .* sin(th);
end

% -------------------------------------------------------------------------
function append_summary(fn, R)
% Append one row to the summary CSV (header written the first time)
fields = {'subject','status','n_bad','bad_channels','n_rescued_ocular','rescued_channels', ...
          'n_bridged','n_bridge_groups', ...
          'largest_group','bridge_note','n_ica_channels','ica_rank','ica_training_s', ...
          'n_ics','n_ics_removed','n_eye_removed','n_muscle_removed','n_brain_ics', ...
          'pct_variance_removed','effective_rank','epochs_EO','epochs_EO_total', ...
          'epochs_EC','epochs_EC_total','pct_epochs_rejected','n_interpolated', ...
          'n_bad_saturated','n_bad_cleanrawdata','n_bridged_copied','n_valid_connectivity', ...
          'pvaf_removed','asr_removed_s','median_sd_pre_uV','median_sd_post_uV', ...
          'minutes','flags','message'};
new = ~isfile(fn);
fid = fopen(fn, 'a');
if new, fprintf(fid, '%s\n', strjoin(fields, ',')); end
vals = cell(1, numel(fields));
for i = 1:numel(fields)
    if isfield(R, fields{i}), v = R.(fields{i}); else, v = ""; end
    if isnumeric(v), v = num2str(v, '%.4g'); end
    v = char(string(v));
    vals{i} = ['"' strrep(v, '"', '''') '"'];
end
fprintf(fid, '%s\n', strjoin(vals, ','));
fclose(fid);
end

