% =========================================================================
% spectral_analysis_rs.m
% =========================================================================
%
% PURPOSE
%   Spectral analysis of the preprocessed resting-state EEG (output of
%   preprocess_rs_eeg.m) for every subject and for three conditions:
%     EO    = the 4 eyes-open periods
%     EC    = the 4 eyes-closed periods
%     total = all 8 periods (mean of the EO and EC spectra, equal weight)
%
%   Measures (all of them per subject x condition x region):
%     - Individual alpha frequency (IAF): peak (PAF) and centre of gravity (CoG)
%     - Theta, alpha (+ lower-1, lower-2, upper alpha) and beta power in
%       IAF-based bands (Klimesch 1999) AND in fixed bands (4-8, 8-13, 13-25 Hz)
%     - Ratios: TBR = theta/beta, ATR = alpha/theta, ABR = alpha/beta
%     - Aperiodic (1/f) component: offset and exponent (spectral
%       parameterization, Donoghue et al. 2020), and every band power and
%       ratio also in its PERIODIC form (power above the 1/f background)
%     - Alpha peak (centre frequency, height, width) and presence of theta
%       and beta peaks above the 1/f background
%     - Berger effect (alpha EC vs EO)
%     - Per-period values (each of the 8 periods separately) and blink rate
%       per period from the VEOG, to check eye closure and vigilance drifts
%       (used by sensitivity_rs.m for the case studies)
%
% VERSION 2 (2026-11-07), changes from version 1:
%   - beta band and aperiodic fit limited to <= 25 / 30 Hz (see point 6)
%   - per-period table and blink rate (table_by_period.csv)
%   - quality flags in table_wide.csv
%   - figures: separate colour scales per map, distinct band colours
% VERSION 2.1: eyes-closed flag based on the alpha increase of every EC
%   period (VEOG events also occur with the eyes closed)
%
% METHOD (reasons in the companion note "Spectral analysis - decisions")
%   1. Clean 4-s windows. The preprocessing produced 2-s epochs (50 %
%      overlap) and marked the rejected ones. A 4-s window is used only if
%      the three 2-s epochs it contains are clean and belong to the same
%      period; windows overlap by 50 %. 4 s gives 0.25 Hz resolution,
%      needed for the IAF (Corcoran et al. 2018).
%   2. Power spectrum: Hann-windowed periodogram of every window, combined
%      with the MEDIAN across windows (Welch-median, corrected for its bias).
%      The median is robust to residual transients (electrode "pops") that
%      ICA does not remove.
%   3. Narrow-band peaks (non-physiological, e.g. ~30 Hz in sub-19) are
%      detected on the mean spectrum and interpolated (+/- 1 Hz) before any
%      measure is computed.
%   4. Regions: 12 clusters of ~5 channels around the 10-20 equivalent
%      positions (Fz, F3, F4, Cz, C3, C4, Pz, P3, P4, Oz, O1, O2), computed on
%      an ideal spherical 10-20 model (channels within 20 deg); regions
%      frontal, central, parietal, occipital (union of 3 clusters each),
%      posterior (parietal + occipital) and global (128 channels).
%   5. IAF on the posterior region (Corcoran et al. 2018): Savitzky-Golay
%      smoothed spectra (11 bins, order 5), search window 7-13 Hz.
%      PAF = median of the channel peaks (>= 3 channels with a clear peak);
%      CoG = power-weighted mean frequency of the alpha peak.
%      Reference IAF for the bands = PAF in EC (fallbacks: CoG EC, PAF total,
%      CoG total, group median; the source is written in the tables).
%   6. Frequency limits: beta = up to 25 Hz, aperiodic fit = 2-30 Hz,
%      relative power = / power 2-30 Hz. Above ~20-25 Hz scalp EEG is
%      contaminated by muscle activity (Whitham et al. 2007,
%      https://doi.org/10.1016/j.clinph.2007.04.027); in this dataset a
%      residual 25-40 Hz bump was found in sub-28 (+4-8 dB) and a smaller
%      29-36 Hz bump in sub-12, 13, 19, 32 (version 1 results).
%   7. Aperiodic fit ('fixed' mode, 2-30 Hz): robust initial fit, iterative
%      Gaussian peak estimates (max 8 peaks, threshold 2 SD and 0.1 log10
%      units, width 1-8 Hz), aperiodic refit on the peak-removed spectrum.
%      This is the FOOOF / specparam algorithm (Donoghue et al. 2020)
%      without the final joint non-linear optimization of the peaks. On the
%      data of sub-01 it agrees with fooof 1.1.1 (exponent difference:
%      median -0.02, 95th percentile 0.05).
%
% OUTPUTS (BIDS/derivatives/spectral/)
%   sub-XX/sub-XX_task-rest_desc-spectra.mat   all spectra and per-channel values
%   sub-XX/sub-XX_task-rest_desc-spectral.png  spectra + topographies
%   table_iaf.csv          IAF per subject x condition (+ reference IAF, bands)
%   table_bandpower.csv    long: subject x condition x region x band set x band
%   table_ratios.csv       long: TBR, ATR, ABR (absolute and periodic)
%   table_aperiodic.csv    long: offset, exponent, peaks per region
%   table_berger.csv       alpha EC vs EO per region
%   table_roi_quality.csv  % interpolated / bridged channels in every region
%   table_wide.csv         ONE ROW PER SUBJECT with the main measures
%                          (ready for SPSS / jamovi / R) + quality flags
%   table_by_period.csv    long: subject x period (1-8) x region; band
%                          powers, aperiodic parameters, blinks per minute
%   table_roi_definitions.csv  channels of every region
%   tables_README.txt      meaning of every column
%   group_spectra.png      grand-average spectra and IAF distribution
%
% REQUIREMENTS
%   MATLAB R2026a, Signal Processing Toolbox (sgolayfilt, hann), EEGLAB
%   (pop_loadset, topoplot). Inputs: BIDS/derivatives/preproc-spectral.
%
% RUN TIME: about 1 min per subject.
% =========================================================================

clear; clc;

%% ========================================================================
%  0. CONFIGURATION
%  ========================================================================
cfg = struct();
base = 'C:\Users\chiar\Documents\MATLAB\FONDECYT\DATOS\Resting State';
cfg.in_root  = fullfile(base, 'BIDS', 'derivatives', 'preproc-spectral');
cfg.out_root = fullfile(base, 'BIDS', 'derivatives', 'spectral');
cfg.overwrite = false;
cfg.n_eeg = 128;

% Spectral estimation
cfg.win_s        = 4;          % window length (s) -> 0.25 Hz resolution
cfg.win_step_s   = 2;          % 50 % overlap
cfg.min_windows  = 20;         % fewer clean windows (80 s) in a condition -> flag
cfg.f_max_store  = 45;         % frequencies stored (Hz)

% Fixed bands (Hz)
cfg.fixed = struct('theta', [4 8], 'alpha', [8 13], 'beta', [13 25]);
% IAF-based bands (Klimesch 1999), offsets from the IAF (Hz)
cfg.iaf_bands = struct('theta', [-6 -4], 'alpha', [-4 2], 'alpha_low1', [-4 -2], ...
                       'alpha_low2', [-2 0], 'alpha_up', [0 2], 'beta', [2 NaN]);
cfg.beta_upper = 25;           % upper edge of the IAF-based beta band (muscle, see header)
cfg.rel_range  = [2 30];       % denominator of the relative power

% IAF (Corcoran et al. 2018)
cfg.iaf_window = [7 13];       % search window (Hz)
cfg.sg_frame   = 11;           % Savitzky-Golay frame (bins)
cfg.sg_order   = 5;
cfg.iaf_cmin   = 3;            % minimum channels with a peak
cfg.cog_limits = [7 14];       % the CoG window can never extend beyond these limits (Hz)

% Aperiodic fit (spectral parameterization)
cfg.fit_range   = [2 30];
cfg.max_peaks   = 8;
cfg.peak_thr_sd = 2;
cfg.min_peak_h  = 0.1;         % log10 units
cfg.peak_width  = [1 8];       % Hz (FWHM limits)

% Narrow-band (non-physiological) peaks
cfg.narrow_range_hz = [24 47];
cfg.narrow_db       = 3;
cfg.narrow_width_hz = 1.5;
cfg.narrow_halfwidth_interp = 1;   % +/- Hz interpolated

% Regions
cfg.cluster_deg = 20;          % channels within 20 deg of the 10-20 position
cfg.regions = struct('frontal', {{'Fz','F3','F4'}}, 'central', {{'Cz','C3','C4'}}, ...
                     'parietal', {{'Pz','P3','P4'}}, 'occipital', {{'Oz','O1','O2'}}, ...
                     'posterior', {{'Pz','P3','P4','Oz','O1','O2'}});
% Regions written to the wide table (plus 'global')
cfg.wide_rois = {'frontal','central','parietal','occipital','posterior','global','Cz','Fz'};

% Per-period values (each of the 8 periods)
cfg.period_rois = {'posterior','frontal','central','parietal','occipital','global','Cz','Fz'};

% Blink detection on the bipolar VEOG (VEOG_U - VEOG_D), 1-10 Hz band-pass
cfg.veog_labels    = {'VEOG_U', 'VEOG_D'};
cfg.blink_thr_uV   = 75;       % peak amplitude of a blink (|uV|)
cfg.blink_min_gap_s = 0.3;     % minimum distance between two blinks (s)

% Quality flags (written in table_wide.csv; they do NOT exclude anybody)
% Suspicious EC period (eyes probably open) = BOTH: posterior log10 alpha
% power not higher than the median of the EO periods + margin, AND VEOG
% events >= the rate below. Neither criterion alone is specific: VEOG
% events also occur with the eyes closed, and some subjects have weak
% alpha reactivity with the eyes really closed.
cfg.flag_ec_alpha_margin   = 0.1; % log10 units (x1.26)
cfg.flag_ec_veog_per_min   = 5;   % VEOG events per minute
cfg.flag_min_r2            = 0.90;% aperiodic fit R^2 (posterior or global) below this

cfg.save_figures = true;
cfg.version = '2.1 (2026-11-07)';

%% ========================================================================
%  1. SETUP
%  ========================================================================
if ~exist('pop_loadset', 'file')
    if exist('eeglab', 'file'), eeglab('nogui'); else, error('EEGLAB not found.'); end
end
assert(exist('sgolayfilt', 'file') == 2, 'Signal Processing Toolbox is required.');
if ~isfolder(cfg.in_root), error('Input folder not found:\n  %s', cfg.in_root); end
if ~isfolder(cfg.out_root), mkdir(cfg.out_root); end
save(fullfile(cfg.out_root, 'spectral_config.mat'), 'cfg');

files = dir(fullfile(cfg.in_root, 'sub-*', 'eeg', '*_desc-clean_eeg.set'));
[~, o] = sort({files.name});  files = files(o);
nS = numel(files);
fprintf('Found %d subject(s).\n', nS);
log_fid = fopen(fullfile(cfg.out_root, 'spectral_log.txt'), 'a');
fprintf(log_fid, '\n%s\nspectral_analysis_rs.m version %s - %s\n', repmat('=', 1, 70), cfg.version, ...
        char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm')));

%% ========================================================================
%  2. PER-SUBJECT SPECTRA (pass 1)
%  ========================================================================
for iS = 1:nS
    sub = extractBefore(files(iS).name, '_task');
    out_dir = fullfile(cfg.out_root, sub);
    out_mat = fullfile(out_dir, [sub '_task-rest_desc-spectra.mat']);
    fprintf('[%d/%d] %s\n', iS, nS, sub);
    if ~cfg.overwrite && isfile(out_mat)
        fprintf('  already done: skipped\n');
        continue
    end
    if ~isfolder(out_dir), mkdir(out_dir); end
    try
        EEG = pop_loadset('filename', files(iS).name, 'filepath', files(iS).folder);
        prov = load(fullfile(files(iS).folder, [sub '_task-rest_desc-provenance.mat']), 'epochs', 'channels');
        S = subject_spectra(EEG, prov, cfg);
        S.subject = sub;
        save(out_mat, '-struct', 'S');
        fprintf(log_fid, '%s: windows EO %d, EC %d; narrow peaks: %s\n', sub, S.n_windows.EO, ...
                S.n_windows.EC, S.narrow_peaks_txt);
    catch ME
        fprintf(2, '  ERROR: %s\n', ME.message);
        fprintf(log_fid, '%s: ERROR %s (%s line %d)\n', sub, ME.message, ME.stack(1).name, ME.stack(1).line);
    end
    clear EEG S prov
end

%% ========================================================================
%  3. MEASURES (pass 2: needs the group median IAF as last fallback)
%  ========================================================================
mats = dir(fullfile(cfg.out_root, 'sub-*', '*_desc-spectra.mat'));
[~, o] = sort({mats.name});  mats = mats(o);
all_S = cell(numel(mats), 1);
for i = 1:numel(mats)
    all_S{i} = load(fullfile(mats(i).folder, mats(i).name));
end

% IAF per subject and condition
iafT = table();
for i = 1:numel(all_S)
    for cond = {'EO', 'EC', 'total'}
        r = all_S{i}.iaf.(cond{1});
        iafT = [iafT; table(string(all_S{i}.subject), string(cond{1}), r.paf, r.n_peak_channels, ...
            r.cog, r.f1, r.f2, 'VariableNames', {'subject','condition','PAF_hz', ...
            'n_channels_with_peak','CoG_hz','cog_f1_hz','cog_f2_hz'})]; %#ok<AGROW>
    end
end
group_iaf = median(iafT.PAF_hz(iafT.condition == "EC"), 'omitnan');

bandT = table(); ratT = table(); apT = table(); bergT = table(); qT = table(); wide = table();
ref = table(); perT = table();
for i = 1:numel(all_S)
    S = all_S{i};
    [iaf_ref, src] = reference_iaf(S.iaf, group_iaf);
    ref = [ref; table(string(S.subject), iaf_ref, string(src), ...
           'VariableNames', {'subject','IAF_ref_hz','IAF_ref_source'})]; %#ok<AGROW>
    bands = make_bands(iaf_ref, cfg);
    [b, r, a, bg, q, w] = subject_measures(S, bands, iaf_ref, src, cfg);
    pt = period_measures(S, bands, cfg);
    w = add_flags(w, S, a, pt, cfg);
    bandT = [bandT; b]; ratT = [ratT; r]; apT = [apT; a]; bergT = [bergT; bg]; %#ok<AGROW>
    qT = [qT; q]; wide = [wide; w]; perT = [perT; pt]; %#ok<AGROW>
    if cfg.save_figures
        try
            subject_figure(S, bands, iaf_ref, cfg, fullfile(cfg.out_root, S.subject, ...
                           [S.subject '_task-rest_desc-spectral.png']));
        catch MEf
            fprintf(log_fid, '%s: figure not saved (%s)\n', S.subject, MEf.message);
        end
    end
end
iafT = join(iafT, ref, 'Keys', 'subject');       % add the reference IAF to every row

writetable(iafT,  fullfile(cfg.out_root, 'table_iaf.csv'));
writetable(bandT, fullfile(cfg.out_root, 'table_bandpower.csv'));
writetable(ratT,  fullfile(cfg.out_root, 'table_ratios.csv'));
writetable(apT,   fullfile(cfg.out_root, 'table_aperiodic.csv'));
writetable(bergT, fullfile(cfg.out_root, 'table_berger.csv'));
writetable(qT,    fullfile(cfg.out_root, 'table_roi_quality.csv'));
writetable(wide,  fullfile(cfg.out_root, 'table_wide.csv'));
writetable(perT,  fullfile(cfg.out_root, 'table_by_period.csv'));
write_roi_definitions(all_S{1}.rois, fullfile(cfg.out_root, 'table_roi_definitions.csv'));
write_dictionary(fullfile(cfg.out_root, 'tables_README.txt'), cfg);
if cfg.save_figures && ~isempty(all_S)
    try
        group_figure(all_S, iafT, cfg, fullfile(cfg.out_root, 'group_spectra.png'));
    catch MEf
        fprintf(log_fid, 'group figure not saved (%s)\n', MEf.message);
    end
end
fprintf(log_fid, 'Group median PAF (EC): %.2f Hz\n', group_iaf);
fclose(log_fid);
fprintf('\nDone. Outputs in:\n  %s\n', cfg.out_root);


%% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================

% -------------------------------------------------------------------------
function S = subject_spectra(EEG, prov, cfg)
% SUBJECT_SPECTRA  Welch-median spectra (EO, EC, total) of the 128 scalp
% channels, narrow-band peak removal, regions and IAF.
fs = EEG.srate;
L  = round(cfg.win_s * fs);
X  = double(EEG.data(1:cfg.n_eeg, :));
w  = hann(L);
U  = sum(w.^2);
f  = (0:floor(L/2))' * fs / L;
keepf = f <= cfg.f_max_store;

nP = 8;                                     % number of periods
S.psd_period = NaN(nP, cfg.n_eeg, sum(keepf), 'single');
S.n_windows_period = zeros(1, nP);
S.period_condition = strings(1, nP);
for cond = {'EO', 'EC'}
    [starts, per] = clean_windows(prov.epochs, cond{1}, cfg);
    n = numel(starts);
    Pw = zeros(cfg.n_eeg, sum(keepf), n);
    for k = 1:n
        s0  = round(starts(k) * fs) + 1;
        seg = X(:, s0 : s0 + L - 1);
        seg = seg - mean(seg, 2);
        F   = fft(seg .* w', [], 2);
        p   = abs(F(:, 1:floor(L/2) + 1)).^2 / (fs * U);
        p(:, 2:end-1) = 2 * p(:, 2:end-1);      % one-sided spectrum
        Pw(:, :, k) = p(:, keepf);
    end
    if n > 0
        S.psd.(cond{1}) = median(Pw, 3) / median_bias(n);
    else
        S.psd.(cond{1}) = NaN(cfg.n_eeg, sum(keepf));
    end
    S.n_windows.(cond{1}) = n;
    % same estimator, separately for each period of this condition
    for ci = unique(per(:))'
        k = per == ci;
        S.psd_period(ci, :, :) = single(median(Pw(:, :, k), 3) / median_bias(sum(k)));
        S.n_windows_period(ci) = sum(k);
        S.period_condition(ci) = string(cond{1});
    end
end
S.f = f(keepf);
S.psd.total = (S.psd.EO + S.psd.EC) / 2;    % equal weight to EO and EC

% Narrow-band peaks: detected on the mean total spectrum, interpolated
[pk_hz, pk_db] = narrow_peaks(S.psd.total, S.f, cfg);
S.narrow_peaks_hz = pk_hz;
S.narrow_peaks_db = pk_db;
S.narrow_peaks_txt = strjoin(arrayfun(@(x) sprintf('%.2f', x), pk_hz, 'UniformOutput', false), ' ');
for cond = {'EO', 'EC', 'total'}
    S.psd.(cond{1}) = interpolate_peaks(S.psd.(cond{1}), S.f, pk_hz, cfg.narrow_halfwidth_interp);
end
for ci = 1:nP
    if S.n_windows_period(ci) > 0
        S.psd_period(ci, :, :) = single(interpolate_peaks(double(squeeze(S.psd_period(ci, :, :))), ...
                                 S.f, pk_hz, cfg.narrow_halfwidth_interp));
    end
end

% Blinks per minute in every period (eye-closure check)
S.blinks_per_min = blink_rate(EEG, cfg);

% Regions and quality
S.chanlocs = EEG.chanlocs(1:cfg.n_eeg);
S.rois = make_rois(S.chanlocs, cfg);
% Channel handling table, re-ordered to match the EEG channel order
H = prov.channels;
[tf, loc] = ismember(string({S.chanlocs.labels}), string(H.channel));
assert(all(tf), 'Channel labels of the provenance table do not match the EEG file.');
H = H(loc, :);
H.handling = string(H.handling);
S.channel_handling = H;

% IAF
for cond = {'EO', 'EC', 'total'}
    S.iaf.(cond{1}) = compute_iaf(S.psd.(cond{1})(S.rois.posterior, :), S.f, cfg);
end
end

% -------------------------------------------------------------------------
function [starts, period] = clean_windows(ep, cond, cfg)
% CLEAN_WINDOWS  Start times of 4-s windows made of 3 clean, consecutive
% 2-s epochs (step 1 s) of the same period, 50 % overlap between windows.
% PERIOD is the period number (1-8) of every window.
starts = []; period = [];
E = ep(string(ep.condition) == string(cond), :);
nSub = round(cfg.win_s / 1) - 1;            % 2-s epochs inside a 4-s window: t, t+1, t+2
for ci = unique(E.condition_index)'
    e  = E(E.condition_index == ci, :);
    on = round(e.onset_s * 1000) / 1000;
    ok = ~logical(e.rejected);
    last = -Inf;
    for j = 1:numel(on)
        need = on(j) + (0:nSub - 1);
        [tf, loc] = ismember(round(need * 1000) / 1000, on);
        if all(tf) && all(ok(loc)) && on(j) - last >= cfg.win_step_s - 1e-6
            starts(end+1) = on(j); %#ok<AGROW>
            period(end+1) = ci;    %#ok<AGROW>
            last = on(j);
        end
    end
end
end

% -------------------------------------------------------------------------
function b = median_bias(n)
% MEDIAN_BIAS  Bias of the median of n chi-square (2 dof) periodogram
% values relative to their mean (as in scipy.signal.welch, average='median').
ii = 1:floor((n - 1) / 2);
b = 1 + sum(1 ./ (2 * ii + 1) - 1 ./ (2 * ii));
end

% -------------------------------------------------------------------------
function [hz, db] = narrow_peaks(P, f, cfg)
% NARROW_PEAKS  Narrow peaks of the mean spectrum above a 4-Hz median background
hz = []; db = [];
m  = 10 * log10(mean(P, 1, 'omitnan') + eps);
df = f(2) - f(1);
bg = movmedian(m, round(4 / df));
r  = m - bg;
sel = f' >= cfg.narrow_range_hz(1) & f' <= cfg.narrow_range_hz(2) & abs(f' - 50) > 2;
r(~sel) = 0;
[pk, loc, wd] = findpeaks(r, 'MinPeakHeight', cfg.narrow_db);
keep = wd * df <= cfg.narrow_width_hz;
hz = f(loc(keep))';
db = pk(keep);
end

% -------------------------------------------------------------------------
function P = interpolate_peaks(P, f, pk_hz, hw)
% INTERPOLATE_PEAKS  Linear interpolation (in log power) across +/- hw Hz
for h = pk_hz
    bad = abs(f - h) <= hw;
    i0 = find(bad, 1, 'first') - 1;
    i1 = find(bad, 1, 'last') + 1;
    if i0 < 1 || i1 > numel(f), continue; end
    lp = log10(P);
    for c = 1:size(P, 1)
        lp(c, bad) = interp1(f([i0 i1]), lp(c, [i0 i1]), f(bad));
    end
    P = 10 .^ lp;
end
end

% -------------------------------------------------------------------------
function R = make_rois(locs, cfg)
% MAKE_ROIS  Channel clusters around 10-20 equivalent positions computed on
% an ideal sphere (EEGLAB axes: +X nose, +Y left ear, +Z vertex).
v = @(th, ph) [sind(th) * cosd(ph), sind(th) * sind(ph), cosd(th)];
P = struct();
P.Fpz = v(90, 0);   P.Fz = v(45, 0);   P.Cz = [0 0 1];
P.Pz  = v(45, 180); P.Oz = v(90, 180);
P.C3  = v(45, 90);  P.C4 = v(45, -90);
P.F7  = v(90, 54);  P.F8 = v(90, -54); P.P7 = v(90, 126); P.P8 = v(90, -126);
P.O1  = v(90, 162); P.O2 = v(90, -162);
nz = @(a) a / norm(a);
P.F3 = nz(P.Fz + P.F7); P.F4 = nz(P.Fz + P.F8);
P.P3 = nz(P.Pz + P.P7); P.P4 = nz(P.Pz + P.P8);
xyz = [[locs.X]' [locs.Y]' [locs.Z]'];
u = xyz ./ vecnorm(xyz, 2, 2);
R = struct();
sites = {'Fz','F3','F4','Cz','C3','C4','Pz','P3','P4','Oz','O1','O2'};
for s = sites
    ang = acosd(max(-1, min(1, u * P.(s{1})')));
    R.(s{1}) = find(ang <= cfg.cluster_deg)';
    [~, R.nearest.(s{1})] = min(ang);
end
for g = fieldnames(cfg.regions)'
    members = cellfun(@(x) R.(x), cfg.regions.(g{1}), 'UniformOutput', false);
    R.(g{1}) = unique([members{:}]);
end
R.global = 1:numel(locs);
R.names = [sites, fieldnames(cfg.regions)', {'global'}];
R.labels = {locs.labels};
end

% -------------------------------------------------------------------------
function r = compute_iaf(P, f, cfg)
% COMPUTE_IAF  PAF and CoG on posterior channels (after Corcoran et al. 2018)
r = struct('paf', NaN, 'n_peak_channels', 0, 'cog', NaN, 'f1', NaN, 'f2', NaN);
if any(isnan(P(:))), return; end
win = f >= cfg.iaf_window(1) & f <= cfg.iaf_window(2);
iw  = find(win);
pafs = [];
for c = 1:size(P, 1)
    s = sgolayfilt(P(c, :), cfg.sg_order, cfg.sg_frame);
    [~, j] = max(s(iw));
    i = iw(j);
    if j > 1 && j < numel(iw) && s(i) > s(i - 1) && s(i) > s(i + 1)
        pafs(end+1) = f(i); %#ok<AGROW>
    end
end
r.n_peak_channels = numel(pafs);
if numel(pafs) >= cfg.iaf_cmin, r.paf = median(pafs); end
% CoG on the mean posterior spectrum, between the minima that bound the
% alpha peak (first-derivative sign changes), clamped to cfg.cog_limits
m = sgolayfilt(mean(P, 1), cfg.sg_order, cfg.sg_frame);
[~, j] = max(m(iw));
i = iw(j);
if j == 1 || j == numel(iw), return; end
d = diff(m);
lo = i; while lo > 1 && f(lo - 1) >= cfg.cog_limits(1) && d(lo - 1) > 0, lo = lo - 1; end
hi = i; while hi < numel(f) && f(hi + 1) <= cfg.cog_limits(2) && d(hi) < 0, hi = hi + 1; end
r.f1 = f(lo); r.f2 = f(hi);
r.cog = sum(f(lo:hi)' .* m(lo:hi)) / sum(m(lo:hi));
end

% -------------------------------------------------------------------------
function [iaf, src] = reference_iaf(I, group_iaf)
% REFERENCE_IAF  IAF used to define the individual bands
cand = {I.EC.paf, 'PAF_EC'; I.EC.cog, 'CoG_EC'; I.total.paf, 'PAF_total'; ...
        I.total.cog, 'CoG_total'; group_iaf, 'group_median'};
for k = 1:size(cand, 1)
    if isfinite(cand{k, 1}), iaf = cand{k, 1}; src = cand{k, 2}; return; end
end
iaf = NaN; src = 'none';
end

% -------------------------------------------------------------------------
function B = make_bands(iaf, cfg)
% MAKE_BANDS  Band edges: set 'iaf' (Klimesch 1999) and set 'fixed'
B = struct();
for b = fieldnames(cfg.iaf_bands)'
    e = iaf + cfg.iaf_bands.(b{1});
    if isnan(cfg.iaf_bands.(b{1})(2)), e(2) = cfg.beta_upper; end
    e(1) = max(e(1), cfg.fit_range(1));
    B.iaf.(b{1}) = e;
end
B.fixed = cfg.fixed;
end

% -------------------------------------------------------------------------
function M = fit_spectrum(p, f, cfg)
% FIT_SPECTRUM  Aperiodic + periodic parameterization of one spectrum
% (FOOOF / specparam algorithm, 'fixed' aperiodic mode, Donoghue et al. 2020):
%   1. robust aperiodic fit: OLS of log10 power on log10 frequency, then a
%      refit on the points at or below the first fit;
%   2. iterative Gaussian peak estimates on the flattened spectrum (peak
%      height > peak_thr_sd * SD and > min_peak_h; FWHM within peak_width);
%   3. aperiodic refit on the spectrum minus the Gaussian peaks.
s  = f >= cfg.fit_range(1) & f <= cfg.fit_range(2);
fx = f(s); x = log10(fx); y = log10(p(s)); y = y(:);
A  = [ones(numel(x), 1), -x(:)];
b  = A \ y;
fl = y - A * b;  fl(fl < 0) = 0;
msk = fl <= prctile_simple(fl, 0.025);
b  = A(msk, :) \ y(msk);
flat = y - A * b;  res = flat;  df = fx(2) - fx(1);
pk = zeros(0, 3);
for k = 1:cfg.max_peaks
    [h, i] = max(res);
    if h <= cfg.peak_thr_sd * std(res) || h <= cfg.min_peak_h, break; end
    l = i; while l > 1 && res(l) > h / 2, l = l - 1; end
    r = i; while r < numel(res) && res(r) > h / 2, r = r + 1; end
    hw = min(i - l, r - i) * df;
    sd = min(max(2 * hw / 2.3548, cfg.peak_width(1) / 2), cfg.peak_width(2) / 2);
    g  = h * exp(-(fx(:) - fx(i)).^2 / (2 * sd^2));
    res = res - g;
    pk(end+1, :) = [fx(i), h, 2 * sd]; %#ok<AGROW>
end
gsum = flat - res;
b = A \ (y - gsum);
M.offset = b(1);
M.exponent = b(2);
M.peaks = pk;                                  % [centre frequency, height (log10), width (Hz)]
ap = M.offset - M.exponent * log10(f(:)');     % aperiodic fit over all frequencies
M.ap_log = ap;
M.flat_log = log10(p(:)') - ap;                % periodic component (log10 units)
yhat = A * b + gsum;
M.r2 = 1 - sum((y - yhat).^2) / sum((y - mean(y)).^2);
end

% -------------------------------------------------------------------------
function q = prctile_simple(v, p)
v = sort(v(:)); n = numel(v);
q = v(max(1, min(n, round(p / 100 * (n - 1)) + 1)));
end

% -------------------------------------------------------------------------
function [bandT, ratT, apT, bergT, qT, W] = subject_measures(S, B, iaf_ref, src, cfg)
% SUBJECT_MEASURES  All measures of one subject (long tables + wide row)
sub = string(S.subject);
f = S.f(:)';
bandT = table(); ratT = table(); apT = table(); bergT = table(); qT = table();
W = table(sub, iaf_ref, string(src), S.iaf.EC.paf, S.iaf.EC.cog, S.iaf.EO.paf, S.iaf.EO.cog, ...
    S.n_windows.EO, S.n_windows.EC, string(S.narrow_peaks_txt), ...
    'VariableNames', {'subject','IAF_ref_hz','IAF_ref_source','PAF_EC','CoG_EC','PAF_EO','CoG_EO', ...
    'n_windows_EO','n_windows_EC','narrow_peaks_removed_hz'});
rel = f >= cfg.rel_range(1) & f <= cfg.rel_range(2);
H = S.channel_handling;

for roi = S.rois.names
    idx = S.rois.(roi{1});
    % quality of the region
    h = H.handling(idx);
    qT = [qT; table(sub, string(roi{1}), numel(idx), 100 * mean(contains(h, "interpolated")), ...
          100 * mean(h == "bridged_copy"), 100 * mean(H.valid_for_connectivity(idx)), ...
          'VariableNames', {'subject','roi','n_channels','pct_interpolated', ...
          'pct_bridged_copy','pct_independent'})]; %#ok<AGROW>
    alpha_store = struct();
    for cond = {'EO', 'EC', 'total'}
        p = mean(S.psd.(cond{1})(idx, :), 1);
        M = fit_spectrum(p, S.f, cfg);
        tot = trapz(f(rel), p(rel));
        % aperiodic and peaks
        [acf, ah, abw] = pick_peak(M.peaks, B.iaf.alpha);
        [tcf, th] = pick_peak(M.peaks, B.iaf.theta);
        [bcf, bh] = pick_peak(M.peaks, B.iaf.beta);
        apT = [apT; table(sub, string(cond{1}), string(roi{1}), M.offset, M.exponent, M.r2, ...
               size(M.peaks, 1), acf, ah, abw, tcf, th, bcf, bh, ...
               'VariableNames', {'subject','condition','roi','ap_offset','ap_exponent','fit_r2', ...
               'n_peaks','alpha_peak_hz','alpha_peak_height','alpha_peak_width_hz', ...
               'theta_peak_hz','theta_peak_height','beta_peak_hz','beta_peak_height'})]; %#ok<AGROW>
        for bset = {'iaf', 'fixed'}
            BB = B.(bset{1});
            vals = struct();
            for b = fieldnames(BB)'
                e = BB.(b{1});
                sel = f >= e(1) & f <= e(2);
                ab = trapz(f(sel), p(sel));
                per_log = mean(M.flat_log(sel));
                per_abs = trapz(f(sel), max(p(sel) - 10 .^ M.ap_log(sel), 0));
                vals.(b{1}) = [ab, per_log];
                bandT = [bandT; table(sub, string(cond{1}), string(roi{1}), string(bset{1}), ...
                    string(b{1}), e(1), e(2), ab, ab / tot, log10(ab), per_log, per_abs, ...
                    'VariableNames', {'subject','condition','roi','band_set','band','f_low_hz', ...
                    'f_high_hz','abs_power_uV2','rel_power','log10_abs_power', ...
                    'periodic_log10','periodic_abs_uV2'})]; %#ok<AGROW>
            end
            % ratios: absolute (classic) and periodic (log10 difference)
            prs = {'TBR', 'theta', 'beta'; 'ATR', 'alpha', 'theta'; 'ABR', 'alpha', 'beta'};
            for k = 1:3
                num = vals.(prs{k, 2}); den = vals.(prs{k, 3});
                ratT = [ratT; table(sub, string(cond{1}), string(roi{1}), string(bset{1}), ...
                    string(prs{k, 1}), num(1) / den(1), log10(num(1) / den(1)), num(2) - den(2), ...
                    'VariableNames', {'subject','condition','roi','band_set','ratio', ...
                    'ratio_abs','log10_ratio_abs','log10_ratio_periodic'})]; %#ok<AGROW>
            end
            alpha_store.(bset{1}).(cond{1}) = vals.alpha;
            % wide table: main measures, IAF bands only
            if strcmp(bset{1}, 'iaf') && any(strcmp(roi{1}, cfg.wide_rois))
                pre = sprintf('%s_%s_', cond{1}, roi{1});
                W.([pre 'theta_abs'])  = vals.theta(1);
                W.([pre 'alpha_abs'])  = vals.alpha(1);
                W.([pre 'beta_abs'])   = vals.beta(1);
                W.([pre 'theta_per'])  = vals.theta(2);
                W.([pre 'alpha_per'])  = vals.alpha(2);
                W.([pre 'beta_per'])   = vals.beta(2);
                W.([pre 'TBR'])        = vals.theta(1) / vals.beta(1);
                W.([pre 'ATR'])        = vals.alpha(1) / vals.theta(1);
                W.([pre 'ABR'])        = vals.alpha(1) / vals.beta(1);
                W.([pre 'TBR_per'])    = vals.theta(2) - vals.beta(2);
                W.([pre 'ATR_per'])    = vals.alpha(2) - vals.theta(2);
                W.([pre 'ABR_per'])    = vals.alpha(2) - vals.beta(2);
                W.([pre 'ap_exponent']) = M.exponent;
                W.([pre 'ap_offset'])   = M.offset;
            end
        end
    end
    for bset = {'iaf', 'fixed'}
        a = alpha_store.(bset{1});
        bergT = [bergT; table(sub, string(roi{1}), string(bset{1}), a.EO(1), a.EC(1), ...
            a.EC(1) / a.EO(1), log10(a.EC(1) / a.EO(1)), a.EC(2) - a.EO(2), ...
            'VariableNames', {'subject','roi','band_set','alpha_EO_uV2','alpha_EC_uV2', ...
            'berger_ratio_EC_EO','log10_berger','berger_periodic_log10'})]; %#ok<AGROW>
        if strcmp(bset{1}, 'iaf') && any(strcmp(roi{1}, cfg.wide_rois))
            W.(sprintf('berger_%s', roi{1})) = a.EC(1) / a.EO(1);
        end
    end
    if any(strcmp(roi{1}, cfg.wide_rois))
        W.(sprintf('pct_interp_%s', roi{1})) = qT.pct_interpolated(end);
    end
end
end

% -------------------------------------------------------------------------
function bpm = blink_rate(EEG, cfg)
% BLINK_RATE  Blinks per minute in each of the 8 periods, from the bipolar
% VEOG (VEOG_U - VEOG_D) band-passed at 1-10 Hz. A blink = peak of |VEOG|
% above cfg.blink_thr_uV, at least cfg.blink_min_gap_s from the previous one.
% Counted from 4 s after the instruction (as the epochs) to the end of the
% period. With the eyes really closed the rate is ~0; a high rate in the
% EC periods means that the eyes were (partly) open or the lids fluttered.
bpm = NaN(1, 8);
lab = {EEG.chanlocs.labels};
iu = find(strcmp(lab, cfg.veog_labels{1}), 1);
id = find(strcmp(lab, cfg.veog_labels{2}), 1);
if isempty(iu) || isempty(id), return; end
fs = EEG.srate;
v  = double(EEG.data(iu, :) - EEG.data(id, :));
[b, a] = butter(2, [1 10] / (fs / 2), 'bandpass');
v  = filtfilt(b, a, v);
[~, pk] = findpeaks(abs(v), 'MinPeakHeight', cfg.blink_thr_uV, ...
                    'MinPeakDistance', round(cfg.blink_min_gap_s * fs));
ev = EEG.event;
for k = 1:numel(ev)
    if ~any(strcmp(ev(k).type, {'EO', 'EC'})) || ~isfield(ev, 'condition_index'), continue; end
    ci = ev(k).condition_index;
    if isempty(ci) || ci < 1 || ci > 8, continue; end
    s0 = ev(k).latency + 4 * fs;                 % skip the first 4 s (as the epochs)
    s1 = ev(k).latency + ev(k).duration;
    if s1 <= s0, continue; end
    bpm(ci) = sum(pk >= s0 & pk < s1) / ((s1 - s0) / fs / 60);
end
end

% -------------------------------------------------------------------------
function T = period_measures(S, B, cfg)
% PERIOD_MEASURES  Band powers (IAF bands) and aperiodic parameters of every
% period (1-8) in the regions cfg.period_rois. Used to check stability over
% time (vigilance), eye closure (blinks) and for the case studies.
T = table();
f = S.f(:)';
for ci = 1:size(S.psd_period, 1)
    if S.n_windows_period(ci) == 0, continue; end
    P = double(squeeze(S.psd_period(ci, :, :)));
    for roi = cfg.period_rois
        p = mean(P(S.rois.(roi{1}), :), 1);
        M = fit_spectrum(p, S.f, cfg);
        v = struct();
        for b = {'theta', 'alpha', 'beta'}
            sel = f >= B.iaf.(b{1})(1) & f <= B.iaf.(b{1})(2);
            v.(b{1}) = [trapz(f(sel), p(sel)), mean(M.flat_log(sel))];
        end
        T = [T; table(string(S.subject), ci, S.period_condition(ci), S.n_windows_period(ci), ...
            S.blinks_per_min(ci), string(roi{1}), v.theta(1), v.alpha(1), v.beta(1), ...
            v.theta(2), v.alpha(2), v.beta(2), M.exponent, M.offset, M.r2, ...
            'VariableNames', {'subject','period','condition','n_windows','blinks_per_min','roi', ...
            'theta_abs_uV2','alpha_abs_uV2','beta_abs_uV2','theta_periodic_log10', ...
            'alpha_periodic_log10','beta_periodic_log10','ap_exponent','ap_offset','fit_r2'})]; %#ok<AGROW>
    end
end
end

% -------------------------------------------------------------------------
function W = add_flags(W, S, apT, perT, cfg)
% ADD_FLAGS  Quality information and flags (they never exclude anybody)
isEC = S.period_condition == "EC";
isEO = S.period_condition == "EO";
W.blinks_per_min_EO = mean(S.blinks_per_min(isEO), 'omitnan');
W.blinks_per_min_EC = mean(S.blinks_per_min(isEC), 'omitnan');
r2 = apT.fit_r2(ismember(apT.roi, ["posterior", "global"]));
W.min_fit_r2 = min(r2);
fl = strings(0);
% Eyes-closed check per EC period (see cfg.flag_ec_alpha_margin)
pp = perT(perT.roi == "posterior", :);
la = log10(pp.alpha_abs_uV2);
ref = median(la(pp.condition == "EO"), 'omitnan') + cfg.flag_ec_alpha_margin;
susp = pp.condition == "EC" & la <= ref & pp.blinks_per_min >= cfg.flag_ec_veog_per_min;
% (build the text first: a table cell must receive exactly one value)
txt = "";
if any(susp), txt = string(strjoin(compose('EC%d', pp.period(susp)'), ' ')); end
W.suspicious_EC_periods = txt;
if any(susp), fl(end+1) = "EC_eyes_maybe_open"; end
if S.n_windows.EO < cfg.min_windows, fl(end+1) = "few_windows_EO"; end
if S.n_windows.EC < cfg.min_windows, fl(end+1) = "few_windows_EC"; end
if W.min_fit_r2 < cfg.flag_min_r2, fl(end+1) = "poor_aperiodic_fit"; end
if W.IAF_ref_source ~= "PAF_EC", fl(end+1) = "IAF_fallback"; end
if ~isempty(S.narrow_peaks_hz), fl(end+1) = "narrow_peak_removed"; end
ftxt = "";
if ~isempty(fl), ftxt = string(strjoin(fl, ' ')); end
W.flags = ftxt;
end

% -------------------------------------------------------------------------
function [cf, h, w] = pick_peak(pk, band)
% PICK_PEAK  Highest parameterized peak inside a band (NaN if none)
cf = NaN; h = NaN; w = NaN;
if isempty(pk), return; end
in = pk(:, 1) >= band(1) & pk(:, 1) <= band(2);
if ~any(in), return; end
c = find(in);
[~, j] = max(pk(c, 2));
cf = pk(c(j), 1); h = pk(c(j), 2); w = pk(c(j), 3);
end

% -------------------------------------------------------------------------
function subject_figure(S, B, iaf_ref, cfg, fn)
% SUBJECT_FIGURE  Posterior and Cz spectra (EO / EC with 1/f fits, IAF
% bands in different colours) and topographies of alpha (EO, EC), Berger
% and TBR (EO). Every map has its own colour scale: absolute maps use the
% 2nd-98th percentile of that map, difference / ratio maps are centred on 0.
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [20 20 1500 800]);
try, theme(fig, 'light'); catch, end
f = S.f(:)';
cols  = struct('EO', [0.1 0.3 0.9], 'EC', [0.85 0.15 0.15]);
bcol  = struct('theta', [0.55 0.80 0.55], 'alpha', [0.95 0.80 0.45], 'beta', [0.65 0.75 0.95]);
inFit = f >= cfg.fit_range(1) & f <= cfg.fit_range(2);
rr = {'posterior', 'Cz'};
for k = 1:2
    subplot(2, 4, [1 2] + 4 * (k - 1)); hold on      % left half: spectra
    for cond = {'EO', 'EC'}
        p = mean(S.psd.(cond{1})(S.rois.(rr{k}), :), 1);
        M = fit_spectrum(p, S.f, cfg);
        plot(f, log10(p), 'Color', cols.(cond{1}), 'LineWidth', 1.6, 'DisplayName', cond{1});
        plot(f(inFit), M.ap_log(inFit), '--', 'Color', cols.(cond{1}), 'LineWidth', 1.1, ...
             'HandleVisibility', 'off');
    end
    yl = ylim;
    for bnd = {'theta', 'alpha', 'beta'}
        e = B.iaf.(bnd{1});
        patch([e fliplr(e)], [yl(1) yl(1) yl(2) yl(2)], bcol.(bnd{1}), 'FaceAlpha', 0.25, ...
              'EdgeColor', 'none', 'DisplayName', bnd{1});
    end
    % frequencies outside the analysed range (above the aperiodic fit) in grey
    patch([cfg.fit_range(2) 40 40 cfg.fit_range(2)], [yl(1) yl(1) yl(2) yl(2)], [0.85 0.85 0.85], ...
          'FaceAlpha', 0.35, 'EdgeColor', 'none', 'DisplayName', 'not analysed');
    ylim(yl);
    xline(iaf_ref, 'k:', sprintf('IAF %.2f', iaf_ref), 'HandleVisibility', 'off');
    xlim([1 40]); xlabel('Hz'); ylabel('log_{10} power (\muV^2/Hz)'); grid on
    legend('Location', 'northeast');
    title(sprintf('%s - %s region (dashed = 1/f fit)', S.subject, rr{k}), 'Interpreter', 'none');
end
% topographies (per channel, IAF bands)
chan = S.chanlocs;
alpha = @(P) log10(band_power(P, f, B.iaf.alpha));
maps = {alpha(S.psd.EO), 'log_{10} alpha EO (\muV^2)', 'abs'; ...
        alpha(S.psd.EC), 'log_{10} alpha EC (\muV^2)', 'abs'; ...
        alpha(S.psd.EC) - alpha(S.psd.EO), 'Berger log_{10}(EC/EO)', 'centred'; ...
        log10(band_power(S.psd.EO, f, B.iaf.theta) ./ band_power(S.psd.EO, f, B.iaf.beta)), ...
        'log_{10} TBR EO', 'centred'};
pos = [3 4 7 8];
for k = 1:4
    subplot(2, 4, pos(k));
    v = maps{k, 1};
    if strcmp(maps{k, 3}, 'abs')
        lim = prctile(v, [2 98]);
        cm = parula(256);
    else
        lim = [-1 1] * max(abs(prctile(v, [2 98])));
        cm = blue_white_red(256);
    end
    if diff(lim) <= 0, lim = lim + [-1 1] * 1e-3; end
    topoplot(v, chan, 'electrodes', 'off', 'maplimits', lim);
    colormap(gca, cm); colorbar; title(maps{k, 2});
end
exportgraphics(fig, fn, 'Resolution', 100, 'BackgroundColor', 'white');
close(fig);
end

% -------------------------------------------------------------------------
function cm = blue_white_red(n)
% BLUE_WHITE_RED  Diverging colour map (blue = negative, white = 0, red = positive)
h = floor(n / 2);
up = linspace(0, 1, h)';
cm = [[up up ones(h, 1)]; [ones(n - h, 1) flipud(linspace(0, 1, n - h)') flipud(linspace(0, 1, n - h)')]];
end

% -------------------------------------------------------------------------
function bp = band_power(P, f, e)
sel = f >= e(1) & f <= e(2);
bp = trapz(f(sel), P(:, sel), 2);
end

% -------------------------------------------------------------------------
function group_figure(all_S, iafT, cfg, fn)
% GROUP_FIGURE  Grand-average posterior spectra (EO / EC) and IAF histogram
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [20 20 1300 450]);
try, theme(fig, 'light'); catch, end
f = all_S{1}.f(:)';
subplot(1, 2, 1); hold on
cols = struct('EO', [0.1 0.3 0.9], 'EC', [0.85 0.15 0.15]);
for cond = {'EO', 'EC'}
    Y = cell2mat(cellfun(@(S) log10(mean(S.psd.(cond{1})(S.rois.posterior, :), 1)), all_S, ...
                 'UniformOutput', false));
    m = mean(Y, 1, 'omitnan'); se = std(Y, 0, 1, 'omitnan') / sqrt(size(Y, 1));
    fill([f fliplr(f)], [m - se fliplr(m + se)], cols.(cond{1}), 'FaceAlpha', 0.2, ...
         'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(f, m, 'Color', cols.(cond{1}), 'LineWidth', 2, 'DisplayName', cond{1});
end
xlim([1 40]); grid on; legend; xlabel('Hz'); ylabel('log_{10} power (\muV^2/Hz)');
title(sprintf('Grand average, posterior region (n = %d, mean \\pm SEM)', numel(all_S)));
subplot(1, 2, 2);
v = iafT.PAF_hz(iafT.condition == "EC");
histogram(v, 7:0.5:14); xlabel('PAF eyes closed (Hz)'); ylabel('subjects'); grid on
title(sprintf('Individual alpha frequency: median %.2f Hz', median(v, 'omitnan')));
exportgraphics(fig, fn, 'Resolution', 110, 'BackgroundColor', 'white');
close(fig);
end

% -------------------------------------------------------------------------
function write_roi_definitions(R, fn)
rows = {};
for n = R.names
    rows(end+1, :) = {n{1}, numel(R.(n{1})), strjoin(R.labels(R.(n{1})), ' ')}; %#ok<AGROW>
end
T = cell2table(rows, 'VariableNames', {'roi','n_channels','channels'});
writetable(T, fn);
end

% -------------------------------------------------------------------------
function write_dictionary(fn, cfg)
txt = {
'SPECTRAL TABLES - spectral_analysis_rs.m'
'========================================='
'Conditions: EO (4 eyes-open periods), EC (4 eyes-closed periods), total (mean of EO and EC spectra).'
'Regions (roi): clusters of ~5 channels around Fz F3 F4 Cz C3 C4 Pz P3 P4 Oz O1 O2;'
'  frontal/central/parietal/occipital = 3 clusters each; posterior = parietal + occipital; global = 128 channels.'
'  Channels of every region: table_roi_definitions.csv.'
'Band sets: iaf = individual bands (Klimesch 1999): theta IAF-6..IAF-4, alpha IAF-4..IAF+2'
'  (alpha_low1 IAF-4..IAF-2, alpha_low2 IAF-2..IAF, alpha_up IAF..IAF+2), beta IAF+2..25 Hz;'
'  fixed = theta 4-8, alpha 8-13, beta 13-25 Hz (upper limit 25 Hz because of muscle activity).'
''
'table_iaf.csv: PAF_hz (peak), n_channels_with_peak, CoG_hz (centre of gravity) and its window;'
'  IAF_ref_hz = IAF used for the bands, IAF_ref_source = PAF_EC | CoG_EC | PAF_total | CoG_total | group_median.'
'table_bandpower.csv: abs_power_uV2 (area under the spectrum), rel_power (/ power 2-30 Hz),'
'  log10_abs_power, periodic_log10 (mean height above the 1/f fit, log10 units: 0.1 = +1 dB),'
'  periodic_abs_uV2 (area above the 1/f fit).'
'table_ratios.csv: TBR = theta/beta, ATR = alpha/theta, ABR = alpha/beta;'
'  ratio_abs (classic), log10_ratio_abs, log10_ratio_periodic (difference of periodic_log10:'
'  ratio of the oscillatory parts, independent of the 1/f slope).'
'table_aperiodic.csv: ap_offset, ap_exponent (slope of the 1/f component, fit 2-30 Hz, fixed mode),'
'  fit_r2, n_peaks, alpha/theta/beta peak: centre frequency, height above 1/f (log10), width (Hz);'
'  NaN = no peak above the 1/f background in that band.'
'table_berger.csv: alpha EO, EC, EC/EO ratio, log10 ratio, periodic difference.'
'table_roi_quality.csv: % interpolated channels, % bridged copies, % independent channels per region.'
'table_wide.csv: one row per subject; columns <condition>_<roi>_<measure> with the IAF bands:'
'  theta/alpha/beta _abs (uV^2) and _per (periodic log10); TBR/ATR/ABR (absolute ratios);'
'  TBR_per/ATR_per/ABR_per (periodic log10 ratios); ap_exponent, ap_offset;'
'  berger_<roi> (alpha EC/EO); pct_interp_<roi> (% interpolated channels);'
'  blinks_per_min_EO / _EC (mean over the 4 periods), min_fit_r2 (worst aperiodic fit, posterior/global),'
'  suspicious_EC_periods: EC periods with BOTH posterior log10 alpha <= EO median + 0.1 AND'
'  >= 5 VEOG events/min (eyes probably open in that period);'
'  flags: EC_eyes_maybe_open (>= 1 suspicious EC period), few_windows_EO/EC (< 20 clean 4-s windows),'
'  poor_aperiodic_fit (R^2 < 0.90), IAF_fallback, narrow_peak_removed. Flags never exclude anybody.'
'table_by_period.csv: one row per subject x period (1-8) x region: condition, n_windows,'
'  blinks_per_min (VEOG peaks > 75 uV, from 4 s after the instruction; blinks AND lid / eye movements,'
'  which also occur with the eyes closed), theta/alpha/beta absolute'
'  power (IAF bands) and periodic power (log10), ap_exponent, ap_offset, fit_r2.'
};
fid = fopen(fn, 'w'); fprintf(fid, '%s\n', txt{:}); fclose(fid);
end
