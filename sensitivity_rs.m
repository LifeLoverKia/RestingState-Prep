% =========================================================================
% sensitivity_rs.m   (version 2, 2026-11-07)
% =========================================================================
%
% PURPOSE
%   General sensitivity analysis of the spectral results (output of
%   spectral_analysis_rs.m). The whole analysis is run more than once, for
%   the subject sets listed in cfg.runs (default: all subjects, and all
%   subjects without sub-04), and the runs are then compared.
%
%   For every run it answers six questions:
%
%   1. SCREENING - Which subjects have extreme values, and in which
%      measures? Robust z-scores (median / MAD) of every main measure, plus
%      data-driven artifact indices computed from the spectra. Artifact
%      indices are ONE-SIDED (only an excess counts as an artifact):
%        - ocular residual  : frontal 1-3 Hz power, EO minus EC (eye
%                             movements only exist with the eyes open)
%        - frontal slow     : frontal minus posterior 1-3 Hz power (EO)
%        - muscle residual  : power excess at 29-36 Hz above a straight
%                             line (log-log) from 22-26 to 38-42 Hz
%        - eye closure      : suspicious EC periods = alpha not higher than
%                             in the EO periods AND many VEOG events (blinks
%                             or lid movements). Each criterion alone is not
%                             specific, so both are required.
%
%   2. CONDITION EFFECT AND INFLUENCE - Does any single subject (or pair of
%      subjects) drive a result? For the within-subject difference EC - EO
%      of every measure (for alpha this is the Berger effect) the script
%      computes the effect on all subjects and then removes EVERY subject
%      (leave-one-out) and EVERY PAIR of subjects (leave-two-out), recording
%      how much the effect size and its significance change (jackknife
%      influence; Efron & Tibshirani 1993). p-values are also obtained by
%      sign-flip permutation and the effect size gets a bootstrap 95 % CI.
%      It also gives the correlation between EO and EC across subjects
%      (r_EO_EC): a high value means that the ranking of the subjects is
%      the same in both conditions (trait-like measure).
%
%   3. ASSOCIATIONS WITH SUBJECT VARIABLES - Pearson and Spearman
%      correlations of every measure, in every condition (EO, EC, total)
%      and of its reactivity (EC - EO, e.g. the Berger effect for alpha),
%      with every variable of participants.tsv listed in cfg.covariates
%      (demographic / cognitive) and cfg.nuisance (technical: they should
%      NOT correlate). Each covariate gets FDR-corrected q-values
%      (Benjamini & Hochberg 1995) and leave-one-out / leave-two-out
%      ranges, so that a correlation driven by one or two subjects is
%      visible.
%
%   4. SUBJECT MULTIVERSE - Every EC - EO effect is recomputed for
%      pre-specified exclusion sets (cfg.exclusion_sets).
%
%   5. ANALYSIS MULTIVERSE (specification curve) - EC - EO difference of
%      alpha, theta, beta, TBR, ATR and ABR for every combination of band
%      set (IAF / fixed) x measure (log absolute, relative, periodic) x
%      region x exclusion set (Steegen et al. 2016; Simonsohn et al. 2020).
%
%   6. CASE STUDIES - One figure + one text summary for each subject in
%      cfg.case_subjects (if it is present in that run).
%
%   Finally, the runs are compared (effects and correlations side by side).
%
% MISSING SUBJECTS
%   The subjects analysed are those of table_wide.csv that also have their
%   spectra file in the spectral folder, minus the exclusions of the run.
%   A subject that appears elsewhere (participants.tsv, exclusion sets,
%   case list) but not in the spectral folder is simply skipped.
%   To remove a subject you do NOT need to move folders: add it to
%   cfg.runs(k).exclude.
%
% OUTPUTS (BIDS/derivatives/spectral/sensitivity/<run name>/)
%   table_screening.csv          long: subject x measure x condition, value, robust z
%   table_artifact_indices.csv   artifact / eye-closure indices and their robust z
%   table_subject_summary.csv    one row per subject: extreme values, artifacts, influence
%   table_effects.csv            EC - EO effects (dz, t, p, sign-flip p, bootstrap CI,
%                                LOO and L2O ranges, r_EO_EC)
%   table_influence.csv          long: subject x effect, change of dz when removed
%   table_covariates.csv         long: measure x condition x covariate: r, p, q (FDR),
%                                rho, LOO and L2O ranges
%   table_exclusion_sets.csv     every effect for every exclusion set
%   table_specification_curve.csv
%   fig_*.png, cases/sub-XX_case.png / .txt, sensitivity_README.txt
% and in BIDS/derivatives/spectral/sensitivity/:
%   table_compare_runs_effects.csv, table_compare_runs_covariates.csv,
%   fig_compare_runs.png
%
% REQUIREMENTS
%   MATLAB R2026a, base MATLAB only (no toolbox needed).
%
% HOW TO ADAPT IT (for the papers)
%   - Covariates: any column of participants.tsv (add your cognitive and
%     questionnaire scores there, one column per score).
%   - Reference effect: change the function effect_matrix.
%   - Runs / exclusion sets: edit cfg.runs and cfg.exclusion_sets.
%
% REFERENCES
%   Benjamini & Hochberg (1995) J R Stat Soc B 57:289-300. doi:10.1111/j.2517-6161.1995.tb02031.x
%   Efron & Tibshirani (1993) An Introduction to the Bootstrap. Chapman & Hall.
%   Iglewicz & Hoaglin (1993) How to detect and handle outliers. ASQC.
%   Leys et al. (2013) J Exp Soc Psychol 49:764-766. doi:10.1016/j.jesp.2013.03.013
%   Steegen et al. (2016) Perspect Psychol Sci 11:702-712. doi:10.1177/1745691616658637
%   Simonsohn, Simmons & Nelson (2020) Nat Hum Behav 4:1208-1214. doi:10.1038/s41562-020-0912-z
%   Clayson et al. (2021) NeuroImage 245:118712. doi:10.1016/j.neuroimage.2021.118712
% =========================================================================

clear; clc;

%% ========================================================================
%  0. CONFIGURATION
%  ========================================================================
cfg = struct();
base = 'C:\Users\chiar\Documents\MATLAB\FONDECYT\DATOS\Resting State';
cfg.spec_root = fullfile(base, 'BIDS', 'derivatives', 'spectral');
cfg.out_root  = fullfile(cfg.spec_root, 'sensitivity');
cfg.participants_file = fullfile(base, 'BIDS', 'participants.tsv');

% Runs: the full analysis is repeated for each subject set, then compared
cfg.runs = struct('name',    {'all_subjects', 'without_sub-04'}, ...
                  'exclude', {{},             {'sub-04'}});

% Subject variables (columns of participants.tsv)
cfg.covariates = {'age', 'sex', 'handedness'};            % of interest
cfg.nuisance   = {'sequence', 'original_srate_hz'};       % technical: should not correlate

% Measures and regions of table_wide.csv
cfg.rois     = {'posterior','frontal','central','Cz','Fz','global'};
cfg.measures = {'theta_abs','alpha_abs','beta_abs','theta_per','alpha_per','beta_per', ...
                'TBR','ATR','ABR','TBR_per','ATR_per','ABR_per','ap_exponent','ap_offset'};
% Measures that are log10-transformed before any statistic (skewed, positive)
cfg.log_measures = {'theta_abs','alpha_abs','beta_abs','TBR','ATR','ABR'};

% Statistics
cfg.z_thr  = 3.5;                    % robust z above this = extreme (Iglewicz & Hoaglin 1993)
cfg.alpha  = 0.05;
cfg.n_perm = 10000;                  % sign-flip permutations
cfg.n_boot = 5000;                   % bootstrap resamples
cfg.seed   = 20261107;               % random seed (results are reproducible)

% Eye-closure check. Suspicious EC period (eyes probably open) = BOTH:
% posterior log10 alpha power <= median of the EO periods + margin, AND
% VEOG events >= rate. Neither criterion alone is specific (VEOG events
% also occur with the eyes closed; some subjects have weak alpha reactivity).
cfg.ec_alpha_margin  = 0.1;          % log10 units (x1.26)
cfg.blink_abs_thr    = 5;            % VEOG events per minute

% Pre-specified exclusion sets (subject multiverse inside each run)
cfg.exclusion_sets = struct( ...
    'name',     {'none', 'EC_eyes_open', 'muscle', 'ocular_residual', 'interpolated', 'bridges'}, ...
    'subjects', {{},     {'sub-04'},     {'sub-28'}, {'sub-07','sub-11','sub-40'}, {'sub-12'}, ...
                 {'sub-05','sub-19','sub-22'}}, ...
    'reason',   {'all subjects of the run', '2 suspicious EC periods (eyes probably open)', ...
                 'residual muscle activity 25-40 Hz', 'frontal 1-3 Hz excess with eyes open', ...
                 '39 bad channels (block B)', 'largest bridged-channel groups'});
cfg.add_all_listed_set = true;

% Specification curve
cfg.spec_constructs = {'alpha', 'theta', 'beta', 'TBR', 'ATR', 'ABR'};
cfg.spec_rois = {'Fz','F3','F4','Cz','C3','C4','Pz','P3','P4','Oz','O1','O2', ...
                 'frontal','central','parietal','occipital','posterior','global'};

% Case studies
cfg.case_subjects = {'sub-04', 'sub-15', 'sub-28', 'sub-40', 'sub-07'};
cfg.case_roi = 'posterior';

%% ========================================================================
%  1. LOAD (once for all runs)
%  ========================================================================
if ~isfolder(cfg.spec_root), error('Spectral folder not found:\n  %s', cfg.spec_root); end
if ~isfolder(cfg.out_root), mkdir(cfg.out_root); end
save(fullfile(cfg.out_root, 'sensitivity_config.mat'), 'cfg');

W  = readtable(fullfile(cfg.spec_root, 'table_wide.csv'), 'TextType', 'string');
BP = readtable(fullfile(cfg.spec_root, 'table_bandpower.csv'), 'TextType', 'string');
RT = readtable(fullfile(cfg.spec_root, 'table_ratios.csv'), 'TextType', 'string');
PT = readtable(fullfile(cfg.spec_root, 'table_by_period.csv'), 'TextType', 'string');
if ~ismember('flags', W.Properties.VariableNames), W.flags = strings(height(W), 1); end
if ~isstring(W.flags), W.flags = string(W.flags); end
W.flags = fillmissing(W.flags, 'constant', "");

% keep only the subjects whose spectra file exists
has = false(height(W), 1);
for i = 1:height(W)
    has(i) = isfile(spectra_file(cfg, W.subject(i)));
end
if any(~has)
    fprintf('Spectra file not found (skipped): %s\n', strjoin(W.subject(~has)', ', '));
end
W = W(has, :);

% subject variables
Pp = table();
if isfile(cfg.participants_file)
    Pp = readtable(cfg.participants_file, 'FileType', 'text', 'Delimiter', '\t', ...
                   'TextType', 'string', 'TreatAsMissing', {'n/a', 'NA'});
else
    fprintf('participants.tsv not found: covariate part skipped.\n');
end

%% ========================================================================
%  2. RUNS
%  ========================================================================
R = cell(numel(cfg.runs), 1);
for k = 1:numel(cfg.runs)
    fprintf('\n===== Run %d/%d: %s =====\n', k, numel(cfg.runs), cfg.runs(k).name);
    R{k} = run_sensitivity(W, BP, RT, PT, Pp, cfg, cfg.runs(k));
end
if numel(R) > 1
    compare_runs(R, cfg);
end
fprintf('\nDone. Outputs in:\n  %s\n', cfg.out_root);


%% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================

% -------------------------------------------------------------------------
function fn = spectra_file(cfg, sub)
fn = fullfile(cfg.spec_root, char(sub), [char(sub) '_task-rest_desc-spectra.mat']);
end

% -------------------------------------------------------------------------
function out = run_sensitivity(W, BP, RT, PT, Pp, cfg, run)
% RUN_SENSITIVITY  The whole sensitivity analysis for one subject set
rng(cfg.seed);
od = fullfile(cfg.out_root, run.name);
if ~isfolder(fullfile(od, 'cases')), mkdir(fullfile(od, 'cases')); end

% subject set of this run
keep = ~ismember(W.subject, string(run.exclude));
W  = W(keep, :);
subs = W.subject;
BP = BP(ismember(BP.subject, subs), :);
RT = RT(ismember(RT.subject, subs), :);
PT = PT(ismember(PT.subject, subs), :);
nS = numel(subs);
fprintf('%d subjects.\n', nS);
SP = cell(nS, 1);
for i = 1:nS
    SP{i} = load(spectra_file(cfg, subs(i)), 'f', 'psd', 'rois', 'subject');
end

% ---- screening ----
[X, names, meas, roi, cond] = measure_matrix(W, cfg);
Z = robust_z(X);
nV = numel(names);
scr = table(repmat(subs, nV, 1), repelem(meas(:), nS), repelem(roi(:), nS), repelem(cond(:), nS), ...
            X(:), Z(:), abs(Z(:)) > cfg.z_thr, 'VariableNames', ...
            {'subject','measure','roi','condition','value','robust_z','extreme'});
writetable(scr, fullfile(od, 'table_screening.csv'));
A = artifact_indices(SP, W, PT, cfg);
writetable(A, fullfile(od, 'table_artifact_indices.csv'));
heatmap_figure(Z, subs, names, cfg, fullfile(od, 'fig_screening_heatmap.png'));

% ---- EC - EO effects and influence ----
[D, enames, dEO, dEC] = effect_matrix(W, cfg);
[E, infl] = effect_stats(D, enames, subs, cfg);
E.r_EO_EC = diag(corr_cols(dEO, dEC));
writetable(E, fullfile(od, 'table_effects.csv'));
writetable(infl, fullfile(od, 'table_influence.csv'));

% ---- associations with subject variables ----
COV = table();
if ~isempty(Pp)
    [C, cn, ctype, cnote] = covariate_matrix(Pp, subs, cfg);
    for j = 1:numel(cn)
        T = covariate_stats(X, meas, roi, cond, C(:, j), cn(j), ctype(j), cnote(j), subs, cfg);
        COV = [COV; T]; %#ok<AGROW>
        covariate_figure(T, cfg, fullfile(od, sprintf('fig_covariate_%s.png', cn(j))));
    end
    writetable(COV, fullfile(od, 'table_covariates.csv'));
end

% ---- subject multiverse ----
sets = cfg.exclusion_sets;
if cfg.add_all_listed_set
    sets(end+1).name = 'all_listed';
    sets(end).subjects = unique([cfg.exclusion_sets.subjects]);
    sets(end).reason = 'all subjects of the sets above removed together';
end
XS = table();
for k = 1:numel(sets)
    kp = ~ismember(subs, string(sets(k).subjects));
    for j = 1:size(D, 2)
        [dz, t, p] = paired_stats(D(kp, j));
        XS = [XS; table(string(sets(k).name), string(sets(k).reason), sum(kp), enames(j), ...
            dz, t, p, p < cfg.alpha, 'VariableNames', {'exclusion_set','reason','n','effect', ...
            'dz','t','p','significant'})]; %#ok<AGROW>
    end
end
writetable(XS, fullfile(od, 'table_exclusion_sets.csv'));
exclusion_figure(XS, sets, fullfile(od, 'fig_exclusion_sets.png'));

% ---- analysis multiverse ----
SC = spec_curve(BP, RT, subs, sets, cfg);
writetable(SC, fullfile(od, 'table_specification_curve.csv'));
for c = cfg.spec_constructs
    speccurve_figure(SC(SC.construct == c{1}, :), c{1}, cfg, ...
                     fullfile(od, sprintf('fig_speccurve_%s.png', c{1})));
end

% ---- subject summary and cases ----
S = subject_summary(subs, scr, A, infl, W, cfg);
writetable(S, fullfile(od, 'table_subject_summary.csv'));
influence_figure(infl, subs, fullfile(od, 'fig_influence.png'));
for c = cfg.case_subjects
    i = find(subs == c{1});
    if isempty(i), continue; end                  % not in this run: skip
    case_study(i, subs, SP, PT, Z, names, A, W, S, cfg, od);
end
write_readme(fullfile(od, 'sensitivity_README.txt'), cfg, run);
out = struct('name', run.name, 'E', E, 'COV', COV, 'n', nS);
end

% -------------------------------------------------------------------------
function [X, names, meas, roi, cond] = measure_matrix(W, cfg)
% MEASURE_MATRIX  Subjects x variables matrix of the main measures (EO, EC,
% total and the reactivity EC - EO), log10 where needed. For alpha_abs the
% reactivity is log10 of the Berger ratio; for alpha_per it is the periodic
% Berger effect.
X = []; names = strings(0); meas = strings(0); roi = strings(0); cond = strings(0);
for r = cfg.rois
    for m = cfg.measures
        v = struct(); ok = true;
        for c = {'EO', 'EC', 'total'}
            col = sprintf('%s_%s_%s', c{1}, r{1}, m{1});
            if ~ismember(col, W.Properties.VariableNames), ok = false; break; end
            x = W.(col);
            if ismember(m{1}, cfg.log_measures), x = log10(x); end
            v.(c{1}) = x;
        end
        if ~ok, continue; end
        v.EC_minus_EO = v.EC - v.EO;
        for c = {'EO', 'EC', 'total', 'EC_minus_EO'}
            X = [X, v.(c{1})]; %#ok<AGROW>
            names(end+1) = sprintf('%s %s %s', m{1}, r{1}, c{1}); %#ok<AGROW>
            meas(end+1) = m{1}; roi(end+1) = r{1}; cond(end+1) = c{1}; %#ok<AGROW>
        end
    end
end
end

% -------------------------------------------------------------------------
function Z = robust_z(X)
% ROBUST_Z  (x - median) / (1.4826 * MAD), column-wise (Leys et al. 2013)
med = median(X, 1, 'omitnan');
mad = 1.4826 * median(abs(X - med), 1, 'omitnan');
mad(mad == 0) = NaN;
Z = (X - med) ./ mad;
end

% -------------------------------------------------------------------------
function A = artifact_indices(SP, W, PT, cfg)
% ARTIFACT_INDICES  Residual-artifact and eye-closure indices (one row per
% subject). Their robust z is one-sided: only an excess is an artifact
% (negative z set to 0).
n = numel(SP);
ocular = NaN(n, 1); slowfront = NaN(n, 1); muscle_EO = NaN(n, 1); muscle_EC = NaN(n, 1);
for i = 1:n
    f = SP{i}.f(:)'; Rr = SP{i}.rois;
    lo = f >= 1 & f <= 3;
    lp = @(P, idx, sel) log10(mean(P(idx, sel), 'all'));
    ocular(i)    = lp(SP{i}.psd.EO, Rr.frontal, lo) - lp(SP{i}.psd.EC, Rr.frontal, lo);
    slowfront(i) = lp(SP{i}.psd.EO, Rr.frontal, lo) - lp(SP{i}.psd.EO, Rr.posterior, lo);
    muscle_EO(i) = hump_db(SP{i}.psd.EO(Rr.global, :), f);
    muscle_EC(i) = hump_db(SP{i}.psd.EC(Rr.global, :), f);
end
% eye closure: suspicious EC periods, VEOG events in EC
ec_blinks = NaN(n, 1); n_noalpha = NaN(n, 1); noalpha_list = strings(n, 1);
for i = 1:n
    t = PT(PT.subject == W.subject(i) & PT.roi == "posterior", :);
    ec = t(t.condition == "EC", :); eo = t(t.condition == "EO", :);
    ec_blinks(i) = mean(ec.blinks_per_min, 'omitnan');
    ref = median(log10(eo.alpha_abs_uV2), 'omitnan') + cfg.ec_alpha_margin;
    bad = log10(ec.alpha_abs_uV2) <= ref & ec.blinks_per_min >= cfg.blink_abs_thr;
    n_noalpha(i) = sum(bad);
    if any(bad), noalpha_list(i) = strjoin(compose('EC%d', ec.period(bad)'), ' '); end
end
A = table(W.subject, ocular, slowfront, muscle_EO, muscle_EC, ec_blinks, n_noalpha, noalpha_list, ...
    'VariableNames', {'subject','ocular_frontal_1_3Hz_EO_minus_EC_log10', ...
    'frontal_minus_posterior_1_3Hz_EO_log10','muscle_29_36Hz_excess_EO_dB', ...
    'muscle_29_36Hz_excess_EC_dB','VEOG_events_per_min_EC','n_suspicious_EC_periods', ...
    'suspicious_EC_periods'});
Zr = robust_z(A{:, 2:6});
Zr(Zr < 0) = 0;                                  % one-sided
Zr(:, 5) = Zr(:, 5) .* (ec_blinks > cfg.blink_abs_thr);   % events: also above an absolute rate
zn = {'z_ocular','z_frontal_slow','z_muscle_EO','z_muscle_EC','z_VEOG_events_EC'};
for k = 1:numel(zn), A.(zn{k}) = Zr(:, k); end
end

% -------------------------------------------------------------------------
function d = hump_db(P, f)
% HUMP_DB  Median over channels of the mean excess (dB) at 29-36 Hz above a
% straight line in log-log coordinates joining 22-26 Hz and 38-42 Hz
L = log10(P);
a = f >= 22 & f <= 26; b = f >= 38 & f <= 42; h = f >= 29 & f <= 36;
xa = mean(log10(f(a))); xb = mean(log10(f(b)));
ya = mean(L(:, a), 2);  yb = mean(L(:, b), 2);
line = ya + (yb - ya) .* (log10(f(h)) - xa) / (xb - xa);
d = median(mean(L(:, h) - line, 2)) * 10;
end

% -------------------------------------------------------------------------
function [D, enames, dEO, dEC] = effect_matrix(W, cfg)
% EFFECT_MATRIX  Reference effects: within-subject difference EC - EO of
% every main measure x region (log10 where needed). One row per subject.
% dEO / dEC are the two conditions (used for the EO-EC correlation).
% REPLACE THIS FUNCTION to test the real hypotheses of a paper.
D = []; dEO = []; dEC = []; enames = strings(0);
for r = cfg.rois
    for m = cfg.measures
        cEO = sprintf('EO_%s_%s', r{1}, m{1}); cEC = sprintf('EC_%s_%s', r{1}, m{1});
        if ~all(ismember({cEO, cEC}, W.Properties.VariableNames)), continue; end
        a = W.(cEO); b = W.(cEC);
        if ismember(m{1}, cfg.log_measures), a = log10(a); b = log10(b); end
        D = [D, b - a]; dEO = [dEO, a]; dEC = [dEC, b]; %#ok<AGROW>
        enames(end+1) = sprintf('%s %s EC-EO', m{1}, r{1}); %#ok<AGROW>
    end
end
end

% -------------------------------------------------------------------------
function R = corr_cols(A, B)
% CORR_COLS  Pearson correlation between every column of A and of B
A = A - mean(A, 1, 'omitnan'); B = B - mean(B, 1, 'omitnan');
R = (A' * B) ./ sqrt(sum(A.^2, 1)' * sum(B.^2, 1));
end

% -------------------------------------------------------------------------
function [dz, t, p] = paired_stats(d)
% PAIRED_STATS  Cohen's dz, t and two-sided p of a one-sample t-test on d
d = d(isfinite(d));
n = numel(d);
if n < 3, dz = NaN; t = NaN; p = NaN; return; end
dz = mean(d) / std(d);
t  = dz * sqrt(n);
p  = t_pvalue(t, n - 1);
end

% -------------------------------------------------------------------------
function p = t_pvalue(t, df)
% T_PVALUE  Two-sided p of Student's t (no Statistics Toolbox needed)
p = betainc(df ./ (df + t.^2), df / 2, 0.5);
end

% -------------------------------------------------------------------------
function [E, infl] = effect_stats(D, enames, subs, cfg)
% EFFECT_STATS  For every effect: full-sample statistics, sign-flip
% permutation p, bootstrap CI of dz, leave-one-out and leave-two-out.
[n, nE] = size(D);
pairs = nchoosek(1:n, 2);
np = size(pairs, 1);
signs = sign(rand(n, cfg.n_perm) - 0.5);          % random sign flips
bidx  = randi(n, n, cfg.n_boot);                  % bootstrap resamples
E = table(); infl = table();
for j = 1:nE
    d = D(:, j);
    [dz, t, p] = paired_stats(d);
    % sign-flip permutation (Monte Carlo, observed value included)
    ds = d .* signs;
    tperm = mean(ds, 1) ./ (std(ds, 0, 1) / sqrt(n));
    p_perm = (sum(abs(tperm) >= abs(t)) + 1) / (cfg.n_perm + 1);
    % bootstrap CI of dz (percentile)
    db = d(bidx);
    ci = prctile_local(mean(db, 1) ./ std(db, 0, 1), [2.5 97.5]);
    % leave-one-out
    dz1 = NaN(n, 1); p1 = NaN(n, 1);
    for i = 1:n
        [dz1(i), ~, p1(i)] = paired_stats(d([1:i-1, i+1:n]));
    end
    % leave-two-out: one column per removed pair, removed values = NaN
    dd = repmat(d, 1, np);
    dd(sub2ind([n np], pairs(:, 1)', 1:np)) = NaN;
    dd(sub2ind([n np], pairs(:, 2)', 1:np)) = NaN;
    dz2 = mean(dd, 1, 'omitnan') ./ std(dd, 0, 1, 'omitnan');
    p2  = t_pvalue(dz2 * sqrt(n - 2), n - 3);
    [~, imax] = max(abs(dz1 - dz));
    E = [E; table(enames(j), n, mean(d), dz, ci(1), ci(2), t, p, p_perm, ...
        min(dz1), max(dz1), mean(p1 < cfg.alpha), subs(imax), dz1(imax) - dz, ...
        min(dz2), max(dz2), mean(p2 < cfg.alpha), ...
        'VariableNames', {'effect','n','mean_diff','dz','dz_boot_ci_low','dz_boot_ci_high', ...
        't','p_ttest','p_signflip','dz_loo_min','dz_loo_max','prop_sig_loo', ...
        'most_influential_subject','dz_change_without_it','dz_l2o_min','dz_l2o_max', ...
        'prop_sig_l2o'})]; %#ok<AGROW>
    ch = dz1 - dz;
    zi = (ch - median(ch)) / (1.4826 * median(abs(ch - median(ch))) + eps);
    infl = [infl; table(subs, repmat(enames(j), n, 1), d, ch, zi, ...
        (p1 < cfg.alpha) ~= (p < cfg.alpha), 'VariableNames', {'subject','effect', ...
        'subject_value','dz_change_if_removed','influence_z','significance_changes'})]; %#ok<AGROW>
end
end

% -------------------------------------------------------------------------
function q = prctile_local(x, p)
% PRCTILE_LOCAL  Percentiles with linear interpolation (no toolbox needed)
x = sort(x(isfinite(x))); n = numel(x);
if n == 0, q = NaN(size(p)); return; end
if n == 1, q = repmat(x, size(p)); return; end
q = interp1(1:n, x, 1 + (n - 1) * p / 100);
end

% -------------------------------------------------------------------------
function [C, cn, ctype, cnote] = covariate_matrix(Pp, subs, cfg)
% COVARIATE_MATRIX  Subject variables aligned with SUBS. Numeric columns
% are used as they are; text columns with 2 levels are coded 0/1 (the
% coding is written in the output); other text columns are skipped.
% Participants that are not in SUBS are ignored.
C = []; cn = strings(0); ctype = strings(0); cnote = strings(0);
[tf, loc] = ismember(subs, Pp.participant_id);
if any(~tf)
    fprintf('Not in participants.tsv (covariates = NaN): %s\n', strjoin(subs(~tf)', ', '));
end
lists = {cfg.covariates, cfg.nuisance}; types = ["interest", "technical"];
for L = 1:2
    for v = lists{L}
        if ~ismember(v{1}, Pp.Properties.VariableNames)
            fprintf('Column "%s" not in participants.tsv: skipped.\n', v{1}); continue
        end
        col = Pp.(v{1});
        x = NaN(numel(subs), 1);
        note = "numeric";
        if isnumeric(col)
            x(tf) = col(loc(tf));
        else
            col = string(col);
            lev = unique(col(~ismissing(col)));
            if numel(lev) ~= 2
                fprintf('Column "%s" has %d levels: skipped.\n', v{1}, numel(lev)); continue
            end
            xx = double(col == lev(2)); xx(ismissing(col)) = NaN;
            x(tf) = xx(loc(tf));
            note = sprintf('coded 0 = %s, 1 = %s', lev(1), lev(2));
        end
        if sum(isfinite(x)) < 6 || numel(unique(x(isfinite(x)))) < 2
            fprintf('Column "%s" has too few values: skipped.\n', v{1}); continue
        end
        C = [C, x]; cn(end+1) = v{1}; ctype(end+1) = types(L); cnote(end+1) = note; %#ok<AGROW>
    end
end
end

% -------------------------------------------------------------------------
function r = r_masked(x, y, K)
% R_MASKED  Pearson r of x and y using, for every column of the logical
% matrix K, only the rows marked true (vectorized leave-k-out)
K = double(K);
nk = sum(K, 1);
mx = sum(x .* K, 1) ./ nk;  my = sum(y .* K, 1) ./ nk;
dx = (x - mx) .* K;         dy = (y - my) .* K;
r = sum(dx .* dy, 1) ./ sqrt(sum(dx.^2, 1) .* sum(dy.^2, 1));
end

% -------------------------------------------------------------------------
function T = covariate_stats(X, meas, roi, cond, c, cname, ctype, cnote, subs, cfg)
% COVARIATE_STATS  Correlation of every measure with one covariate:
% Pearson r and p, FDR q (Benjamini-Hochberg over all the measures of this
% covariate), Spearman rho, leave-one-out and leave-two-out ranges.
nV = size(X, 2);
r = NaN(nV, 1); p = r; rho = r; n = r; r1min = r; r1max = r; ps1 = r; r2min = r; r2max = r; ps2 = r;
infl_sub = strings(nV, 1); dr = r;
for j = 1:nV
    good = isfinite(X(:, j)) & isfinite(c);
    x = X(good, j); y = c(good); s = subs(good); m = numel(x);
    if m < 6, continue; end
    n(j) = m;
    r(j) = r_masked(x, y, true(m, 1));
    p(j) = t_pvalue(r(j) * sqrt((m - 2) / max(1 - r(j)^2, eps)), m - 2);
    rho(j) = r_masked(ranks(x), ranks(y), true(m, 1));
    K1 = ~eye(m);                                         % leave-one-out
    rr1 = r_masked(x, y, K1);
    pp1 = t_pvalue(rr1 .* sqrt((m - 3) ./ max(1 - rr1.^2, eps)), m - 3);
    pr = nchoosek(1:m, 2); K2 = true(m, size(pr, 1));     % leave-two-out
    K2(sub2ind(size(K2), pr(:, 1)', 1:size(pr, 1))) = false;
    K2(sub2ind(size(K2), pr(:, 2)', 1:size(pr, 1))) = false;
    rr2 = r_masked(x, y, K2);
    pp2 = t_pvalue(rr2 .* sqrt((m - 4) ./ max(1 - rr2.^2, eps)), m - 4);
    r1min(j) = min(rr1); r1max(j) = max(rr1); ps1(j) = mean(pp1 < cfg.alpha);
    r2min(j) = min(rr2); r2max(j) = max(rr2); ps2(j) = mean(pp2 < cfg.alpha);
    [~, k] = max(abs(rr1 - r(j))); infl_sub(j) = s(k); dr(j) = rr1(k) - r(j);
end
q = fdr_bh(p);
T = table(meas(:), roi(:), cond(:), repmat(cname, nV, 1), repmat(ctype, nV, 1), ...
    repmat(cnote, nV, 1), n, r, p, q, rho, r1min, r1max, ps1, infl_sub, dr, r2min, r2max, ps2, ...
    'VariableNames', {'measure','roi','condition','covariate','covariate_type','coding','n', ...
    'pearson_r','p','q_fdr','spearman_rho','r_loo_min','r_loo_max','prop_sig_loo', ...
    'most_influential_subject','r_change_without_it','r_l2o_min','r_l2o_max','prop_sig_l2o'});
end

% -------------------------------------------------------------------------
function q = fdr_bh(p)
% FDR_BH  Benjamini-Hochberg adjusted p-values (q), NaN ignored
q = NaN(size(p));
ok = find(isfinite(p)); m = numel(ok);
if m == 0, return; end
[ps, o] = sort(p(ok));
qs = ps .* m ./ (1:m)';
qs = min(1, flipud(cummin(flipud(qs))));
q(ok(o)) = qs;
end

% -------------------------------------------------------------------------
function rk = ranks(x)
% RANKS  Ranks with ties averaged (for Spearman's rho)
[~, o] = sort(x(:)); rk = zeros(numel(x), 1); rk(o) = 1:numel(x);
[u, ~, g] = unique(x(:));
for k = 1:numel(u)
    if sum(g == k) > 1, rk(g == k) = mean(rk(g == k)); end
end
end

% -------------------------------------------------------------------------
function SC = spec_curve(BP, RT, subs, sets, cfg)
% SPEC_CURVE  EC - EO effect of every construct for every combination of
% band set x measure type x region x exclusion set.
SC = table();
for c = cfg.spec_constructs
    if ismember(c{1}, {'TBR', 'ATR', 'ABR'})
        T = RT(RT.ratio == c{1}, :); tl = {'log10_ratio_abs', 'log10_ratio_periodic'};
    else
        T = BP(BP.band == c{1}, :);  tl = {'log10_abs_power', 'rel_power', 'periodic_log10'};
    end
    for bs = ["iaf", "fixed"]
        for mt = tl
            for r = cfg.spec_rois
                sel = T.band_set == bs & T.roi == r{1};
                a = value_by_subject(T(sel & T.condition == "EO", :), mt{1}, subs);
                b = value_by_subject(T(sel & T.condition == "EC", :), mt{1}, subs);
                if strcmp(mt{1}, 'rel_power'), a = log10(a); b = log10(b); end
                d = b - a;
                for k = 1:numel(sets)
                    kp = ~ismember(subs, string(sets(k).subjects));
                    [dz, ~, p] = paired_stats(d(kp));
                    SC = [SC; table(string(c{1}), bs, string(mt{1}), string(r{1}), ...
                        string(sets(k).name), sum(kp & isfinite(d)), dz, p, p < cfg.alpha, ...
                        'VariableNames', {'construct','band_set','measure','roi', ...
                        'exclusion_set','n','dz','p','significant'})]; %#ok<AGROW>
                end
            end
        end
    end
end
end

% -------------------------------------------------------------------------
function v = value_by_subject(T, col, subs)
v = NaN(numel(subs), 1);
[tf, loc] = ismember(subs, T.subject);
v(tf) = T.(col)(loc(tf));
end

% -------------------------------------------------------------------------
function S = subject_summary(subs, scr, A, infl, W, cfg)
% SUBJECT_SUMMARY  One row per subject with everything that may make it "special"
n = numel(subs);
n_ext = zeros(n, 1); pct_ext = zeros(n, 1); top = strings(n, 1);
inf_med = zeros(n, 1); n_sigchange = zeros(n, 1); n_top_infl = zeros(n, 1);
eff = unique(infl.effect);
most = strings(numel(eff), 1);
for e = 1:numel(eff)
    t = infl(infl.effect == eff(e), :);
    [~, k] = max(abs(t.dz_change_if_removed)); most(e) = t.subject(k);
end
for i = 1:n
    s = scr(scr.subject == subs(i), :);
    n_ext(i) = sum(s.extreme); pct_ext(i) = 100 * mean(s.extreme);
    [~, o] = sort(abs(s.robust_z), 'descend', 'MissingPlacement', 'last');
    o = o(1:min(5, numel(o)));
    top(i) = strjoin(compose('%s %s %s (z=%.1f)', s.measure(o), s.roi(o), s.condition(o), s.robust_z(o))', '; ');
    t = infl(infl.subject == subs(i), :);
    inf_med(i) = median(abs(t.influence_z), 'omitnan');
    n_sigchange(i) = sum(t.significance_changes);
    n_top_infl(i) = sum(most == subs(i));
end
zn = A.Properties.VariableNames(startsWith(A.Properties.VariableNames, 'z_'));
za = A{:, zn};
art = strings(n, 1);
for i = 1:n
    hit = za(i, :) > cfg.z_thr;
    if any(hit), art(i) = strjoin(string(zn(hit)), ' '); end
end
S = table(subs, n_ext, pct_ext, top, art, A.n_suspicious_EC_periods, ...
    A.suspicious_EC_periods, inf_med, n_top_infl, n_sigchange, W.flags, ...
    'VariableNames', {'subject', ...
    'n_extreme_values','pct_extreme_values','most_extreme_values','artifact_index_extreme', ...
    'n_suspicious_EC_periods','suspicious_EC_periods', ...
    'median_abs_influence_z','n_effects_most_influential','n_effects_significance_changes', ...
    'spectral_flags'});
end

% -------------------------------------------------------------------------
function heatmap_figure(Z, subs, names, cfg, fn)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [20 20 1800 900]);
try, theme(fig, 'light'); catch, end
imagesc(max(min(Z, 6), -6)'); colormap(blue_white_red(256)); clim([-6 6]);
cb = colorbar; cb.Label.String = 'robust z';
set(gca, 'XTick', 1:numel(subs), 'XTickLabel', subs, 'XTickLabelRotation', 90, ...
    'YTick', 1:numel(names), 'YTickLabel', names, 'FontSize', 5, 'TickLabelInterpreter', 'none');
hold on
[r, c] = find(abs(Z') > cfg.z_thr);
plot(c, r, 'k.', 'MarkerSize', 6);
title(sprintf('Robust z of every measure (dot = |z| > %.1f)', cfg.z_thr));
exportgraphics(fig, fn, 'Resolution', 150, 'BackgroundColor', 'white');
close(fig);
end

% -------------------------------------------------------------------------
function covariate_figure(T, cfg, fn)
% COVARIATE_FIGURE  Pearson r of every measure x region (rows) with the
% covariate in every condition (columns). o = p < .05; * = q(FDR) < .05.
conds = ["EO", "EC", "total", "EC_minus_EO"];
key = unique(T.measure + " " + T.roi, 'stable');
M = NaN(numel(key), 4); P = M; Q = M;
for k = 1:4
    t = T(T.condition == conds(k), :);
    [~, loc] = ismember(key, t.measure + " " + t.roi);
    ok = loc > 0;
    M(ok, k) = t.pearson_r(loc(ok)); P(ok, k) = t.p(loc(ok)); Q(ok, k) = t.q_fdr(loc(ok));
end
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [20 20 700 1500]);
try, theme(fig, 'light'); catch, end
imagesc(M); colormap(blue_white_red(256)); clim([-0.7 0.7]); cb = colorbar; cb.Label.String = 'Pearson r';
hold on
[r, c] = find(P < cfg.alpha); plot(c, r, 'ko', 'MarkerSize', 3);
[r, c] = find(Q < cfg.alpha); plot(c, r, 'k*', 'MarkerSize', 6);
set(gca, 'XTick', 1:4, 'XTickLabel', {'EO', 'EC', 'total', 'EC - EO'}, 'YTick', 1:numel(key), ...
    'YTickLabel', key, 'FontSize', 6, 'TickLabelInterpreter', 'none');
title({sprintf('Correlation with %s (%s)', T.covariate(1), T.coding(1)), ...
       'o = p < .05 (uncorrected), * = q < .05 (FDR)'}, 'Interpreter', 'none');
exportgraphics(fig, fn, 'Resolution', 130, 'BackgroundColor', 'white');
close(fig);
end

% -------------------------------------------------------------------------
function influence_figure(infl, subs, fn)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [20 20 1400 500]);
try, theme(fig, 'light'); catch, end
hold on
for i = 1:numel(subs)
    v = infl.influence_z(infl.subject == subs(i));
    scatter(i + 0.15 * randn(size(v)), v, 8, [0.4 0.4 0.4], 'filled', 'MarkerFaceAlpha', 0.4);
    plot(i, median(abs(v), 'omitnan'), 'rd', 'MarkerFaceColor', 'r');
end
yline([-3.5 3.5], 'k--');
set(gca, 'XTick', 1:numel(subs), 'XTickLabel', subs, 'XTickLabelRotation', 90, ...
    'TickLabelInterpreter', 'none');
ylabel('influence z (change of dz when removed)'); grid on
title('Influence of every subject on every EC - EO effect (red = median |z|)');
exportgraphics(fig, fn, 'Resolution', 120, 'BackgroundColor', 'white');
close(fig);
end

% -------------------------------------------------------------------------
function exclusion_figure(XS, sets, fn)
eff = unique(XS.effect, 'stable');
M = NaN(numel(eff), numel(sets)); S = false(size(M));
for k = 1:numel(sets)
    t = XS(XS.exclusion_set == sets(k).name, :);
    [~, loc] = ismember(eff, t.effect);
    M(:, k) = t.dz(loc); S(:, k) = t.significant(loc);
end
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [20 20 900 1600]);
try, theme(fig, 'light'); catch, end
lim = max(abs(M(:)));
imagesc(M); colormap(blue_white_red(256)); clim([-lim lim]); cb = colorbar; cb.Label.String = 'dz (EC - EO)';
hold on
[r, c] = find(~S); plot(c, r, 'kx', 'MarkerSize', 4);
set(gca, 'XTick', 1:numel(sets), 'XTickLabel', {sets.name}, 'XTickLabelRotation', 45, ...
    'YTick', 1:numel(eff), 'YTickLabel', eff, 'FontSize', 5, 'TickLabelInterpreter', 'none');
title('EC - EO effect for every exclusion set (x = not significant)');
exportgraphics(fig, fn, 'Resolution', 150, 'BackgroundColor', 'white');
close(fig);
end

% -------------------------------------------------------------------------
function speccurve_figure(T, construct, cfg, fn)
% SPECCURVE_FIGURE  Specification curve (Simonsohn et al. 2020): effects
% sorted by size (top) and the choices that produced each one (bottom).
T = sortrows(T, 'dz');
n = height(T);
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [20 20 1500 900]);
try, theme(fig, 'light'); catch, end
subplot(4, 1, 1); hold on
col = repmat([0.6 0.6 0.6], n, 1);
col(T.significant & T.dz > 0, :) = repmat([0.85 0.15 0.15], sum(T.significant & T.dz > 0), 1);
col(T.significant & T.dz < 0, :) = repmat([0.1 0.3 0.9], sum(T.significant & T.dz < 0), 1);
scatter(1:n, T.dz, 6, col, 'filled'); yline(0, 'k-');
xlim([0 n + 1]); ylabel('dz (EC - EO)'); grid on
title(sprintf('%s: %d specifications, %.0f %% significant (red/blue = p < %.2f)', construct, n, ...
      100 * mean(T.significant), cfg.alpha), 'Interpreter', 'none');
dims = {'band_set', 'measure', 'exclusion_set'};
for d = 1:3
    subplot(4, 1, d + 1);
    [lev, ~, g] = unique(T.(dims{d}));
    scatter(1:n, g, 4, 'k', 'filled');
    set(gca, 'YTick', 1:numel(lev), 'YTickLabel', lev, 'TickLabelInterpreter', 'none', 'FontSize', 7);
    ylim([0.5 numel(lev) + 0.5]); xlim([0 n + 1]); ylabel(dims{d}, 'Interpreter', 'none'); grid on
end
exportgraphics(fig, fn, 'Resolution', 120, 'BackgroundColor', 'white');
close(fig);
end

% -------------------------------------------------------------------------
function case_study(i, subs, SP, PT, Z, names, A, W, S, cfg, od)
% CASE_STUDY  Figure + text summary of one subject compared with the group
sub = subs(i);
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [20 20 1600 950]);
try, theme(fig, 'light'); catch, end
f = SP{1}.f(:)';
cols = struct('EO', [0.1 0.3 0.9], 'EC', [0.85 0.15 0.15]);
subplot(2, 3, 1); hold on
for c = {'EO', 'EC'}
    G = cell2mat(cellfun(@(x) log10(mean(x.psd.(c{1})(x.rois.(cfg.case_roi), :), 1)), SP, ...
                 'UniformOutput', false));
    m = mean(G, 1, 'omitnan'); s = std(G, 0, 1, 'omitnan');
    fill([f fliplr(f)], [m - 2 * s fliplr(m + 2 * s)], cols.(c{1}), 'FaceAlpha', 0.12, ...
         'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(f, G(i, :), 'Color', cols.(c{1}), 'LineWidth', 1.8, 'DisplayName', c{1});
end
xlim([1 40]); grid on; legend('Location', 'northeast');
xlabel('Hz'); ylabel('log_{10} power');
title({sprintf('Spectra, %s region', cfg.case_roi), 'band = group mean \pm 2 SD'});
vars = {'alpha_periodic_log10', 'theta_abs_uV2', 'beta_abs_uV2', 'ap_exponent', 'blinks_per_min'};
ttl  = {'periodic alpha (log_{10})', 'log_{10} theta', 'log_{10} beta', 'aperiodic exponent', ...
        'VEOG events / min'};
P = PT(PT.roi == cfg.case_roi, :);
for k = 1:numel(vars)
    subplot(2, 3, k + 1); hold on
    g = NaN(numel(subs), 8);
    for j = 1:numel(subs)
        t = P(P.subject == subs(j), :);
        g(j, t.period) = t.(vars{k})';
    end
    if endsWith(vars{k}, '_uV2'), g = log10(g); end
    tp = P(P.subject == sub, :);
    pc = strings(1, 8); pc(tp.period) = tp.condition';
    for q = 1:8                       % group reference: same condition in that period
        same = false(numel(subs), 1);
        for j = 1:numel(subs)
            t = P(P.subject == subs(j) & P.period == q, :);
            same(j) = ~isempty(t) && t.condition(1) == pc(q);
        end
        qq = prctile_local(g(same, q), [25 50 75]);
        plot([q q], qq([1 3]), '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 6);
        plot(q, qq(2), 'o', 'Color', [0.3 0.3 0.3]);
    end
    plot(1:8, g(i, :), 'k-');
    for q = 1:8
        if pc(q) == "", continue; end
        plot(q, g(i, q), 's', 'MarkerSize', 9, 'MarkerFaceColor', cols.(char(pc(q))), ...
             'MarkerEdgeColor', 'k');
    end
    set(gca, 'XTick', 1:8, 'XTickLabel', compose('%d %s', (1:8)', pc'));
    xlim([0.5 8.5]); grid on
    title({ttl{k}, 'grey = group IQR (same condition)'});
end
sgtitle(sprintf('Case study %s', sub), 'Interpreter', 'none');
exportgraphics(fig, fullfile(od, 'cases', [char(sub) '_case.png']), ...
               'Resolution', 110, 'BackgroundColor', 'white');
close(fig);
% text summary
fid = fopen(fullfile(od, 'cases', [char(sub) '_case.txt']), 'w');
fprintf(fid, 'CASE STUDY %s\n%s\n', sub, repmat('=', 1, 40));
fprintf(fid, 'Spectral flags: %s\n', W.flags(i));
fprintf(fid, 'Suspicious EC periods (eyes maybe open): %d (%s)\n', A.n_suspicious_EC_periods(i), ...
        A.suspicious_EC_periods(i));
fprintf(fid, 'IAF reference: %.2f Hz (%s)\n', W.IAF_ref_hz(i), W.IAF_ref_source(i));
fprintf(fid, '\nArtifact indices (value, one-sided robust z):\n');
vn = A.Properties.VariableNames;
zc = find(startsWith(vn, 'z_'));
for k = 1:numel(zc)
    fprintf(fid, '  %-45s %8.3f   z = %5.1f\n', vn{k + 1}, A{i, k + 1}, A{i, zc(k)});
end
fprintf(fid, '\nExtreme values (|robust z| > %.1f): %d of %d measures\n', cfg.z_thr, ...
        S.n_extreme_values(i), size(Z, 2));
[~, o] = sort(abs(Z(i, :)), 'descend', 'MissingPlacement', 'last');
for k = o(1:min(15, numel(o)))
    fprintf(fid, '  %-40s z = %6.2f\n', names(k), Z(i, k));
end
fprintf(fid, '\nInfluence: median |influence z| %.2f; most influential subject in %d effects; ', ...
        S.median_abs_influence_z(i), S.n_effects_most_influential(i));
fprintf(fid, 'removing it changes significance in %d effects.\n', S.n_effects_significance_changes(i));
fclose(fid);
end

% -------------------------------------------------------------------------
function compare_runs(R, cfg)
% COMPARE_RUNS  Effects and correlations of the first run (reference)
% side by side with every other run
ref = R{1};
CE = table(); CC = table();
for k = 2:numel(R)
    a = ref.E(:, {'effect','dz','p_ttest'});  b = R{k}.E(:, {'effect','dz','p_ttest'});
    a.Properties.VariableNames = {'effect', 'dz_ref', 'p_ref'};
    b.Properties.VariableNames = {'effect', 'dz_run', 'p_run'};
    J = join(a, b, 'Keys', 'effect');
    J.run = repmat(string(R{k}.name), height(J), 1);
    J.dz_change = J.dz_run - J.dz_ref;
    J.significance_changes = (J.p_ref < cfg.alpha) ~= (J.p_run < cfg.alpha);
    CE = [CE; J]; %#ok<AGROW>
    if ~isempty(ref.COV) && ~isempty(R{k}.COV)
        key = @(T) T.measure + "|" + T.roi + "|" + T.condition + "|" + T.covariate;
        A = ref.COV; B = R{k}.COV;
        [tf, loc] = ismember(key(A), key(B));
        A = A(tf, :); B = B(loc(tf), :);
        CC = [CC; table(A.measure, A.roi, A.condition, A.covariate, repmat(string(R{k}.name), height(A), 1), ...
            A.pearson_r, A.p, A.q_fdr, B.pearson_r, B.p, B.q_fdr, B.pearson_r - A.pearson_r, ...
            (A.p < cfg.alpha) ~= (B.p < cfg.alpha), (A.q_fdr < cfg.alpha) ~= (B.q_fdr < cfg.alpha), ...
            'VariableNames', {'measure','roi','condition','covariate','run','r_ref','p_ref','q_ref', ...
            'r_run','p_run','q_run','r_change','p_significance_changes','q_significance_changes'})]; %#ok<AGROW>
    end
end
writetable(CE, fullfile(cfg.out_root, 'table_compare_runs_effects.csv'));
if ~isempty(CC), writetable(CC, fullfile(cfg.out_root, 'table_compare_runs_covariates.csv')); end
% figure
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [20 20 1200 550]);
try, theme(fig, 'light'); catch, end
subplot(1, 2, 1); hold on
scatter(CE.dz_ref, CE.dz_run, 12, double(CE.significance_changes), 'filled');
lim = [min([CE.dz_ref; CE.dz_run]) max([CE.dz_ref; CE.dz_run])];
plot(lim, lim, 'k--'); grid on; colormap([0.5 0.5 0.5; 0.85 0.15 0.15]);
xlabel(sprintf('dz, %s', ref.name), 'Interpreter', 'none'); ylabel('dz, other run');
title(sprintf('EC - EO effects (red = significance changes: %d)', sum(CE.significance_changes)));
if ~isempty(CC)
    subplot(1, 2, 2); hold on
    scatter(CC.r_ref, CC.r_run, 8, double(CC.p_significance_changes), 'filled');
    plot([-1 1], [-1 1], 'k--'); grid on; axis([-0.8 0.8 -0.8 0.8]);
    xlabel(sprintf('r, %s', ref.name), 'Interpreter', 'none'); ylabel('r, other run');
    title(sprintf('Correlations (red = p < .05 changes: %d)', sum(CC.p_significance_changes)));
end
exportgraphics(fig, fullfile(cfg.out_root, 'fig_compare_runs.png'), 'Resolution', 120, ...
               'BackgroundColor', 'white');
close(fig);
end

% -------------------------------------------------------------------------
function cm = blue_white_red(n)
h = floor(n / 2);
up = linspace(0, 1, h)';
cm = [[up up ones(h, 1)]; [ones(n - h, 1) flipud(linspace(0, 1, n - h)') flipud(linspace(0, 1, n - h)')]];
end

% -------------------------------------------------------------------------
function write_readme(fn, cfg, run)
txt = {
sprintf('SENSITIVITY TABLES - sensitivity_rs.m - run "%s"', run.name)
'=========================================================='
sprintf('Subjects excluded in this run: %s', char(strjoin(["", string(run.exclude)], ' ')))
sprintf('Extreme value: |robust z| > %.1f; robust z = (x - median) / (1.4826 * MAD).', cfg.z_thr)
'Measures: from table_wide.csv (IAF bands); *_abs and TBR/ATR/ABR are log10-transformed.'
'Condition EC_minus_EO = reactivity: for alpha_abs it is log10 of the Berger ratio (EC/EO),'
'  for alpha_per the periodic Berger effect; for the other measures the same idea.'
''
'table_screening.csv: subject, measure, roi, condition (EO, EC, total, EC_minus_EO), value, robust_z, extreme.'
'table_artifact_indices.csv: ocular (frontal 1-3 Hz log10 power EO - EC), frontal_slow (frontal -'
'  posterior 1-3 Hz, EO), muscle (dB excess at 29-36 Hz above the 22-26 / 38-42 Hz line, global),'
'  VEOG events per minute in EC, number and list of suspicious EC periods (posterior log10 alpha'
sprintf('  <= median of EO periods + %.2f AND >= %d VEOG events/min). z_* = ONE-SIDED robust z.', cfg.ec_alpha_margin, cfg.blink_abs_thr)
'table_effects.csv: one row per EC - EO effect: mean_diff, dz (Cohen), bootstrap 95% CI of dz,'
'  t, p_ttest, p_signflip (random sign flips), dz_loo_min/max and prop_sig_loo (leave-one-out),'
'  most_influential_subject, dz_change_without_it, dz_l2o_min/max and prop_sig_l2o (leave-two-out,'
'  all pairs), r_EO_EC (correlation of EO and EC across subjects: high = trait-like ranking).'
'table_influence.csv: subject x effect: change of dz when that subject is removed, influence_z,'
'  significance_changes (removing it flips p < .05).'
'table_covariates.csv: measure x roi x condition x covariate: n, Pearson r, p, q_fdr (Benjamini-'
'  Hochberg within each covariate), Spearman rho, r_loo_min/max, prop_sig_loo, most influential'
'  subject and r change without it, r_l2o_min/max, prop_sig_l2o. covariate_type: interest /'
'  technical (technical variables should not correlate). coding: 0/1 coding of 2-level variables.'
'table_exclusion_sets.csv: every effect recomputed without each pre-specified set of subjects.'
'table_specification_curve.csv: construct x band_set x measure x roi x exclusion_set: n, dz, p.'
'table_subject_summary.csv: one row per subject: extreme values, artifact indices above threshold,'
'  suspicious EC periods, influence summary, spectral flags.'
'fig_covariate_<name>.png: r of every measure x region (rows) in EO, EC, total and EC - EO.'
'../table_compare_runs_*.csv: the first run side by side with the others.'
};
fid = fopen(fn, 'w'); fprintf(fid, '%s\n', txt{:}); fclose(fid);
end
