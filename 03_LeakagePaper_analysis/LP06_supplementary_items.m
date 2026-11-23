function LP06_supplementary_items()
%LP06_SUPPLEMENTARY_ITEMS  Leakage paper - supplementary tables and figures (frozen outputs only).
%
% Tables (one XLSX workbook, one sheet each, plus CSVs):
%   S1  predictors and preprocessing per pipeline (summary, full exclusion list, overlap)
%   S2  thresholds and validation metrics
%   S3  full candidate-level estimates (LP05a) with 95% CIs
%   S4  supporting burden agreement (volume, SUVmean, TL-PSMA)
%   S5  supplementary ablation (candidate level and segmentation)
%   S6  reference-informed reconstruction (descriptive; circular internally)
%   S7  calibration (summary and deciles)
%   S8  deviations register and prespecified analyses not performed
%   S9  frozen runs, hashes and safeguards
%   S10 source data for Figure 4 (acceptance by occupancy; lesion outcomes)
% Figures (greyscale only; PNG/TIFF 600 dpi, PDF, FIG):
%   S1 ROC curves with locked operating points      S4 Dice per study
%   S2 calibration plots                            S5 reconstruction internal vs external
%   S3 volume agreement (scatter, Bland-Altman)     S6 validation threshold curves
% No new statistical analysis: tables re-format frozen estimates; figures plot
% frozen per-candidate / per-study values and frozen estimates.

clc;
V2 = string(fullfile(psma_root("PSMA_PROJECT_ROOT"),"PSMA_AI_Project_V2_PET_Derived"));
R = fullfile(V2,"06_results"); EV = fullfile(R,"LeakagePaper_Evaluation"); MD = fullfile(R,"LeakagePaper_ModelData");
L4 = fullfile(MD,"LP04_20261109_204900"); A5 = fullfile(EV,"LP05a_20261109_210731"); B5 = fullfile(EV,"LP05b_20261109_220537");
f = struct( ...
 "l4sum",  fullfile(L4,"LP04_summary.csv"), ...
 "l4ovl",  fullfile(L4,"LP04_predictor_overlap.csv"), ...
 "l4miss", fullfile(L4,"LP04_ALL_pipelines_missingness_and_exclusions.csv"), ...
 "aEst",   fullfile(A5,"LP05a_estimates_with_bootstrap_CI.csv"), ...
 "aAcc",   fullfile(A5,"LP05a_acceptance_by_occupancy.csv"), ...
 "aCal",   fullfile(A5,"LP05a_calibration_summary.csv"), ...
 "aBins",  fullfile(A5,"LP05a_calibration_bins.csv"), ...
 "aSafe",  fullfile(A5,"LP05a_safeguards.csv"), ...
 "aScI",   fullfile(A5,"LP05a_scores_internal_test_B.csv"), ...
 "aScE",   fullfile(A5,"LP05a_scores_external_B.csv"), ...
 "bEst",   fullfile(B5,"LP05b_estimates_with_bootstrap_CI.csv"), ...
 "bStudy", fullfile(B5,"LP05b_study_level_results.csv"), ...
 "bFail",  fullfile(B5,"LP05b_failure_decomposition.csv"), ...
 "bC1",    fullfile(B5,"LP05b_C1_only_detection.csv"), ...
 "bSafe",  fullfile(B5,"LP05b_safeguards.csv"));
h = struct( ...
 "l4sum","2b8a4ee6d9d2dc8310e2e0adf92e1bc446c5992f51ea2df92e6fe4d9f7e32127", ...
 "l4ovl","523a9943af8b125c140a1f48a6e2b842231a153da7da7ba57c71cad990720044", ...
 "l4miss","8fa9e321620cab99fc590fc0d99668c604096b48fcd0c52d923f5d615df25f8e", ...
 "aEst","2a7c792ca3a96443e523c460eb2ce9569f09e0aa91585755942fc0cd67ab1e56", ...
 "aAcc","2f36134662846e11931ecfb5de550589dda7c9e4c4995d1b324d21acf901a4b8", ...
 "aCal","71b5731c1b12c79a0220b6d5999f2aaf26bed554b5f99b7e19aadeb50eac4337", ...
 "aBins","e5bfea7b015614a8b763c7f72472dc4f2411f81d5e9f876f2152ceccdebb8326", ...
 "aSafe","abf403ebeca5968f6ed665ef00c0fe244c7c4454c26c0d9809daca426999d6e7", ...
 "aScI","7c8932007fa45826fbbb3847a339f43f3a22ff273d753310d4231448f9440470", ...
 "aScE","b33bb7bc80e73c6621532ad7d4f4f57c34bdfbba39f70bda208c78dd44d5c121", ...
 "bEst","88131b7d4f102c70a4d1acc9bf5085b6d6c6d471b7219ebe3df10b58e596fb89", ...
 "bStudy","75260b1e26cbf9c4589391a93130777452a642572c560d34a8ef00e5c6bcb5da", ...
 "bFail","52980c24443a48bd0d9cd299de069ed2c80356b4b10e141bef5275e9a219c6e4", ...
 "bC1","4b02e63e9878d9ea2374e1365ff80e7317a989a53256be56ef396ad0e9259787", ...
 "bSafe","76d78047aab6a411b46e45f5a4ca0a43b648dd6f3e076b87f8e47f717cfd7e93");
for k = string(fieldnames(f))'
    if fileSha256(f.(k)) ~= h.(k), error("Frozen input hash mismatch: %s", f.(k)); end
end
outDir = fullfile(V2,"LeakagePaper","Manuscript","LP06_supp_" + string(datetime("now","Format","yyyyMMdd_HHmmss")));
mkdir(outDir);
diary(char(fullfile(outDir,"LP06_supp_log.txt"))); cleanupDiary = onCleanup(@() diary("off"));
fprintf("LP06 supplementary items %s\nMATLAB %s\nAll %d frozen inputs verified (SHA256).\nOutput: %s\n", ...
    string(datetime("now")), version, numel(fieldnames(f)), outDir);
gridF = struct();
for m = ["B","Astd","Ahist","Ahist_ablation"]
    gridF.(m) = fullfile(L4, "LP04_" + m + "_validation_threshold_grid.csv");
    fprintf("  LP04 grid %-15s SHA256 %s (frozen LP04 folder; logged)\n", m, fileSha256(gridF.(m)));
end

rd = @(p) readtable(p, "TextType","string", "VariableNamingRule","preserve");
xl = fullfile(outDir, "Supplementary_Tables_S1-S10.xlsx");
mods = ["B","Astd","Ahist"]; modLab = ["B (PET-only)","A-std","A-hist"]; coh = ["internal","external"];
AE = rd(f.aEst); BE = rd(f.bEst);
cap = strings(0,1);

%% ---------------- S1 predictors and preprocessing ----------------------
MS = rd(f.l4miss);
S1a = table();
for p = unique(MS.pipeline, "stable")'
    q = MS(MS.pipeline == p & MS.type ~= "one-hot", :);
    o = MS(MS.pipeline == p & MS.type == "one-hot", :);
    S1a = [S1a; table(p, height(q), sum(q.train_NaN_or_missing + q.train_Inf > 0), sum(q.validation_NaN_or_missing + q.validation_Inf > 0), ...
        sum(contains(q.reason,"non-finite")), sum(contains(q.reason,"zero variance") & ~contains(q.reason,"near-zero")), sum(contains(q.reason,"near-zero")), ...
        sum(contains([q.reason; o.reason], "|r|")), sum(q.decision=="retained") + sum(o.decision=="retained"), ...
        'VariableNames', {'pipeline','candidate_source_predictors','with_NaN_or_Inf_training','with_NaN_or_Inf_validation', ...
        'removed_nonfinite','removed_zero_variance','removed_near_zero_variance','removed_correlation','final_predictors'})]; %#ok<AGROW>
end
wsheet(S1a, xl, "S1a_preprocessing_summary", outDir);
wsheet(MS(:, {'pipeline','feature','type','decision','reason','train_NaN_or_missing','validation_NaN_or_missing','validation_values_imputed_with_training_median'}), xl, "S1b_predictor_decisions", outDir);
wsheet(rd(f.l4ovl), xl, "S1c_predictor_overlap", outDir);
cap(end+1) = "Supplementary Table S1. Training-only preprocessing for each pipeline: (a) summary of exclusions; (b) decision and reason for every candidate predictor, including one-hot columns; (c) number of final predictors shared between pipelines.";

%% ---------------- S2 thresholds -----------------------------------------
S2 = rd(f.l4sum);
S2 = S2(:, {'pipeline','n_train','n_train_pos','n_validation','n_validation_pos','n_final_predictors','validation_auc', ...
    'threshold_raw','threshold_doublelogit','fallback_used','val_sensitivity','val_specificity','val_precision','val_MCC'});
wsheet(S2, xl, "S2_thresholds_validation", outDir);
cap(end+1) = "Supplementary Table S2. Model development: training and validation sizes, final predictors, validation AUC on each pipeline's own validation candidates, selected threshold on the raw margin and its doublelogit probability, and validation metrics at the selected threshold. Thresholds were selected on a grid of raw-margin cut-offs s = ln(p/(1 - p)), p = 0.001-0.999 (Amendment 02).";

%% ---------------- S3 full candidate-level estimates -----------------------
wsheet(fmtEst(AE), xl, "S3_candidate_level_estimates", outDir);
cap(end+1) = "Supplementary Table S3. All candidate-level estimates on identical Pipeline B candidates with 95% patient-clustered bootstrap CIs (2,000 resamples): AUC under three label conventions, paired differences, locked-threshold metrics, Brier score and acceptance by occupancy class, and the descriptive A-hist -> A-hist reconstruction. Sequential contrasts are descriptive, not causal.";

%% ---------------- S4 burden ---------------------------------------------
tok = regexp(BE.metric, "^(bias|LoAlow|LoAhigh|MAE|RMSE|Pearson|Spearman|CCC|n)_(volume|TLPSMA|SUVmean)_(B|Astd|Ahist|Ahist_ablation)$", "tokens", "once");
k = ~cellfun(@isempty, tok);
S4 = fmtEst(BE(k,:));
tk = vertcat(tok{k});
S4.statistic = string(tk(:,1)); S4.quantity = string(tk(:,2)); S4.model = string(tk(:,3));
S4 = movevars(S4, {'quantity','model','statistic'}, "After", "cohort");
wsheet(S4, xl, "S4_burden_agreement", outDir);
cap(end+1) = "Supplementary Table S4. Agreement between predicted and reference tumour burden (all studies, zeros included): total tumour volume, SUVmean (studies with both masks non-empty; n given) and TL-PSMA. Bias = predicted - reference; limits of agreement = bias +/- 1.96 SD.";

%% ---------------- S5 ablation -------------------------------------------
k = contains(BE.metric, "Ahist_ablation") | startsWith(BE.metric, ["cand","d_cand"]);
wsheet(fmtEst(BE(k,:)), xl, "S5_ablation", outDir);
cap(end+1) = "Supplementary Table S5. Supplementary ablation: A-hist retrained without PET first-order minimum, 10th percentile and range, evaluated on Pipeline B candidates (candidate level, segmentation, burden and detection) with paired differences against A-hist and B. The model was trained independently; results support but do not establish a causal role of the removed predictors.";

%% ---------------- S6 reconstruction -------------------------------------
S6 = [fmtEst(AE(contains(AE.metric,"Ahist_on_Ahist"),:)); fmtEst(BE(endsWith(BE.metric,"_Ahist_recon"),:))];
wsheet(S6, xl, "S6_reconstruction", outDir);
cap(end+1) = "Supplementary Table S6. Descriptive reference-informed reconstruction (A-hist applied to A-hist candidates). Internally the candidates are derived from the reference mask, so the near-perfect agreement is circular and is not deployable performance.";

%% ---------------- S7 calibration ----------------------------------------
wsheet(rd(f.aCal), xl, "S7a_calibration_summary", outDir);
wsheet(rd(f.aBins), xl, "S7b_calibration_deciles", outDir);
cap(end+1) = "Supplementary Table S7. Calibration of doublelogit probabilities P = 1/(1 + exp(-2s)) on Pipeline B candidates: Brier score, calibration intercept (calibration-in-the-large) and slope (point estimates), and decile bins.";

%% ---------------- S8 deviations and not performed ------------------------
S8 = cell2table({ ...
 "Deviation","Amendment 01 (9 Nov 2026, before any model was trained)","Preprocessing cut-offs written as implemented (<= 0.01, >= 95, |r| >= 0.95); predictor list clarified; relationship to the earlier pipeline clarified","Wording; no analysis change";
 "Deviation","Amendment 02 (11 Nov 2026, after freezing)","Threshold grid described as raw-margin cut-offs s = ln(p/(1 - p)); tumour-presence label and one-voxel-overlap recall terminology; importance definition","Wording; no analysis change";
 "Deviation","LP04 run 20261109_204348 superseded by 20261109_204900","Rerun with added leakage guards and preprocessing-invariance checks","Identical results";
 "Deviation","LP05a check S4b added before the single LP05a execution","Reproduction of all 999 validation threshold-grid rows and the selected thresholds","Quality control only";
 "Deviation","Post hoc diagnostics LP05a-D and LP05a-D2","Read-only mechanism diagnostics; D2 band edges set after viewing summary distributions, before band-specific acceptance; stopping rule recorded in advance","Exploratory; reported as such";
 "Deviation","LP05b attempt 1 (20261109_220050) aborted","Integrity check stopped the run: voxel volume taken from the reference-mask header; corrected to the PET grid with stricter checks (revision 1)","No results produced by attempt 1";
 "Not performed","Single-feature AUCs of construction-sensitive features (Methods 10.3)","Not performed; mechanism investigation completed by LP05a-D/D2 and ended by the prespecified stopping rule","Partial coverage: model-based importance (LP05a-D)";
 "Not performed","Fraction of voxels below SUV 3 by pipeline and class","Not performed (same reason)","Partial coverage: SUVmin < 3 share (audit) and SUVmin distributions (LP05a-D)";
 "Not performed","Surface-to-volume ratio and compactness by pipeline and class","Not performed (same reason)","Partial coverage: selected standardised mean differences (LP05a-D)";
 "Partly performed","Negative-construction audit","Retention counts and per-candidate 5-mm flags reported; volume and geometry comparison of A-hist negatives not performed","Counts in LP03 and pre-analysis audit"}, ...
 'VariableNames', {'category','item','description','consequence'});
wsheet(S8, xl, "S8_deviations_not_performed", outDir);
cap(end+1) = "Supplementary Table S8. Deviations from the locked analysis plan and prespecified analyses not performed. Items not performed are listed explicitly; no results are implied for them.";

%% ---------------- S9 frozen runs and safeguards --------------------------
S9a = cell2table({ ...
 "Locked methods v1.0","LeakagePaper_Methods_LOCKED_v1.0_2026-11-09.md","d9ca33362c8459f29bd8d870e5df55c1e63f1d8cfcfe88edf427e4b3695c0fe4";
 "Amendment 01","LeakagePaper_Methods_AMENDMENT_01_2026-11-09.md","cba53eab04c465e8f205cc15b2cc257656770387ef2197250cf63f2fe76b747d";
 "Amendment 02","LeakagePaper_Methods_AMENDMENT_02_2026-11-11.md","8eb4c1667800aa437f6b78cd8744abe98f2bbd623616e5f347664b7db08f605c";
 "LP03 candidate sets (official)","LP03_20261109_203442/LP03_candidate_sets.mat","6ce98177fb7dffda4804d1a0717849be331571dfe252b842514b0f15583f3eb2";
 "LP04 models and thresholds (official)","LP04_20261109_204900/LP04_models_and_thresholds.mat","bc87e01862d59086d3134aa73e61382276aaec7bce4dcd584e109da5a209a4d6";
 "LP05a results (official)","LP05a_20261109_210731/LP05a_results.mat","caa6b7cc06a2fb68f158e98ea1df3c61b0276e5add39cd30f2ca17e6720f1fd3";
 "LP05a-D results","LP05aD_20261109_212046/LP05aD_results.mat","98dd3a83ed79a4a9db3a56a3e411c4f822cdaf45ff2832ae4c6463c3cccf8eff";
 "LP05a-D2 results","LP05aD2_20261109_214757/LP05aD2_results.mat","c6e50567f5669f877887e8a62c0c6d1761da0dbd5164e616a4acff8148b094c5";
 "LP05a-D2 pre-run commitment","LP05aD2_PRE_RUN_COMMITMENT.md","3142c748f907a28988c5587cac94f5ef47184416c6c7ccdc586d731064ec9c0c";
 "LP05b results (official)","LP05b_20261109_220537/LP05b_results.mat","d5eb35988f7885c06f85ca3c5e6c6ea6e77c681b087cf169f4bc2fa232164374";
 "Internal test candidate ID list","LP05a S1 (2,673 IDs)","0ac8b41344f89df2dcddada653b23e3e5aec1d79d34e11393b0924862b01d697";
 "External candidate ID list (locked 95)","LP05a S1 (4,173 IDs)","4a2761baea90644a49c52e5a5293939bb0c65ead6124a5ed3ed834f6469d4638"}, ...
 'VariableNames', {'item','file','SHA256'});
wsheet(S9a, xl, "S9a_frozen_runs_hashes", outDir);
% read by column position (header detection of these short CSVs is unreliable)
sa = readtable(f.aSafe, "TextType","string", "Delimiter",",", "NumHeaderLines",1, "ReadVariableNames",false);
sb = readtable(f.bSafe, "TextType","string", "Delimiter",",", "NumHeaderLines",1, "ReadVariableNames",false);
S9b = [table(repmat("LP05a",height(sa),1), string(sa{:,1}) + ": " + string(sa{:,2}) + " -> " + string(sa{:,3}), 'VariableNames', {'stage','result'}); ...
       table(repmat("LP05b",height(sb),1), string(sb{:,1}), 'VariableNames', {'stage','result'})];
wsheet(S9b, xl, "S9b_safeguards", outDir);
cap(end+1) = "Supplementary Table S9. Frozen analysis runs with SHA256 hashes, and the safeguards verified in the evaluation scripts before results were produced.";

%% ---------------- S10 Figure 4 source data -------------------------------
wsheet(rd(f.aAcc), xl, "S10a_acceptance_by_occupancy", outDir);
FD = rd(f.bFail); C1 = rd(f.bC1);
wsheet(join(FD, C1, "Keys", {'cohort','model'}), xl, "S10b_lesion_outcomes", outDir);
cap(end+1) = "Supplementary Table S10. Source data for Figure 4: acceptance, accepted volume and reference-tumour / non-tumour volume by occupancy class for each model at its locked threshold; and outcomes of reference lesions (one-voxel-overlap detection, classifier and candidate-generation misses, detection exclusively through accepted C1 candidates).";

%% ================= FIGURES ================================================
sty = struct("B", {{"-", [0 0 0], "o", [0 0 0], 1.6}}, "Astd", {{"--", [0.45 0.45 0.45], "s", [0.75 0.75 0.75], 1.4}}, ...
             "Ahist", {{":", [0 0 0], "^", [1 1 1], 1.8}}, "Ahist_ablation", {{"-.", [0.3 0.3 0.3], "d", [0.5 0.5 0.5], 1.2}});
ge = @(E,c,m) E(E.cohort==c & E.metric==m, :);

% ---- S-Fig 1 ROC ---------------------------------------------------------
fig = figure("Color","w","Units","centimeters","Position",[2 2 24 11]);
t = tiledlayout(fig,1,2,"TileSpacing","compact","Padding","compact");
for ci = 1:2
    c = coh(ci); sfx = ["I","E"]; S = rd(f.("aSc" + sfx(ci)));
    ax = nexttile(t); hold(ax,"on"); hh = gobjects(0); lab = strings(0);
    for j = 1:3
        m = mods(j); s = sty.(m); y = double(S.label_any);
        [X, Y] = perfcurve(y, double(S.("s_" + m)), 1);
        a = ge(AE, c, "AUC_" + m + "_any");
        hh(end+1) = plot(ax, X, Y, "LineStyle", s{1}, "Color", s{2}, "LineWidth", s{5}); %#ok<AGROW>
        lab(end+1) = sprintf("%s: AUC %.3f (%.3f-%.3f)", modLab(j), a.estimate, a.ci_low, a.ci_high); %#ok<AGROW>
        acc = logical(S.("accept_" + m));
        plot(ax, mean(acc(y==0)), mean(acc(y==1)), s{3}, "MarkerSize", 8, "MarkerEdgeColor", "k", "MarkerFaceColor", s{4}, "HandleVisibility","off");
    end
    plot(ax, [0 1], [0 1], "Color", [0.7 0.7 0.7], "LineStyle", "-", "HandleVisibility","off");
    axis(ax, "square"); xlim(ax,[0 1]); ylim(ax,[0 1]); grid(ax,"on");
    xlabel(ax, "1 - specificity"); ylabel(ax, "Sensitivity");
    title(ax, sprintf("(%s) %s Pipeline B candidates (n = %s)", char('a'+ci-1), c, fmtN(height(S))));
    legend(ax, hh, lab, "Location","southeast", "FontSize",7);
end
saveFig(fig, outDir, "SFig1_ROC_locked_operating_points");
cap(end+1) = "Supplementary Figure S1. ROC curves of the three models on identical Pipeline B candidates (any-overlap labels). Markers show each model's locked source threshold (circle B, square A-std, triangle A-hist). AUCs with 95% patient-clustered CIs.";

% ---- S-Fig 2 calibration --------------------------------------------------
BI = rd(f.aBins); AC = rd(f.aCal);
fig = figure("Color","w","Units","centimeters","Position",[2 2 24 11]);
t = tiledlayout(fig,1,2,"TileSpacing","compact","Padding","compact");
for ci = 1:2
    c = coh(ci); ax = nexttile(t); hold(ax,"on"); hh = gobjects(0); lab = strings(0);
    plot(ax, [0 1], [0 1], "Color", [0.7 0.7 0.7], "HandleVisibility","off");
    for j = 1:3
        m = mods(j); s = sty.(m); b = BI(BI.cohort==c & BI.model==m, :); a = AC(AC.cohort==c & AC.model==m, :);
        hh(end+1) = plot(ax, b.mean_predicted, b.observed_rate, "LineStyle", s{1}, "Color", s{2}, "LineWidth", s{5}, ...
            "Marker", s{3}, "MarkerEdgeColor","k", "MarkerFaceColor", s{4}); %#ok<AGROW>
        lab(end+1) = sprintf("%s: Brier %.3f, intercept %.2f, slope %.2f", modLab(j), a.Brier, a.calibration_intercept_CITL, a.calibration_slope); %#ok<AGROW>
    end
    axis(ax,"square"); xlim(ax,[0 1]); ylim(ax,[0 1]); grid(ax,"on");
    xlabel(ax, "Mean predicted probability (doublelogit)"); ylabel(ax, "Observed proportion positive");
    title(ax, sprintf("(%s) %s", char('a'+ci-1), c)); legend(ax, hh, lab, "Location","northwest", "FontSize",7);
end
saveFig(fig, outDir, "SFig2_calibration");
cap(end+1) = "Supplementary Figure S2. Calibration of doublelogit probabilities by decile of predicted probability. The reference-informed models show marked miscalibration with systematic underprediction and slopes below 1.";

% ---- S-Fig 3 volume agreement -------------------------------------------
ST = rd(f.bStudy);
fig = figure("Color","w","Units","centimeters","Position",[2 2 32 16]);
t = tiledlayout(fig,2,4,"TileSpacing","compact","Padding","compact");
for ci = 1:2
    c = coh(ci);
    ax = nexttile(t); hold(ax,"on"); hh = gobjects(0);
    for j = 1:3
        m = mods(j); s = sty.(m); q = ST(ST.cohort==c & ST.model==m, :);
        hh(end+1) = scatter(ax, q.ref_volume_mL + 1, q.pred_volume_mL + 1, 22, "Marker", s{3}, "MarkerEdgeColor","k", "MarkerFaceColor", s{4}); %#ok<AGROW>
    end
    lim = [1 2e4]; plot(ax, lim, lim, "k--", "HandleVisibility","off");
    set(ax, "XScale","log", "YScale","log"); xlim(ax, lim); ylim(ax, lim); axis(ax,"square"); grid(ax,"on");
    xlabel(ax, "Reference volume + 1 (mL)"); ylabel(ax, "Predicted volume + 1 (mL)");
    title(ax, sprintf("(%s) %s: predicted vs reference", char('a'+4*(ci-1)), c));
    legend(ax, hh, modLab, "Location","northwest", "FontSize",7);
    for j = 1:3
        m = mods(j); s = sty.(m); q = ST(ST.cohort==c & ST.model==m, :);
        ax = nexttile(t); hold(ax,"on");
        mn = (q.ref_volume_mL + q.pred_volume_mL)/2; df = q.pred_volume_mL - q.ref_volume_mL;
        scatter(ax, mn + 1, df, 18, "Marker", s{3}, "MarkerEdgeColor","k", "MarkerFaceColor", s{4});
        bi = ge(BE,c,"bias_volume_"+m); lo = ge(BE,c,"LoAlow_volume_"+m); hi = ge(BE,c,"LoAhigh_volume_"+m);
        yline(ax, bi.estimate, "k-", "LineWidth",1.2); yline(ax, lo.estimate, "k--"); yline(ax, hi.estimate, "k--"); yline(ax, 0, "Color",[0.7 0.7 0.7]);
        set(ax, "XScale","log"); grid(ax,"on");
        xlabel(ax, "Mean of predicted and reference + 1 (mL)"); ylabel(ax, "Predicted - reference (mL)");
        title(ax, sprintf("(%s) %s, %s\nbias %+.0f mL; LoA %+.0f to %+.0f", char('a'+4*(ci-1)+j), c, modLab(j), bi.estimate, lo.estimate, hi.estimate), "FontSize",8);
    end
end
saveFig(fig, outDir, "SFig3_volume_agreement");
cap(end+1) = "Supplementary Figure S3. Agreement of total tumour volume per study. Left: predicted versus reference volume (log scales; 1 mL added to show zero values). Right: Bland-Altman plots for each model (solid line, bias; dashed lines, 95% limits of agreement, frozen patient-clustered estimates).";

% ---- S-Fig 4 Dice per study ---------------------------------------------
fig = figure("Color","w","Units","centimeters","Position",[2 2 24 10]);
t = tiledlayout(fig,1,2,"TileSpacing","compact","Padding","compact");
m4 = ["B","Astd","Ahist","Ahist_ablation"]; l4 = ["B","A-std","A-hist","Ablation"];
for ci = 1:2
    c = coh(ci); ax = nexttile(t); hold(ax,"on");
    q = ST(ST.cohort==c & ST.tumour_positive==1 & ismember(ST.model, m4), :);
    g = categorical(q.model, m4, l4);
    boxchart(ax, g, q.Dice, "BoxFaceColor",[0.75 0.75 0.75], "BoxEdgeColor","k", "WhiskerLineColor","k", "MarkerStyle","none");
    rng(1); swarmchart(ax, g, q.Dice, 8, "k", "filled", "XJitterWidth", 0.4);
    ylim(ax,[0 1]); ylabel(ax, "Dice per tumour-positive study"); grid(ax,"on");
    title(ax, sprintf("(%s) %s (n = %d studies)", char('a'+ci-1), c, numel(unique(q.case_id))));
end
saveFig(fig, outDir, "SFig4_Dice_per_study");
cap(end+1) = "Supplementary Figure S4. Dice per tumour-positive study for each model at its locked threshold (box: median and interquartile range; points: individual studies; horizontal jitter for display only).";

% ---- S-Fig 5 reconstruction ---------------------------------------------
fig = figure("Color","w","Units","centimeters","Position",[2 2 24 11]);
t = tiledlayout(fig,1,2,"TileSpacing","compact","Padding","compact");
for ci = 1:2
    c = coh(ci); ax = nexttile(t); hold(ax,"on");
    q = ST(ST.cohort==c & ST.model=="Ahist_recon", :);
    scatter(ax, q.ref_volume_mL + 1, q.pred_volume_mL + 1, 22, "MarkerEdgeColor","k", "MarkerFaceColor",[0.6 0.6 0.6]);
    lim = [1 2e4]; plot(ax, lim, lim, "k--");
    set(ax,"XScale","log","YScale","log"); xlim(ax,lim); ylim(ax,lim); axis(ax,"square"); grid(ax,"on");
    cc = ge(BE,c,"CCC_volume_Ahist_recon"); dc = ge(BE,c,"Dice_Ahist_recon");
    xlabel(ax, "Reference volume + 1 (mL)"); ylabel(ax, "Reconstructed volume + 1 (mL)");
    ttl = ["internal: circular (candidates derived from the reference)", "external"];
    title(ax, sprintf("(%s) %s\nCCC %.3f (%.3f-%.3f); Dice %.2f", char('a'+ci-1), ttl(ci), cc.estimate, cc.ci_low, cc.ci_high, dc.estimate), "FontSize",9);
end
saveFig(fig, outDir, "SFig5_reconstruction_internal_external");
cap(end+1) = "Supplementary Figure S5. Descriptive reference-informed reconstruction (A-hist applied to A-hist candidates). Internally, near-perfect agreement is circular because the candidates are derived from the reference mask; externally, agreement collapses. Not deployable performance.";

% ---- S-Fig 6 threshold curves -------------------------------------------
fig = figure("Color","w","Units","centimeters","Position",[2 2 24 18]);
t = tiledlayout(fig,2,2,"TileSpacing","compact","Padding","compact");
ml = ["B (PET-only)","A-std","A-hist","A-hist ablation"];
for j = 1:4
    m = m4(j); G = rd(gridF.(m)); ax = nexttile(t); hold(ax,"on");
    h1 = plot(ax, G.threshold_raw, G.sensitivity, "k-", "LineWidth",1.4);
    h2 = plot(ax, G.threshold_raw, G.specificity, "--", "Color",[0.45 0.45 0.45], "LineWidth",1.4);
    h3 = plot(ax, G.threshold_raw, G.MCC, "k:", "LineWidth",1.6);
    th = G.threshold_raw(logical(G.selected)); xline(ax, th, "k-.", "LineWidth",1);
    yline(ax, 0.90, "Color",[0.7 0.7 0.7]);
    ylim(ax,[0 1.02]); xlim(ax,[-7 7]); grid(ax,"on");
    xlabel(ax, "Threshold on raw ensemble margin s"); ylabel(ax, "Validation metric");
    title(ax, sprintf("(%s) %s: selected s = %.3f", char('a'+j-1), ml(j), th));
    if j == 1, legend(ax, [h1 h2 h3], ["Sensitivity","Specificity","MCC"], "Location","southwest", "FontSize",7); end
end
saveFig(fig, outDir, "SFig6_validation_threshold_curves");
cap(end+1) = "Supplementary Figure S6. Validation threshold curves on each pipeline's own validation candidates over the raw-margin grid (999 cut-offs). Dash-dot line: selected threshold (sensitivity >= 0.90, then maximum MCC); grey line: sensitivity 0.90.";

writelines(cap, fullfile(outDir, "LP06_supp_captions.txt"));
fprintf("\nLP06 supplementary items COMPLETED. Frozen outputs only; no new statistical analysis.\nOutput: %s\n", outDir);
end

%% ======================================================================
function T = fmtEst(E)
    s = strings(height(E),1); p = strings(height(E),1);
    for i = 1:height(E)
        if isfinite(E.estimate(i))
            s(i) = sprintf("%.4g (%.4g to %.4g)", E.estimate(i), E.ci_low(i), E.ci_high(i));
        end
        if isfinite(E.bootstrap_p_two_sided(i))
            if E.bootstrap_p_two_sided(i) < 0.001, p(i) = "< 0.001"; else, p(i) = sprintf("%.3f", E.bootstrap_p_two_sided(i)); end
        end
    end
    T = table(E.cohort, E.metric, E.family, E.estimate, E.ci_low, E.ci_high, s, p, ...
        'VariableNames', {'cohort','metric','family','estimate','ci_low','ci_high','estimate_95CI','bootstrap_p_two_sided'});
end

function wsheet(T, xl, sheet, outDir)
    writetable(T, xl, "Sheet", sheet);
    writetable(T, fullfile(outDir, "SuppTable_" + sheet + ".csv"));
    fprintf("  %-32s %5d rows\n", sheet, height(T));
end

function s = fmtN(x)
    s = string(sprintf('%d', x));
    if x >= 1000, s = string(sprintf('%d,%03d', floor(x/1000), mod(x,1000))); end
end

function saveFig(fig, outDir, name)
    exportgraphics(fig, fullfile(outDir, name + ".png"), "Resolution", 600);
    exportgraphics(fig, fullfile(outDir, name + ".tif"), "Resolution", 600);
    exportgraphics(fig, fullfile(outDir, name + ".pdf"), "ContentType", "vector");
    savefig(fig, fullfile(outDir, name + ".fig"));
    fprintf("  %s saved (PNG/TIFF 600 dpi, PDF, FIG)\n", name);
    close(fig);
end

function h = fileSha256(f)
    if ~isfile(f), h = "missing"; return; end
    md = java.security.MessageDigest.getInstance("SHA-256");
    fid = fopen(f, "r"); c = onCleanup(@() fclose(fid));
    while true
        buf = fread(fid, 2^24, "*uint8");
        if isempty(buf), break; end
        md.update(typecast(buf, "int8"));
    end
    d = md.digest();
    h = string(lower(reshape(dec2hex(typecast(int8(d), "uint8"), 2)', 1, [])));
end
