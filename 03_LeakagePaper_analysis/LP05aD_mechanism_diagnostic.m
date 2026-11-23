function LP05aD_mechanism_diagnostic()
%LP05AD_MECHANISM_DIAGNOSTIC  Leakage paper - READ-ONLY mechanism diagnostic.
%
% Purpose: examine whether the transfer failure of the reference-informed
% models (LP05a, frozen) is consistent with a construction-specific shortcut
% (minimum SUV / 10th-percentile SUV) and screen for other distribution shifts.
% DIAGNOSTIC AND DESCRIPTIVE ONLY: no retraining, no preprocessing change,
% no threshold change, LP05a outputs are not modified. The ablation model is
% the already-trained LP04 model; its formal supplementary estimates (with CI)
% remain part of LP05b.
%
% Parts
%   0  input integrity: LP03/LP04 SHA256; scores reproduce frozen LP05a scores
%   1  feature distributions (SUVmin, SUV 10th percentile, SUVmean, volume) in
%      training, internal-test and external populations
%   2  raw-score distributions and acceptance per model and population
%   3  KEY TEST: acceptance of reference-derived and PET-only positives,
%      stratified by SUVmin < 3 vs >= 3, for every model incl. ablation
%   4  ablation model on Pipeline B candidates (point estimates; discrimination
%      and threshold/calibration reported separately)
%   5  predictor importance of A-hist and ablation models
%   6  distribution-shift screen: standardized mean differences of all A-hist
%      model inputs between training and transfer populations
% Interpretation rule: results can support or weaken the SUVmin hypothesis;
% they cannot by themselves establish causality.

clc;
t0 = tic;
V2 = string(fullfile(psma_root("PSMA_PROJECT_ROOT"),"PSMA_AI_Project_V2_PET_Derived"));
MD = fullfile(V2,"06_results","LeakagePaper_ModelData");
EV = fullfile(V2,"06_results","LeakagePaper_Evaluation");
lp03Mat = fullfile(MD,"LP03_20261109_203442","LP03_candidate_sets.mat");
lp04Mat = fullfile(MD,"LP04_20261109_204900","LP04_models_and_thresholds.mat");
lp05aDir = fullfile(EV,"LP05a_20261109_210731");
sha03 = "6ce98177fb7dffda4804d1a0717849be331571dfe252b842514b0f15583f3eb2";
sha04 = "bc87e01862d59086d3134aa73e61382276aaec7bce4dcd584e109da5a209a4d6";

outDir = fullfile(EV, "LP05aD_" + string(datetime("now","Format","yyyyMMdd_HHmmss")));
mkdir(outDir);
diary(char(fullfile(outDir,"LP05aD_log.txt"))); cleanupDiary = onCleanup(@() diary("off"));
fprintf("LP05a-D (read-only diagnostic) started %s\nMATLAB %s\nOutput: %s\n\n", string(datetime("now")), version, outDir);

%% Part 0 integrity -----------------------------------------------------
if fileSha256(lp03Mat) ~= sha03, error("LP03 hash mismatch."); end
if fileSha256(lp04Mat) ~= sha04, error("LP04 hash mismatch."); end
fprintf("Part 0: LP03 and LP04 hashes verified.\n");
D = load(lp03Mat, "data"); data = D.data;
M = load(lp04Mat, "models"); models = M.models;
mods = ["B","Astd","Ahist","Ahist_ablation"];
FMIN = "PET_original_firstorder_Minimum"; F10 = "PET_original_firstorder_10Percentile";

% populations
P = struct();
P.Ahist_train_pos = sel(data.internal_Ahist, "train", 1);
P.Ahist_train_neg = sel(data.internal_Ahist, "train", 0);
P.Astd_train_neg  = sel(data.internal_Astd,  "train", 0);
P.B_train_pos     = sel(data.internal_B, "train", 1);
P.B_train_neg     = sel(data.internal_B, "train", 0);
P.B_test_pos      = sel(data.internal_B, "test", 1);
P.B_test_neg      = sel(data.internal_B, "test", 0);
P.Ahist_test_pos  = sel(data.internal_Ahist, "test", 1);
P.Ahist_test_neg  = sel(data.internal_Ahist, "test", 0);
P.ext_A_pos       = data.external_Ahist(data.external_Ahist.y==1, :);
P.ext_A_neg       = data.external_Ahist(data.external_Ahist.y==0, :);
P.ext_B_pos       = data.external_B(data.external_B.y==1, :);
P.ext_B_neg       = data.external_B(data.external_B.y==0, :);
pn = string(fieldnames(P))';

% scores for every model x population
SC = struct();
for p = pn
    for m = mods
        SC.(p).(m) = scoreTable(P.(p), models.(m));
    end
end

% consistency with frozen LP05a scores
L5i = readtable(fullfile(lp05aDir,"LP05a_scores_internal_test_B.csv"), "VariableNamingRule","preserve");
L5e = readtable(fullfile(lp05aDir,"LP05a_scores_external_B.csv"), "VariableNamingRule","preserve");
maxd = 0;
for m = ["B","Astd","Ahist"]
    maxd = max(maxd, cmpScores(data.internal_B(data.internal_B.split=="test",:), models.(m), L5i, "s_" + m));
    maxd = max(maxd, cmpScores(data.external_B, models.(m), L5e, "s_" + m));
end
if maxd > 1e-9, error("Scores differ from frozen LP05a (max %.3g).", maxd); end
fprintf("Part 0: scores reproduce frozen LP05a scores (max abs diff %.2g).\n\n", maxd);

%% Part 1 feature distributions ---------------------------------------
fd = table();
for p = pn
    T = P.(p);
    fd = [fd; table(p, height(T), ...
        median(T.(FMIN)), q(T.(FMIN),.25), q(T.(FMIN),.75), 100*mean(T.(FMIN) < 3), ...
        median(T.(F10)), q(T.(F10),.25), q(T.(F10),.75), 100*mean(T.(F10) < 3), ...
        median(T.SUVmean), median(T.candidate_volume_mL), ...
        'VariableNames', {'population','n','SUVmin_median','SUVmin_Q1','SUVmin_Q3','pct_SUVmin_lt3', ...
        'SUV10_median','SUV10_Q1','SUV10_Q3','pct_SUV10_lt3','SUVmean_median','volume_mL_median'})]; %#ok<AGROW>
end
fprintf("Part 1: feature distributions\n"); disp(fd);

%% Part 2 score distributions and acceptance --------------------------
sd = table();
for p = pn
    for m = mods
        s = SC.(p).(m); thr = models.(m).threshold_raw;
        sd = [sd; table(p, m, numel(s), median(s), q(s,.25), q(s,.75), thr, 100*mean(s >= thr), ...
            'VariableNames', {'population','model','n','score_median','score_Q1','score_Q3','threshold_raw','pct_accepted'})]; %#ok<AGROW>
    end
end
fprintf("Part 2: raw-score distributions and acceptance (pct_accepted = %% at or above the model's locked threshold)\n");
disp(sd(ismember(sd.population, ["B_test_pos","B_test_neg","Ahist_test_pos","ext_A_pos","ext_B_pos","ext_B_neg"]), :));

%% Part 3 KEY TEST: acceptance stratified by SUVmin ---------------------
st = table();
for p = ["Ahist_train_pos","Ahist_test_pos","ext_A_pos","B_test_pos","ext_B_pos"]
    T = P.(p);
    for stratum = ["SUVmin<3","SUVmin>=3"]
        if stratum == "SUVmin<3", k = T.(FMIN) < 3; else, k = T.(FMIN) >= 3; end
        for m = mods
            s = SC.(p).(m)(k); thr = models.(m).threshold_raw;
            st = [st; table(p, stratum, m, sum(k), 100*mean(s >= thr), median(s), ...
                'VariableNames', {'positive_population','stratum','model','n','pct_accepted','score_median'})]; %#ok<AGROW>
        end
    end
end
fprintf("Part 3 (KEY): acceptance of POSITIVE candidates by SUVmin stratum\n");
fprintf("  Hypothesis-consistent pattern: A-hist accepts SUVmin<3 positives but rejects SUVmin>=3 positives,\n");
fprintf("  and the gap narrows in the ablation model. Other patterns weaken the hypothesis.\n");
disp(st);

%% Part 4 ablation on Pipeline B candidates ----------------------------
ab = table();
cohortT = struct("internal_test", data.internal_B(data.internal_B.split=="test",:), "external", data.external_B);
for c = ["internal_test","external"]
    T = cohortT.(c); y = double(T.label_any);
    for m = ["B","Ahist","Ahist_ablation"]
        s = scoreTable(T, models.(m)); thr = models.(m).threshold_raw;
        yh = s >= thr; Pp = 1./(1+exp(-2*s));
        TP = sum(yh & y==1); FP = sum(yh & y==0); TN = sum(~yh & y==0); FN = sum(~yh & y==1);
        b0 = glmfit(ones(size(y)), y, "binomial", "Constant","off", "Offset", 2*s);
        ab = [ab; table(c, m, aucMW(y,s), TP/(TP+FN), TN/(TN+FP), mean((Pp-y).^2), b0(1), mean(Pp), mean(y), ...
            'VariableNames', {'cohort','model','AUC_discrimination','sensitivity_at_locked_threshold', ...
            'specificity_at_locked_threshold','Brier','calibration_intercept_CITL','mean_predicted','observed_rate'})]; %#ok<AGROW>
    end
end
fprintf("Part 4: ablation on Pipeline B candidates (diagnostic point estimates; formal CIs in LP05b)\n");
fprintf("  Discrimination (AUC) and threshold/calibration (sensitivity, CITL) are reported separately.\n");
disp(ab);

%% Part 5 predictor importance ------------------------------------------
imp = table();
for m = ["Ahist","Ahist_ablation","B"]
    mdl = models.(m).model;
    v = predictorImportance(mdl);
    names = models.(m).predictors(:);
    [vs, o] = sort(v(:), "descend");
    r = (1:numel(vs))';
    imp = [imp; table(repmat(m,numel(vs),1), r, names(o), vs, 100*vs/sum(vs), ...
        'VariableNames', {'model','rank','predictor','importance','pct_of_total'})]; %#ok<AGROW>
end
fprintf("Part 5: predictor importance (top 15 per model; rank of SUVmin / SUV10 shown)\n");
for m = ["Ahist","Ahist_ablation","B"]
    I = imp(imp.model == m, :);
    disp(I(1:min(15,height(I)), :));
    for f = [FMIN, F10]
        rr = I.rank(I.predictor == f);
        if isempty(rr), fprintf("  %s: %s not a model input\n", m, f);
        else, fprintf("  %s: %s rank %d of %d (%.1f%% of total importance)\n", m, f, rr, height(I), I.pct_of_total(I.predictor==f)); end
    end
end

%% Part 6 distribution-shift screen ------------------------------------
pairs = ["Ahist_train_pos","ext_A_pos"; "Ahist_train_pos","B_test_pos"; "Ahist_train_pos","ext_B_pos"; ...
         "Ahist_train_neg","B_test_neg"; "Ahist_train_neg","ext_B_neg"];
prep = models.Ahist.preprocessing;
Iah = imp(imp.model=="Ahist", :);
sh = table();
for k = 1:size(pairs,1)
    X1 = applyPreprocessing(P.(pairs(k,1)), prep); X2 = applyPreprocessing(P.(pairs(k,2)), prep);
    smd = (mean(X1) - mean(X2)) ./ sqrt((var(X1) + var(X2))/2);
    smd(~isfinite(smd)) = 0;
    [~,o] = sort(abs(smd), "descend");
    for j = 1:numel(o)
        f = prep.finalNames(o(j));
        sh = [sh; table(pairs(k,1) + " vs " + pairs(k,2), j, f, smd(o(j)), Iah.rank(Iah.predictor==f), ...
            'VariableNames', {'comparison','shift_rank','predictor','SMD','Ahist_importance_rank'})]; %#ok<AGROW>
    end
end
fprintf("Part 6: largest standardized mean differences in A-hist model inputs (top 10 per comparison)\n");
for k = 1:size(pairs,1)
    cmp = pairs(k,1) + " vs " + pairs(k,2);
    S6 = sh(sh.comparison == cmp, :);
    disp(S6(1:10, :));
    for f = [FMIN, F10]
        rr = S6.shift_rank(S6.predictor == f);
        if ~isempty(rr), fprintf("  %s: %s shift rank %d, SMD %.2f\n", cmp, f, rr, S6.SMD(S6.predictor==f)); end
    end
    fprintf("  predictors with |SMD| >= 0.5: %d of %d\n\n", sum(abs(S6.SMD) >= 0.5), height(S6));
end

%% Save --------------------------------------------------------------------
writetable(fd, fullfile(outDir,"LP05aD_part1_feature_distributions.csv"));
writetable(sd, fullfile(outDir,"LP05aD_part2_score_distributions_acceptance.csv"));
writetable(st, fullfile(outDir,"LP05aD_part3_acceptance_by_SUVmin_stratum.csv"));
writetable(ab, fullfile(outDir,"LP05aD_part4_ablation_on_B_candidates.csv"));
writetable(imp, fullfile(outDir,"LP05aD_part5_predictor_importance.csv"));
writetable(sh, fullfile(outDir,"LP05aD_part6_distribution_shift_SMD.csv"));
save(fullfile(outDir,"LP05aD_results.mat"), "fd","sd","st","ab","imp","sh","-v7.3");
fprintf("LP05a-D COMPLETED (%.1f min). Read-only; no model, threshold or LP05a output changed.\nResults: %s\n", toc(t0)/60, outDir);
end

%% ======================================================================
function T = sel(T, split, y)
    T = T(T.split == split & T.y == y, :);
end

function v = q(x, p)
    v = quantile(x, p);
end

function d = cmpScores(T, mm, L5, col)
    s = scoreTable(T, mm);
    [tf, loc] = ismember(T.row_id, string(L5.row_id));
    if ~all(tf), error("Row IDs missing from frozen LP05a scores."); end
    d = max(abs(s - L5.(col)(loc)));
end

function a = aucMW(y, s)
    y = double(y(:)); s = double(s(:));
    n1 = sum(y==1); n0 = sum(y==0);
    if n1 == 0 || n0 == 0, a = NaN; return; end
    r = tiedrank(s);
    a = (sum(r(y==1)) - n1*(n1+1)/2) / (n1*n0);
end

function s = scoreTable(T, mm)
    X = applyPreprocessing(T, mm.preprocessing);
    [~, sc] = predict(mm.model, X);
    s = double(sc(:, mm.model.ClassNames == 1));
end

function X = applyPreprocessing(T, prep)
    Xnum = double(table2array(T(:, cellstr(prep.keepNumeric))));
    for j = 1:size(Xnum,2)
        bad = ~isfinite(Xnum(:,j));
        if any(bad), Xnum(bad,j) = prep.trainingMedians(j); end
    end
    Xenc = zeros(height(T), numel(prep.encodedSource));
    for k = 1:numel(prep.encodedSource)
        Xenc(:,k) = double(normalizeCategory(T.(prep.encodedSource(k))) == prep.encodedCategory(k));
    end
    Xall = [Xnum, Xenc];
    [tf, loc] = ismember(prep.finalNames, prep.allPreparedNames);
    if ~all(tf), error("Final predictor missing."); end
    X = Xall(:, loc);
end

function values = normalizeCategory(x)
    values = strtrim(string(x));
    miss = ismissing(values) | strlength(values) == 0 | lower(values) == "nan" | ...
        lower(values) == "na" | lower(values) == "none";
    values(miss) = "<missing>";
end

function h = fileSha256(f)
    md = java.security.MessageDigest.getInstance("SHA-256");
    fid = fopen(f, "r"); if fid < 0, error("Cannot open %s", f); end
    c = onCleanup(@() fclose(fid));
    while true
        buf = fread(fid, 2^24, "*uint8");
        if isempty(buf), break; end
        md.update(typecast(buf, "int8"));
    end
    d = md.digest();
    h = string(lower(reshape(dec2hex(typecast(int8(d), "uint8"), 2)', 1, [])));
end
