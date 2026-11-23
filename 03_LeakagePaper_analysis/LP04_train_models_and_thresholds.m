function LP04_train_models_and_thresholds(lp03Mat)
%LP04_TRAIN_MODELS_AND_THRESHOLDS  Leakage paper - preprocessing, models, thresholds.
%
% Locked methods: LeakagePaper_Methods_LOCKED_v1.0_2026-11-09.md (Sections 6-8, 10.6)
%
% For each pipeline (B, Astd, Ahist) plus the supplementary ablation
% (Ahist_ablation = A-hist without PET_original_firstorder_Minimum,
% _10Percentile, _Range):
%   1. Training-only preprocessing, replicating V2 Objective 9G exactly:
%      non-finite (any) -> remove; zero variance -> remove; near-zero variance
%      (unique fraction <= 0.01 AND dominant/second frequency ratio >= 95)
%      -> remove; location zones one-hot encoded from training categories;
%      greedy correlation filter |r| >= 0.95 in column order (later column
%      removed); training medians lock imputation for unexpected gaps.
%   2. LogitBoost (fitcensemble), 200 cycles, LearnRate 0.05, MaxNumSplits 31,
%      MinLeafSize 5, Surrogate off, ClassNames [0 1], rng(20260717,'twister').
%   3. Threshold on the model's raw margin s, selected on that pipeline's own
%      internal validation candidates: grid s = logit(p), p = 0.001:0.001:0.999;
%      sensitivity >= 0.90 then max MCC; ties: precision, specificity,
%      balanced accuracy, higher threshold; fallback max MCC; undefined MCC
%      excluded (V2 9H2 REVISED_v2 rule).
%
% Test and external candidates are NOT scored here (that is LP05).
% Usage:  LP04_train_models_and_thresholds()            % newest LP03 output
%         LP04_train_models_and_thresholds("...\LP03_candidate_sets.mat")

clc;
t0 = tic;
V2 = string(fullfile(psma_root("PSMA_PROJECT_ROOT"),"PSMA_AI_Project_V2_PET_Derived"));
dataRoot = fullfile(V2,"06_results","LeakagePaper_ModelData");
if nargin < 1 || strlength(string(lp03Mat)) == 0
    d = dir(fullfile(dataRoot,"LP03_*","LP03_candidate_sets.mat"));
    if isempty(d), error("No LP03 output found under %s", dataRoot); end
    [~,k] = max([d.datenum]);
    lp03Mat = fullfile(d(k).folder, d(k).name);
end
lp03Mat = string(lp03Mat);

outDir = fullfile(dataRoot, "LP04_" + string(datetime("now","Format","yyyyMMdd_HHmmss")));
mkdir(outDir);
diary(char(fullfile(outDir,"LP04_log.txt"))); cleanupDiary = onCleanup(@() diary("off"));
fprintf("LP04 started %s\nMATLAB %s\nInput : %s\nOutput: %s\n\n", string(datetime("now")), version, lp03Mat, outDir);

S = load(lp03Mat, "data", "predictors");
data = S.data; pred = S.predictors;

% ---- locked settings (V2 9G / 9H1 / 9H2) ----
cfg = struct();
cfg.maxTrainMissingFraction = 0.00;
cfg.nearZeroUniqueFraction  = 0.01;
cfg.nearZeroFrequencyRatio  = 95;
cfg.correlationThreshold    = 0.95;
cfg.seed        = 20260717;
cfg.cycles      = 200;
cfg.learnRate   = 0.05;
cfg.maxSplits   = 31;
cfg.minLeaf     = 5;
cfg.minSensitivity = 0.90;
cfg.pGrid = (0.001:0.001:0.999)';
cfg.sGrid = log(cfg.pGrid ./ (1 - cfg.pGrid));

jobs = struct( ...
    "name",   {"B", "Astd", "Ahist", "Ahist_ablation"}, ...
    "table",  {"internal_B", "internal_Astd", "internal_Ahist", "internal_Ahist"}, ...
    "dropNumeric", {strings(1,0), strings(1,0), strings(1,0), pred.ablation_removed});

models = struct();
summary = table();
allReports = table();

for j = 1:numel(jobs)
    J = jobs(j);
    fprintf("\n============================================================\n");
    fprintf("PIPELINE %s (training table %s)\n", J.name, J.table);
    fprintf("============================================================\n");
    T = data.(J.table);
    Ttr = T(T.split == "train", :);
    Tva = T(T.split == "validation", :);
    if any(ismember(["test","external"], unique([Ttr.split; Tva.split])))
        error("Test/external rows leaked into training or validation.");
    end
    % Leakage guards: no shared rows, studies or patients between train and validation
    if any(ismember(Ttr.row_id, Tva.row_id)), error("Row overlap between train and validation."); end
    if any(ismember(Ttr.case_id, Tva.case_id)), error("Study overlap between train and validation."); end
    if any(ismember(Ttr.patient_id, Tva.patient_id)), error("Patient overlap between train and validation."); end
    fprintf("Leakage guard OK: 0 shared rows, studies or patients between train and validation.\n");
    numVars = setdiff(pred.numeric, J.dropNumeric, "stable");
    catVars = pred.categorical;
    fprintf("Train: %d (%d pos)  Validation: %d (%d pos)  Candidate predictors: %d numeric + %d zones\n", ...
        height(Ttr), sum(Ttr.y==1), height(Tva), sum(Tva.y==1), numel(numVars), numel(catVars));

    % Preprocessing is learned from the TRAINING table only (Tva is never passed in)
    prep = learnPreprocessing(Ttr, numVars, catVars, cfg);
    prepLocked = prep;   % copy taken before validation is touched
    rep =missingnessAndExclusionReport(Ttr, Tva, prep, string(J.name));
    writetable(rep, fullfile(outDir, "LP04_" + J.name + "_missingness_and_exclusions.csv"));
    allReports = [allReports; rep]; %#ok<AGROW>
    Xtr = applyPreprocessing(Ttr, prep);
    Xva = applyPreprocessing(Tva, prep);
    ytr = double(Ttr.y); yva = double(Tva.y);
    % Invariance check: re-learning on training data alone gives identical
    % decisions, and the stored decisions were not altered by transforming validation.
    prep2 = learnPreprocessing(Ttr, numVars, catVars, cfg);
    if ~isequaln(prep2, prepLocked) || ~isequaln(prep, prepLocked)
        error("Preprocessing decisions changed after validation was transformed.");
    end
    if ~isequal(Xtr, applyPreprocessing(Ttr, prep))
        error("Training transform not reproducible.");
    end
    fprintf("Preprocessing check OK: decisions learned on training only; validation transformed with locked training decisions (%d columns).\n", size(Xva,2));

    rng(cfg.seed, "twister");
    tree = templateTree("MaxNumSplits", cfg.maxSplits, "MinLeafSize", cfg.minLeaf, "Surrogate", "off");
    mdl = fitcensemble(Xtr, ytr, "Method","LogitBoost", "Learners",tree, ...
        "NumLearningCycles",cfg.cycles, "LearnRate",cfg.learnRate, "ClassNames",[0 1]);
    if string(mdl.ScoreTransform) ~= "none"
        error("Unexpected ScoreTransform %s", string(mdl.ScoreTransform));
    end
    mdl = compact(mdl);

    sTr = rawMargin(mdl, Xtr);
    sVa = rawMargin(mdl, Xva);
    [~,~,~,aucTr] = perfcurve(ytr, sTr, 1);
    [~,~,~,aucVa] = perfcurve(yva, sVa, 1);

    [thrTable, sel, fallbackUsed] = selectThreshold(sVa, yva, cfg);
    thrRaw = sel.threshold_raw;

    m = struct();
    m.name = J.name;
    m.training_table = J.table;
    m.model = mdl;
    m.preprocessing = prep;
    m.predictors = prep.finalNames;
    m.threshold_raw = thrRaw;
    m.threshold_p_logistic_grid = sel.p_grid;
    m.threshold_doublelogit = 1 / (1 + exp(-2*thrRaw));
    m.selected_row = sel;
    m.fallback_used = fallbackUsed;
    m.train_auc = aucTr;
    m.validation_auc = aucVa;
    m.n_train = height(Ttr); m.n_train_pos = sum(ytr==1);
    m.n_validation = height(Tva); m.n_validation_pos = sum(yva==1);
    m.cfg = cfg;
    models.(J.name) = m;

    writetable(prep.audit, fullfile(outDir, "LP04_" + J.name + "_feature_audit.csv"));
    writetable(prep.corrPairs, fullfile(outDir, "LP04_" + J.name + "_high_correlation_pairs.csv"));
    writetable(table(prep.finalNames(:), 'VariableNames', {'final_predictor'}), ...
        fullfile(outDir, "LP04_" + J.name + "_final_predictors.csv"));
    writetable(thrTable, fullfile(outDir, "LP04_" + J.name + "_validation_threshold_grid.csv"));

    fprintf("Predictors: %d candidate -> removed missing %d, zero-var %d, near-zero %d; zone columns %d; correlation removed %d; FINAL %d\n", ...
        numel(numVars), prep.nRemovedMissing, prep.nRemovedZero, prep.nRemovedNearZero, ...
        numel(prep.encodedNames), sum(prep.corrRemoved), numel(prep.finalNames));
    fprintf("Train AUC %.4f | Validation AUC %.4f\n", aucTr, aucVa);
    fprintf("Threshold: raw s = %.6f (logistic-grid p = %.3f; doublelogit P = %.6f); fallback = %d\n", ...
        thrRaw, sel.p_grid, m.threshold_doublelogit, fallbackUsed);
    fprintf("Validation at threshold: sens %.4f spec %.4f prec %.4f MCC %.4f\n", ...
        sel.sensitivity, sel.specificity, sel.precision, sel.MCC);

    summary = [summary; table(string(J.name), string(J.table), m.n_train, m.n_train_pos, m.n_validation, m.n_validation_pos, ...
        numel(numVars), numel(prep.finalNames), aucTr, aucVa, thrRaw, sel.p_grid, m.threshold_doublelogit, fallbackUsed, ...
        sel.sensitivity, sel.specificity, sel.precision, sel.MCC, ...
        'VariableNames', {'pipeline','training_table','n_train','n_train_pos','n_validation','n_validation_pos', ...
        'n_candidate_numeric_predictors','n_final_predictors','train_auc','validation_auc','threshold_raw', ...
        'threshold_p_logistic_grid','threshold_doublelogit','fallback_used', ...
        'val_sensitivity','val_specificity','val_precision','val_MCC'})]; %#ok<AGROW>
end

% predictor overlap between pipelines (reported, Methods 7)
names = string(fieldnames(models));
ov = table();
for a = 1:numel(names)
    for b = a+1:numel(names)
        pa = models.(names(a)).predictors; pb = models.(names(b)).predictors;
        ov = [ov; table(names(a), names(b), numel(pa), numel(pb), numel(intersect(pa,pb)), ...
            'VariableNames', {'pipeline_1','pipeline_2','n_1','n_2','n_shared'})]; %#ok<AGROW>
    end
end

writetable(allReports, fullfile(outDir, "LP04_ALL_pipelines_missingness_and_exclusions.csv"));
writetable(summary, fullfile(outDir, "LP04_summary.csv"));
writetable(ov, fullfile(outDir, "LP04_predictor_overlap.csv"));
vs = ver("stats");
info = struct("created", string(datetime("now")), "matlab_version", string(version), ...
    "stats_toolbox", string(vs(1).Version), "lp03_mat", lp03Mat);
save(fullfile(outDir, "LP04_models_and_thresholds.mat"), "models", "summary", "ov", "info", "-v7.3");

disp(summary(:, {'pipeline','n_train','n_final_predictors','validation_auc','threshold_raw','threshold_doublelogit','fallback_used','val_sensitivity','val_MCC'}));
disp(ov);
fprintf("\nLP04 COMPLETED (%.1f min). Test and external candidates were not scored.\n", toc(t0)/60);
fprintf("Models: %s\n", fullfile(outDir, "LP04_models_and_thresholds.mat"));
end

%% ======================================================================
function rep = missingnessAndExclusionReport(Ttr, Tva, prep, pipelineName)
% One row per candidate source predictor (numeric and zone) for this
% training pipeline: NaN / Inf counts in training and validation, how many
% validation values were imputed with the locked training median, and the
% final decision with its reason. One row per generated one-hot column.
    srcNames = [prep.numericCandidates(:); prep.categoricalCandidates(:)];
    n = numel(srcNames);
    kind = [repmat("numeric", numel(prep.numericCandidates),1); repmat("zone", numel(prep.categoricalCandidates),1)];
    trNaN = zeros(n,1); trInf = zeros(n,1); vaNaN = zeros(n,1); vaInf = zeros(n,1);
    vaImputed = zeros(n,1); decision = strings(n,1); reason = strings(n,1);
    for i = 1:n
        v = srcNames(i);
        if kind(i) == "numeric"
            xt = double(Ttr.(v)); xv = double(Tva.(v));
            trNaN(i) = sum(isnan(xt)); trInf(i) = sum(isinf(xt));
            vaNaN(i) = sum(isnan(xv)); vaInf(i) = sum(isinf(xv));
            a = prep.audit(prep.audit.feature_name == v, :);
            fd = a.final_decision(1);
            switch fd
                case "retained_final"
                    decision(i) = "retained"; reason(i) = "passed all training-only filters";
                    if ismember(v, prep.finalNames), vaImputed(i) = vaNaN(i) + vaInf(i); end
                case "removed_training_missing_or_nonfinite"
                    decision(i) = "excluded";
                    reason(i) = sprintf("non-finite in training (%d NaN, %d Inf of %d)", trNaN(i), trInf(i), height(Ttr));
                case "removed_zero_variance"
                    decision(i) = "excluded"; reason(i) = "zero variance in training";
                case "removed_near_zero_variance"
                    decision(i) = "excluded";
                    reason(i) = sprintf("near-zero variance in training (unique fraction %.4g <= 0.01; frequency ratio %.4g >= 95)", ...
                        a.train_unique_fraction(1), a.dominant_frequency_ratio(1));
                case "removed_high_training_correlation"
                    decision(i) = "excluded";
                    p = prep.corrPairs(prep.corrPairs.removed_feature == v, :);
                    reason(i) = sprintf("|r| = %.4f >= 0.95 with retained earlier predictor %s", abs(p.pearson_r(1)), p.feature_1(1));
                otherwise
                    decision(i) = "check"; reason(i) = fd;
            end
        else
            st = normalizeCategory(Ttr.(v)); sv = normalizeCategory(Tva.(v));
            trNaN(i) = sum(st == "<missing>"); vaNaN(i) = sum(sv == "<missing>");
            cats = prep.categories.(v);
            unseen = sum(~ismember(sv, [cats(:); "<missing>"]));
            if numel(cats) <= 1
                decision(i) = "excluded"; reason(i) = "zero or one category in training";
            else
                decision(i) = "one-hot encoded";
                reason(i) = sprintf("%d training categories; %d validation values in categories unseen in training (all one-hot columns 0)", numel(cats), unseen);
            end
        end
    end
    rep = table(repmat(pipelineName,n,1), srcNames, kind, repmat(height(Ttr),n,1), trNaN, trInf, ...
        repmat(height(Tva),n,1), vaNaN, vaInf, vaImputed, decision, reason, ...
        'VariableNames', {'pipeline','feature','type','n_train','train_NaN_or_missing','train_Inf', ...
        'n_validation','validation_NaN_or_missing','validation_Inf','validation_values_imputed_with_training_median', ...
        'decision','reason'});

    % one-hot columns
    for k = 1:numel(prep.encodedNames)
        e = prep.encodedNames(k);
        if ismember(e, prep.finalNames)
            d = "retained"; r = "one-hot column from " + prep.encodedSource(k);
        else
            p = prep.corrPairs(prep.corrPairs.removed_feature == e, :);
            d = "excluded"; r = sprintf("one-hot column; |r| = %.4f >= 0.95 with retained earlier predictor %s", abs(p.pearson_r(1)), p.feature_1(1));
        end
        rep = [rep; table(pipelineName, e, "one-hot", height(Ttr), 0, 0, height(Tva), 0, 0, 0, d, string(r), ...
            'VariableNames', rep.Properties.VariableNames)]; %#ok<AGROW>
    end

    % console report
    fprintf("\n--- %s: missing/non-finite values and exclusions ---\n", pipelineName);
    m = rep.type ~= "one-hot";
    fprintf("Source predictors with any NaN/Inf in TRAINING  : %d\n", sum(m & (rep.train_NaN_or_missing + rep.train_Inf) > 0));
    fprintf("Source predictors with any NaN/Inf in VALIDATION: %d (values imputed with training medians: %d)\n", ...
        sum(m & (rep.validation_NaN_or_missing + rep.validation_Inf) > 0), sum(rep.validation_values_imputed_with_training_median));
    ex = rep(rep.decision == "excluded", :);
    fprintf("Excluded: %d\n", height(ex));
    for k = 1:height(ex)
        fprintf("  - %-55s %s\n", ex.feature(k), ex.reason(k));
    end
    im = rep(rep.validation_values_imputed_with_training_median > 0, :);
    for k = 1:height(im)
        fprintf("  * imputed in validation: %-40s %d value(s)\n", im.feature(k), im.validation_values_imputed_with_training_median(k));
    end
    fprintf("Retained (final model predictors): %d\n", sum(rep.decision == "retained"));
end

function prep = learnPreprocessing(Ttr, numVars, catVars, cfg)
    nNum = numel(numVars);
    missFrac = zeros(nNum,1); uFrac = zeros(nNum,1); uCount = zeros(nNum,1);
    vr = zeros(nNum,1); domRatio = zeros(nNum,1);
    isZero = false(nNum,1); isNear = false(nNum,1); rmMiss = false(nNum,1);
    for i = 1:nNum
        x = double(Ttr.(numVars(i)));
        bad = ~isfinite(x);
        missFrac(i) = mean(bad);
        xv = x(~bad);
        if isempty(xv)
            isZero(i) = true; isNear(i) = true; domRatio(i) = Inf; vr(i) = NaN;
        else
            [u,~,g] = unique(xv);
            c = sort(accumarray(g,1), "descend");
            uCount(i) = numel(u); uFrac(i) = numel(u)/numel(xv); vr(i) = var(xv,0);
            if numel(c) == 1, domRatio(i) = Inf; else, domRatio(i) = c(1)/c(2); end
            isZero(i) = numel(u) <= 1 || ~isfinite(vr(i)) || vr(i) <= eps;
            isNear(i) = uFrac(i) <= cfg.nearZeroUniqueFraction && domRatio(i) >= cfg.nearZeroFrequencyRatio;
        end
        rmMiss(i) = missFrac(i) > cfg.maxTrainMissingFraction;
    end
    rmNum = rmMiss | isZero | isNear;
    keepNum = numVars(~rmNum);
    decision = repmat("retained_initially", nNum,1);
    decision(rmMiss) = "removed_training_missing_or_nonfinite";
    decision(~rmMiss & isZero) = "removed_zero_variance";
    decision(~rmMiss & ~isZero & isNear) = "removed_near_zero_variance";

    % one-hot zones, training categories in stable order (9G)
    cats = struct(); encNames = strings(1,0); encSrc = strings(1,0); encCat = strings(1,0);
    catDecision = strings(numel(catVars),1);
    for i = 1:numel(catVars)
        v = catVars(i);
        tv = normalizeCategory(Ttr.(v));
        tc = unique(tv, "stable"); tc(tc == "<missing>") = [];
        cats.(v) = tc;
        if numel(tc) <= 1
            catDecision(i) = "removed_zero_or_single_training_category";
            continue;
        end
        catDecision(i) = "retained_and_one_hot_encoded";
        for c = 1:numel(tc)
            nm = string(matlab.lang.makeValidName(char(v + "__" + string(matlab.lang.makeValidName(char(tc(c)))))));
            encNames(end+1) = nm; encSrc(end+1) = v; encCat(end+1) = tc(c); %#ok<AGROW>
        end
    end

    Xnum = double(table2array(Ttr(:, cellstr(keepNum))));
    medians = median(Xnum, 1, "omitnan");
    Xenc = encode(Ttr, encSrc, encCat);
    allNames = [keepNum, encNames];
    X = [Xnum, Xenc];
    if size(X,2) > 1
        Rm = corr(X, "Rows","pairwise", "Type","Pearson");
    else
        Rm = 1;
    end
    [corrRemoved, corrPairs] = greedyCorrelationFilter(Rm, allNames, cfg.correlationThreshold);
    finalNames = allNames(~corrRemoved);

    audit = table(numVars(:), missFrac, uCount, uFrac, vr, domRatio, isZero, isNear, decision, ...
        'VariableNames', {'feature_name','train_missing_fraction','train_unique_count','train_unique_fraction', ...
        'train_variance','dominant_frequency_ratio','is_zero_variance','is_near_zero_variance','initial_decision'});
    audit.final_decision = audit.initial_decision;
    audit.final_decision(ismember(audit.feature_name, allNames(corrRemoved))) = "removed_high_training_correlation";
    audit.final_decision(ismember(audit.feature_name, finalNames)) = "retained_final";
    catAudit = table(catVars(:), NaN(numel(catVars),1), NaN(numel(catVars),1), NaN(numel(catVars),1), ...
        NaN(numel(catVars),1), NaN(numel(catVars),1), false(numel(catVars),1), false(numel(catVars),1), catDecision, catDecision, ...
        'VariableNames', audit.Properties.VariableNames);
    encAudit = table(encNames(:), zeros(numel(encNames),1), NaN(numel(encNames),1), NaN(numel(encNames),1), ...
        NaN(numel(encNames),1), NaN(numel(encNames),1), false(numel(encNames),1), false(numel(encNames),1), ...
        repmat("one_hot_column", numel(encNames),1), repmat("retained_final", numel(encNames),1), ...
        'VariableNames', audit.Properties.VariableNames);
    encAudit.final_decision(ismember(encAudit.feature_name, allNames(corrRemoved))) = "removed_high_training_correlation";
    audit = [audit; catAudit; encAudit];

    prep = struct();
    prep.numericCandidates = numVars;
    prep.categoricalCandidates = catVars;
    prep.keepNumeric = keepNum;
    prep.trainingMedians = medians;
    prep.categories = cats;
    prep.encodedNames = encNames;
    prep.encodedSource = encSrc;
    prep.encodedCategory = encCat;
    prep.allPreparedNames = allNames;
    prep.corrRemoved = corrRemoved;
    prep.corrPairs = corrPairs;
    prep.finalNames = finalNames;
    prep.audit = audit;
    prep.nRemovedMissing = sum(rmMiss);
    prep.nRemovedZero = sum(isZero & ~rmMiss);
    prep.nRemovedNearZero = sum(isNear & ~isZero & ~rmMiss);
    prep.cfg = cfg;
end

function X = applyPreprocessing(T, prep)
    Xnum = double(table2array(T(:, cellstr(prep.keepNumeric))));
    for j = 1:size(Xnum,2)
        bad = ~isfinite(Xnum(:,j));
        if any(bad), Xnum(bad,j) = prep.trainingMedians(j); end
    end
    Xenc = encode(T, prep.encodedSource, prep.encodedCategory);
    Xall = [Xnum, Xenc];
    [~, loc] = ismember(prep.finalNames, prep.allPreparedNames);
    X = Xall(:, loc);
    if any(~isfinite(X), "all"), error("Non-finite values remain after preprocessing."); end
end

function Xenc = encode(T, encSrc, encCat)
    Xenc = zeros(height(T), numel(encSrc));
    for k = 1:numel(encSrc)
        Xenc(:,k) = double(normalizeCategory(T.(encSrc(k))) == encCat(k));
    end
end

function values = normalizeCategory(x)
    values = strtrim(string(x));
    miss = ismissing(values) | strlength(values) == 0 | lower(values) == "nan" | ...
        lower(values) == "na" | lower(values) == "none";
    values(miss) = "<missing>";
end

function [removeMask, pairTable] = greedyCorrelationFilter(R, names, threshold)
% Identical to V2 9G: scan columns left to right; a column is removed when
% it reaches |r| >= threshold with any earlier column that is still kept.
    p = numel(names);
    removeMask = false(1,p);
    f1 = strings(0,1); f2 = strings(0,1); rr = zeros(0,1);
    for j = 2:p
        if removeMask(j), continue; end
        for i = 1:(j-1)
            if removeMask(i), continue; end
            r = R(i,j);
            if isfinite(r) && abs(r) >= threshold
                removeMask(j) = true;
                f1(end+1,1) = names(i); f2(end+1,1) = names(j); rr(end+1,1) = r; %#ok<AGROW>
                break;
            end
        end
    end
    pairTable = table(f1, f2, rr, f2, 'VariableNames', {'feature_1','feature_2','pearson_r','removed_feature'});
end

function s = rawMargin(mdl, X)
    [~, sc] = predict(mdl, X);
    s = double(sc(:, mdl.ClassNames == 1));
end

function [tbl, sel, fallbackUsed] = selectThreshold(s, y, cfg)
    n = numel(cfg.sGrid);
    TP = zeros(n,1); FP = TP; TN = TP; FN = TP;
    sens = NaN(n,1); spec = sens; prec = sens; npv = sens; acc = sens; bal = sens; f1 = sens; mcc = sens;
    for i = 1:n
        yhat = s >= cfg.sGrid(i);
        TP(i) = sum(yhat & y==1); FP(i) = sum(yhat & y==0);
        TN(i) = sum(~yhat & y==0); FN(i) = sum(~yhat & y==1);
        sens(i) = sdiv(TP(i), TP(i)+FN(i));
        spec(i) = sdiv(TN(i), TN(i)+FP(i));
        prec(i) = sdiv(TP(i), TP(i)+FP(i));
        npv(i)  = sdiv(TN(i), TN(i)+FN(i));
        acc(i)  = sdiv(TP(i)+TN(i), TP(i)+FP(i)+TN(i)+FN(i));
        bal(i)  = mean([sens(i), spec(i)], "omitnan");
        f1(i)   = sdiv(2*TP(i), 2*TP(i)+FP(i)+FN(i));
        den = sqrt((TP(i)+FP(i))*(TP(i)+FN(i))*(TN(i)+FP(i))*(TN(i)+FN(i)));
        mcc(i)  = sdiv(TP(i)*TN(i) - FP(i)*FN(i), den);
    end
    tbl = table(cfg.pGrid, cfg.sGrid, TP, FP, TN, FN, sens, spec, prec, npv, acc, bal, f1, mcc, ...
        'VariableNames', {'p_grid','threshold_raw','TP','FP','TN','FN','sensitivity','specificity', ...
        'precision','NPV','accuracy','balanced_accuracy','F1','MCC'});
    validMcc = isfinite(tbl.MCC);
    eligible = tbl.sensitivity >= cfg.minSensitivity & validMcc;
    if any(eligible)
        pool = tbl(eligible,:); fallbackUsed = false;
    else
        pool = tbl(validMcc,:); fallbackUsed = true;
    end
    if isempty(pool), error("No threshold produced a finite MCC."); end
    pool = sortrows(pool, {'MCC','precision','specificity','balanced_accuracy','threshold_raw'}, ...
        {'descend','descend','descend','descend','descend'});
    sel = table2struct(pool(1,:));
    tbl.MCC_is_finite = validMcc;
    tbl.meets_sensitivity_constraint = tbl.sensitivity >= cfg.minSensitivity;
    tbl.selected = tbl.threshold_raw == sel.threshold_raw;
end

function r = sdiv(a, b)
    if b == 0 || ~isfinite(b), r = NaN; else, r = a/b; end
end
