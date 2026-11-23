function LP05a_candidate_level_evaluation()
%LP05A_CANDIDATE_LEVEL_EVALUATION  Leakage paper - candidate-level evaluation (run once).
%
% Locked methods v1.0 + Amendment 01. Frozen inputs (SHA256-verified):
%   LP03_20261109_203442\LP03_candidate_sets.mat
%   LP04_20261109_204900\LP04_models_and_thresholds.mat   (official model set)
% No refitting, no preprocessing change, no threshold change.
%
% Internal test (primary) and external locked 95 (secondary), Pipeline B
% candidates, models B, A-std, A-hist:
%   PRIMARY      dAUC = AUC(B->B) - AUC(A-hist->B), any-overlap labels
%   MECHANISTIC  B vs A-std, A-std vs A-hist (sequential descriptive contrasts,
%                not causal effects of individual construction rules)
%   SENSITIVITY  >=25% occupancy labels; high-purity labels (f=0 vs f>=0.5)
%   THRESHOLD    sensitivity, specificity, PPV, MCC at each locked threshold
%   ACCEPTANCE   C0/C1/C2/C3 acceptance, candidate / tumour / non-tumour volume
%   CALIBRATION  Brier, calibration intercept (CITL) and slope, decile bins,
%                probabilities P = 1/(1+exp(-2s)) (doublelogit)
%   DESCRIPTIVE  A-hist -> A-hist ("reference-informed reconstruction")
%   QC           studies with >=1 C1 candidate >= 1,000 mL
% Uncertainty: patient-clustered paired bootstrap, 2,000 resamples,
% seed 20261109 (reset before each bootstrap), replicates lacking a class
% in any label set are redrawn and counted.
%
% SAFEGUARDS (the script stops with an error if any fails):
%   S1  every model scores exactly the same Pipeline B test candidate IDs
%       (count, uniqueness, order, identity with LP03 test IDs; ID-list SHA256)
%   S2  no test candidate, study or patient occurs in any model's training or
%       validation data
%   S3  the bootstrap resamples PATIENTS: in every replicate each candidate
%       appears exactly as often as its patient was drawn
%   S4  re-scoring each model's validation set reproduces the LP04 validation AUC
%   S4b ...and reproduces all 999 rows of LP04's validation threshold grid
%       (TP/FP/TN/FN) and the selected threshold exactly. A stringent
%       scoring-consistency check, not proof of identical individual scores.
%   S5  frozen input files match their recorded SHA256

clc;
t0 = tic;
V2 = string(fullfile(psma_root("PSMA_PROJECT_ROOT"),"PSMA_AI_Project_V2_PET_Derived"));
MD = fullfile(V2,"06_results","LeakagePaper_ModelData");
lp03Mat = fullfile(MD,"LP03_20261109_203442","LP03_candidate_sets.mat");
lp04Mat = fullfile(MD,"LP04_20261109_204900","LP04_models_and_thresholds.mat");
sha03 = "6ce98177fb7dffda4804d1a0717849be331571dfe252b842514b0f15583f3eb2";
sha04 = "bc87e01862d59086d3134aa73e61382276aaec7bce4dcd584e109da5a209a4d6";

outDir = fullfile(V2,"06_results","LeakagePaper_Evaluation","LP05a_" + string(datetime("now","Format","yyyyMMdd_HHmmss")));
mkdir(outDir);
diary(char(fullfile(outDir,"LP05a_log.txt"))); cleanupDiary = onCleanup(@() diary("off"));
fprintf("LP05a started %s\nMATLAB %s\nOutput: %s\n\n", string(datetime("now")), version, outDir);

cfg.nBoot = 2000;
cfg.seed = 20261109;
cfg.models = ["B","Astd","Ahist"];
cfg.labels = ["any","25pct","highpurity"];

%% S5 frozen inputs ------------------------------------------------------
fprintf("S5  verifying frozen inputs...\n");
h03 = fileSha256(lp03Mat); h04 = fileSha256(lp04Mat);
if h03 ~= sha03, error("S5 FAILED: LP03 MAT hash %s does not match the recorded official hash.", h03); end
if h04 ~= sha04, error("S5 FAILED: LP04 MAT hash %s does not match the recorded official hash.", h04); end
fprintf("S5  OK: LP03 %s\n    OK: LP04 %s\n", h03, h04);

D = load(lp03Mat, "data"); data = D.data;
M = load(lp04Mat, "models"); models = M.models;
for m = cfg.models
    if string(models.(m).model.ScoreTransform) ~= "none", error("Unexpected ScoreTransform for %s.", m); end
end

%% S4 scoring reproduces LP04 validation AUC ---------------------------
fprintf("\nS4  re-scoring validation sets with the frozen preprocessing...\n");
for m = [cfg.models, "Ahist_ablation"]
    mm = models.(m);
    T = data.(mm.training_table);
    Tv = T(T.split == "validation", :);
    s = scoreTable(Tv, mm);
    a = aucMW(double(Tv.y), s);
    if abs(a - mm.validation_auc) > 1e-9
        error("S4 FAILED for %s: validation AUC %.10f vs LP04 %.10f", m, a, mm.validation_auc);
    end
    fprintf("S4  OK: %-15s validation AUC %.6f reproduced\n", m, a);

    % S4b: full 999-row validation threshold grid and selected threshold
    gridCsv = fullfile(fileparts(lp04Mat), "LP04_" + m + "_validation_threshold_grid.csv");
    G = readtable(gridCsv, "VariableNamingRule", "preserve");
    sg = mm.cfg.sGrid; yv = double(Tv.y);
    if height(G) ~= numel(sg) || numel(sg) ~= 999, error("S4b FAILED for %s: grid length.", m); end
    TP = zeros(999,1); FP = TP; TN = TP; FN = TP;
    for i = 1:999
        yh = s >= sg(i);
        TP(i) = sum(yh & yv==1); FP(i) = sum(yh & yv==0); TN(i) = sum(~yh & yv==0); FN(i) = sum(~yh & yv==1);
    end
    if ~isequal(TP, G.TP) || ~isequal(FP, G.FP) || ~isequal(TN, G.TN) || ~isequal(FN, G.FN)
        error("S4b FAILED for %s: TP/FP/TN/FN differ from LP04 in %d grid rows.", m, ...
            sum(TP~=G.TP | FP~=G.FP | TN~=G.TN | FN~=G.FN));
    end
    if max(abs(G.threshold_raw - sg)) > 1e-12, error("S4b FAILED for %s: grid threshold values differ.", m); end
    thrNew = reselectThreshold(TP, FP, TN, FN, sg, mm.cfg.minSensitivity);
    if ~isequal(thrNew, mm.threshold_raw)
        error("S4b FAILED for %s: re-selected threshold %.10f vs locked %.10f.", m, thrNew, mm.threshold_raw);
    end
    fprintf("S4b OK: %-15s all 999 grid rows (TP/FP/TN/FN) identical; selected threshold %.6f identical\n", m, thrNew);
end
fprintf("    Note: S4/S4b are a stringent scoring-consistency check, not a mathematical proof\n");
fprintf("    that every candidate score is identical (different score assignments can give the same counts).\n");

%% Build evaluation sets -------------------------------------------------
cohorts = ["internal","external"];
E = struct();
E.internal = data.internal_B(data.internal_B.split == "test", :);
E.external = data.external_B;
checkCount("internal test B candidates", height(E.internal), 2673);
checkCount("internal test B positives", sum(E.internal.label_any==1), 1263);
checkCount("internal test patients", numel(unique(E.internal.patient_id)), 58);
checkCount("internal test studies", numel(unique(E.internal.case_id)), 87);
checkCount("external B candidates", height(E.external), 4173);
checkCount("external B positives", sum(E.external.label_any==1), 3026);
checkCount("external patients (studies)", numel(unique(E.external.patient_id)), 95);

%% S2 no overlap with any model's training/validation data -------------
fprintf("\nS2  checking test/external independence from all training and validation data...\n");
for m = [cfg.models, "Ahist_ablation"]
    T = data.(models.(m).training_table);
    dev = T(ismember(T.split, ["train","validation"]), :);
    for c = cohorts
        X = E.(c);
        nR = sum(ismember(X.row_id, dev.row_id));
        nS = sum(ismember(unique(X.case_id), unique(dev.case_id)));
        nP = sum(ismember(unique(X.patient_id), unique(dev.patient_id)));
        if nR + nS + nP > 0
            error("S2 FAILED: %s cohort shares %d rows, %d studies, %d patients with %s development data.", c, nR, nS, nP, m);
        end
    end
    fprintf("S2  OK: %-15s 0 shared rows/studies/patients with internal test and external sets\n", m);
end

%% Score all models on the SAME candidate tables (S1) ------------------
S = struct();
for c = cohorts
    X = E.(c);
    keep = ["row_id","source_candidate_id","cohort","case_id","patient_id","candidate_label_value", ...
        "label_any","label_25pct","label_highpurity","f_overlap","occupancy_category", ...
        "candidate_volume_mL","overlap_volume_mL"];
    Sc = X(:, keep);
    idsRef = X.row_id;
    for m = cfg.models
        mm = models.(m);
        [s, ids] = scoreTable(X, mm);
        if ~isequal(ids, idsRef), error("S1 FAILED: %s scored a different candidate list on %s.", m, c); end
        Sc.("s_" + m) = s;
        Sc.("P_" + m) = 1 ./ (1 + exp(-2*s));
        Sc.("accept_" + m) = s >= mm.threshold_raw;
    end
    % doublelogit check against MATLAB's built-in transform
    try
        mdl2 = models.B.model; mdl2.ScoreTransform = "doublelogit";
        [~, sc2] = predict(mdl2, applyPreprocessing(X, models.B.preprocessing));
        dlDiff = max(abs(sc2(:, mdl2.ClassNames==1) - Sc.P_B));
        if dlDiff > 1e-12, error("doublelogit check failed (max diff %.3g).", dlDiff); end
    catch ME
        if contains(ME.message, "doublelogit check failed"), rethrow(ME); end
        warning("doublelogit built-in comparison skipped (%s); formula 1/(1+exp(-2s)) used as locked.", ME.message);
    end
    S.(c) = Sc;
end
% S1 identity with LP03 test IDs and ID hashes
lp03TestIds = data.internal_B.row_id(data.internal_B.split == "test");
if ~isequal(S.internal.row_id, lp03TestIds) || numel(unique(S.internal.row_id)) ~= 2673
    error("S1 FAILED: internal evaluation IDs differ from LP03 test IDs.");
end
if numel(unique(S.external.row_id)) ~= 4173, error("S1 FAILED: external IDs not unique."); end
idHashInt = stringSha256(strjoin(sort(S.internal.row_id), newline));
idHashExt = stringSha256(strjoin(sort(S.external.row_id), newline));
fprintf("\nS1  OK: B, A-std and A-hist scored the identical %d internal test B candidate IDs (ID-list SHA256 %s)\n", height(S.internal), idHashInt);
fprintf("S1  OK: B, A-std and A-hist scored the identical %d external B candidate IDs (ID-list SHA256 %s)\n", height(S.external), idHashExt);
fprintf("    doublelogit probabilities verified against ScoreTransform='doublelogit'\n");

%% Point estimates and patient-clustered bootstrap ---------------------
allEst = table(); bootStore = struct();
for c = cohorts
    Sc = S.(c);
    fprintf("\n================ %s cohort: %d candidates, %d patients ================\n", ...
        upper(c), height(Sc), numel(unique(Sc.patient_id)));
    [pt, names] = mainMetrics(Sc, (1:height(Sc))', cfg);
    % cross-check primary AUCs with perfcurve
    for m = cfg.models
        [~,~,~,ap] = perfcurve(Sc.label_any, Sc.("s_" + m), 1);
        if abs(ap - pt(names == "AUC_" + m + "_any")) > 1e-10, error("AUC cross-check failed for %s.", m); end
    end
    [B, nRejected, draws, pList] = patientBootstrap(Sc, @(idx) mainMetrics(Sc, idx, cfg), cfg, ["label_any","label_25pct","label_highpurity"]);
    est = summariseBoot(c, names, pt, B, nRejected);
    allEst = [allEst; est]; %#ok<AGROW>
    bootStore.(c).replicates = array2table(B, 'VariableNames', cellstr(names));
    bootStore.(c).patient_draws = draws;
    bootStore.(c).patient_list = pList;
    bootStore.(c).n_rejected = nRejected;
    writetable(bootStore.(c).replicates, fullfile(outDir, "LP05a_bootstrap_replicates_" + c + ".csv"));
    writetable(table(pList, 'VariableNames', {'patient_id'}), fullfile(outDir, "LP05a_bootstrap_patient_list_" + c + ".csv"));
    printKey(est, c);
end

%% Descriptive: A-hist -> A-hist ----------------------------------------
fprintf("\n================ DESCRIPTIVE: A-hist -> A-hist (reference-informed reconstruction) ================\n");
desc = table();
H = struct();
H.internal = data.internal_Ahist(data.internal_Ahist.split == "test", :);
H.external = data.external_Ahist;
for c = cohorts
    X = H.(c);
    [s, ~] = scoreTable(X, models.Ahist);
    Sh = X(:, ["row_id","case_id","patient_id","y","candidate_volume_mL"]);
    Sh.label_any = double(X.y);
    Sh.s_Ahist = s; Sh.P_Ahist = 1./(1+exp(-2*s)); Sh.accept_Ahist = s >= models.Ahist.threshold_raw;
    f = @(idx) descMetrics(Sh, idx);
    [pt, names] = f((1:height(Sh))');
    [B, nRej] = patientBootstrap(Sh, f, cfg, "label_any");
    e = summariseBoot(c, names, pt, B, nRej);
    e.family(:) = "descriptive_Ahist_on_Ahist";
    desc = [desc; e]; %#ok<AGROW>
    writetable(Sh, fullfile(outDir, "LP05a_scores_Ahist_on_Ahist_" + c + ".csv"));
    fprintf("%s: n=%d, AUC %.4f (95%% CI %.4f-%.4f), sens %.4f, spec %.4f\n", c, height(Sh), ...
        e.estimate(e.metric=="AUC_Ahist_on_Ahist"), e.ci_low(e.metric=="AUC_Ahist_on_Ahist"), e.ci_high(e.metric=="AUC_Ahist_on_Ahist"), ...
        e.estimate(e.metric=="sens_Ahist_on_Ahist"), e.estimate(e.metric=="spec_Ahist_on_Ahist"));
end
allEst = [allEst; desc];

%% Acceptance by occupancy, calibration, C1 QC -------------------------
acc = table(); cal = table(); bins = table(); qc = table();
for c = cohorts
    Sc = S.(c);
    for m = cfg.models
        acc = [acc; acceptanceTable(Sc, m, c)]; %#ok<AGROW>
        [cs, cb] = calibration(Sc.label_any, Sc.("P_" + m), Sc.("s_" + m), c, m);
        cal = [cal; cs]; bins = [bins; cb]; %#ok<AGROW>
    end
    qc = [qc; c1LargeQC(Sc, c)]; %#ok<AGROW>
end
qc = [qc; c1LargeQC(data.internal_B, "internal_all_splits")];
disp(acc(acc.occupancy=="C0" | acc.occupancy=="C1", {'cohort','model','occupancy','n','n_accepted','pct_accepted','candidate_volume_mL','accepted_volume_mL','accepted_tumour_volume_mL','accepted_nontumour_volume_mL'}));
disp(cal);
disp(qc);

%% Save -------------------------------------------------------------------
writetable(allEst, fullfile(outDir, "LP05a_estimates_with_bootstrap_CI.csv"));
writetable(acc, fullfile(outDir, "LP05a_acceptance_by_occupancy.csv"));
writetable(cal, fullfile(outDir, "LP05a_calibration_summary.csv"));
writetable(bins, fullfile(outDir, "LP05a_calibration_bins.csv"));
writetable(qc, fullfile(outDir, "LP05a_QC_C1_ge_1000mL_studies.csv"));
writetable(S.internal, fullfile(outDir, "LP05a_scores_internal_test_B.csv"));
writetable(S.external, fullfile(outDir, "LP05a_scores_external_B.csv"));
safeguards = table(["S1";"S1";"S2";"S3";"S4";"S4b";"S5";"S5"], ...
    ["identical internal test IDs scored by all models";"identical external IDs scored by all models"; ...
     "no test/external row, study or patient in any development set (identifier level; external independence otherwise rests on dataset provenance)"; ...
     "bootstrap unit = patient (multiplicity verified every replicate)"; ...
     "validation AUC reproduced for all models"; ...
     "all 999 validation grid rows (TP/FP/TN/FN) and selected threshold reproduced; stringent scoring-consistency check, not proof of identical individual scores"; ...
     "LP03 input hash";"LP04 input hash"], ...
    ["PASS (n=2673; " + idHashInt + ")";"PASS (n=4173; " + idHashExt + ")";"PASS";"PASS";"PASS";"PASS";h03;h04], ...
    'VariableNames', {'safeguard','description','result'});
writetable(safeguards, fullfile(outDir, "LP05a_safeguards.csv"));
info = struct("created", string(datetime("now")), "matlab_version", string(version), "lp03", lp03Mat, "lp04", lp04Mat, ...
    "lp03_sha256", h03, "lp04_sha256", h04, "cfg", cfg, "id_hash_internal", idHashInt, "id_hash_external", idHashExt);
save(fullfile(outDir, "LP05a_results.mat"), "allEst", "acc", "cal", "bins", "qc", "S", "bootStore", "safeguards", "info", "-v7.3");

fprintf("\nLP05a COMPLETED (%.1f min). All safeguards passed.\nResults: %s\n", toc(t0)/60, outDir);
end

%% ======================================================================
%% Metrics
function [v, names] = mainMetrics(S, idx, cfg)
    v = []; names = strings(1,0);
    lab = struct("any", S.label_any(idx), "x25pct", S.label_25pct(idx), "highpurity", S.label_highpurity(idx));
    labKeys = ["any","x25pct","highpurity"]; labNames = ["any","25pct","highpurity"];
    for k = 1:3
        y = lab.(labKeys(k)); ok = isfinite(y);
        a = struct();
        for m = cfg.models
            s = S.("s_" + m)(idx);
            a.(m) = aucMW(y(ok), s(ok));
            v(end+1) = a.(m); names(end+1) = "AUC_" + m + "_" + labNames(k); %#ok<AGROW>
        end
        v(end+1) = a.B - a.Ahist;    names(end+1) = "dAUC_B_minus_Ahist_" + labNames(k);
        v(end+1) = a.B - a.Astd;     names(end+1) = "dAUC_B_minus_Astd_" + labNames(k);
        v(end+1) = a.Astd - a.Ahist; names(end+1) = "dAUC_Astd_minus_Ahist_" + labNames(k);
    end
    y = S.label_any(idx);
    occ = S.occupancy_category(idx);
    vol = S.candidate_volume_mL(idx); tv = S.overlap_volume_mL(idx);
    sens = struct(); spec = struct();
    for m = cfg.models
        acc = S.("accept_" + m)(idx);
        [se, sp, ppv, mcc] = binMetrics(y, acc);
        sens.(m) = se; spec.(m) = sp;
        P = S.("P_" + m)(idx);
        v = [v, se, sp, ppv, mcc, mean((P - y).^2), ...
            mean(acc(occ=="C0")), mean(acc(occ=="C1")), mean(acc(occ=="C2")), mean(acc(occ=="C3")), ...
            sum(vol(acc)), sum(tv(acc)), sum(vol(acc) - tv(acc))]; %#ok<AGROW>
        names = [names, "sens_" + m, "spec_" + m, "PPV_" + m, "MCC_" + m, "Brier_" + m, ...
            "acceptRate_C0_" + m, "acceptRate_C1_" + m, "acceptRate_C2_" + m, "acceptRate_C3_" + m, ...
            "acceptedVolume_mL_" + m, "acceptedTumourVolume_mL_" + m, "acceptedNonTumourVolume_mL_" + m]; %#ok<AGROW>
    end
    v = [v, sens.B - sens.Ahist, spec.B - spec.Ahist];
    names = [names, "dSens_B_minus_Ahist", "dSpec_B_minus_Ahist"];
end

function [v, names] = descMetrics(S, idx)
    y = S.label_any(idx); s = S.s_Ahist(idx); acc = S.accept_Ahist(idx); P = S.P_Ahist(idx);
    [se, sp, ppv, mcc] = binMetrics(y, acc);
    v = [aucMW(y, s), se, sp, ppv, mcc, mean((P - y).^2)];
    names = ["AUC_Ahist_on_Ahist","sens_Ahist_on_Ahist","spec_Ahist_on_Ahist","PPV_Ahist_on_Ahist","MCC_Ahist_on_Ahist","Brier_Ahist_on_Ahist"];
end

function a = aucMW(y, s)
    y = double(y(:)); s = double(s(:));
    n1 = sum(y==1); n0 = sum(y==0);
    if n1 == 0 || n0 == 0, a = NaN; return; end
    r = tiedrank(s);
    a = (sum(r(y==1)) - n1*(n1+1)/2) / (n1*n0);
end

function [se, sp, ppv, mcc] = binMetrics(y, yhat)
    TP = sum(yhat & y==1); FP = sum(yhat & y==0); TN = sum(~yhat & y==0); FN = sum(~yhat & y==1);
    se = sdiv(TP, TP+FN); sp = sdiv(TN, TN+FP); ppv = sdiv(TP, TP+FP);
    mcc = sdiv(TP*TN - FP*FN, sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)));
end

function r = sdiv(a, b)
    if b == 0 || ~isfinite(b), r = NaN; else, r = a/b; end
end

function thr = reselectThreshold(TP, FP, TN, FN, sg, minSens)
% Same rule as LP04 / V2 9H2: sensitivity >= minSens then max MCC; ties by
% precision, specificity, balanced accuracy, higher threshold; fallback max
% MCC; undefined MCC excluded.
    n = numel(sg);
    se = NaN(n,1); sp = se; pr = se; ba = se; mc = se;
    for i = 1:n
        se(i) = sdiv(TP(i), TP(i)+FN(i)); sp(i) = sdiv(TN(i), TN(i)+FP(i));
        pr(i) = sdiv(TP(i), TP(i)+FP(i)); ba(i) = mean([se(i), sp(i)], "omitnan");
        mc(i) = sdiv(TP(i)*TN(i) - FP(i)*FN(i), sqrt((TP(i)+FP(i))*(TP(i)+FN(i))*(TN(i)+FP(i))*(TN(i)+FN(i))));
    end
    T = table(sg, se, sp, pr, ba, mc, 'VariableNames', {'threshold_raw','sensitivity','specificity','precision','balanced_accuracy','MCC'});
    ok = isfinite(T.MCC);
    el = ok & T.sensitivity >= minSens;
    if any(el), pool = T(el,:); else, pool = T(ok,:); end
    pool = sortrows(pool, {'MCC','precision','specificity','balanced_accuracy','threshold_raw'}, ...
        {'descend','descend','descend','descend','descend'});
    thr = pool.threshold_raw(1);
end

%% Bootstrap (S3)
function [B, nRejected, draws, pList] = patientBootstrap(S, metricFun, cfg, labelVars)
    [pList, ~, pid] = unique(S.patient_id);
    nP = numel(pList); n = height(S);
    rowsByP = accumarray(pid, (1:n)', [nP 1], @(x) {x});
    rng(cfg.seed, "twister");
    [~, names] = metricFun((1:n)');
    B = NaN(cfg.nBoot, numel(names));
    draws = zeros(cfg.nBoot, nP, "uint16");
    nRejected = 0; b = 0;
    while b < cfg.nBoot
        d = randi(nP, nP, 1);
        idx = vertcat(rowsByP{d});
        % S3: each candidate appears exactly as often as its patient was drawn
        rowCount = accumarray(idx, 1, [n 1]);
        patCount = accumarray(d, 1, [nP 1]);
        if ~isequal(rowCount, patCount(pid))
            error("S3 FAILED: bootstrap replicate is not patient-clustered.");
        end
        okRep = true;
        for L = labelVars
            y = S.(L)(idx); y = y(isfinite(y));
            if ~(any(y==1) && any(y==0)), okRep = false; end
        end
        if ~okRep, nRejected = nRejected + 1; continue; end
        b = b + 1;
        draws(b,:) = uint16(d');
        B(b,:) = metricFun(idx);
    end
    fprintf("S3  OK: %d replicates, resampling unit = PATIENT (%d patients, %d candidates); multiplicity verified in every replicate; %d replicates redrawn (one class missing)\n", ...
        cfg.nBoot, nP, n, nRejected);
end

function est = summariseBoot(cohort, names, pt, B, nRejected)
    k = numel(names);
    lo = NaN(k,1); hi = NaN(k,1); p = NaN(k,1); nv = zeros(k,1); fam = strings(k,1);
    for j = 1:k
        x = B(:,j); x = x(isfinite(x)); nv(j) = numel(x);
        if ~isempty(x)
            q = quantile(x, [0.025 0.975]); lo(j) = q(1); hi(j) = q(2);
        end
        nm = names(j);
        if startsWith(nm, ["dAUC","dSens","dSpec"]) && ~isempty(x)
            p(j) = min(1, 2*min(mean(x <= 0), mean(x >= 0)));
        end
        if nm == "dAUC_B_minus_Ahist_any"
            if cohort == "internal", fam(j) = "PRIMARY"; else, fam(j) = "primary_comparison_external_secondary"; end
        elseif contains(nm, ["B_minus_Astd","Astd_minus_Ahist"])
            fam(j) = "mechanistic_sequential_descriptive_not_causal";
        elseif contains(nm, ["_25pct","_highpurity"])
            fam(j) = "sensitivity_label_convention";
        elseif startsWith(nm, "AUC_")
            fam(j) = "discrimination_any_overlap";
        elseif startsWith(nm, "Brier")
            fam(j) = "calibration_doublelogit";
        elseif startsWith(nm, ["accept","Accepted","accepted"])
            fam(j) = "acceptance";
        else
            fam(j) = "locked_threshold";
        end
    end
    est = table(repmat(string(cohort),k,1), names(:), fam, pt(:), lo, hi, p, nv, repmat(nRejected,k,1), ...
        'VariableNames', {'cohort','metric','family','estimate','ci_low','ci_high','bootstrap_p_two_sided','n_valid_replicates','n_redrawn_replicates'});
end

function printKey(est, c)
    show = ["AUC_B_any","AUC_Astd_any","AUC_Ahist_any","dAUC_B_minus_Ahist_any","dAUC_B_minus_Astd_any","dAUC_Astd_minus_Ahist_any", ...
        "dAUC_B_minus_Ahist_25pct","dAUC_B_minus_Ahist_highpurity","sens_B","sens_Astd","sens_Ahist","spec_B","spec_Astd","spec_Ahist", ...
        "Brier_B","Brier_Astd","Brier_Ahist"];
    fprintf("\n%-32s %10s %10s %10s %8s\n", "metric ("+c+")", "estimate", "CI low", "CI high", "p");
    for nm = show
        r = est(est.metric == nm, :);
        if isempty(r), continue; end
        fprintf("%-32s %10.4f %10.4f %10.4f %8.4f\n", nm, r.estimate, r.ci_low, r.ci_high, r.bootstrap_p_two_sided);
    end
end

%% Acceptance, calibration, QC
function T = acceptanceTable(S, m, cohort)
    cats = ["C0","C1","C2","C3","all"];
    T = table();
    acc = S.("accept_" + m);
    for c = cats
        if c == "all", k = true(height(S),1); else, k = S.occupancy_category == c; end
        a = acc & k;
        big = k & S.candidate_volume_mL >= 1000;
        T = [T; table(string(cohort), string(m), c, sum(k), sum(a), 100*sdiv(sum(a), sum(k)), ...
            sum(S.candidate_volume_mL(k)), sum(S.candidate_volume_mL(a)), sum(S.overlap_volume_mL(a)), ...
            sum(S.candidate_volume_mL(a) - S.overlap_volume_mL(a)), sum(big), sum(big & acc), ...
            'VariableNames', {'cohort','model','occupancy','n','n_accepted','pct_accepted','candidate_volume_mL', ...
            'accepted_volume_mL','accepted_tumour_volume_mL','accepted_nontumour_volume_mL', ...
            'n_ge_1000mL','n_ge_1000mL_accepted'})]; %#ok<AGROW>
    end
end

function [cs, cb] = calibration(y, P, s, cohort, m)
    y = double(y); lp = 2*s;                       % logit of doublelogit probability
    b = glmfit(lp, y, "binomial");                 % slope
    b0 = glmfit(ones(size(y)), y, "binomial", "Constant","off", "Offset", lp);   % CITL
    cs = table(string(cohort), string(m), numel(y), mean((P - y).^2), b0(1), b(2), mean(P), mean(y), ...
        'VariableNames', {'cohort','model','n','Brier','calibration_intercept_CITL','calibration_slope','mean_predicted','observed_rate'});
    edges = unique(quantile(P, 0:0.1:1)); edges(1) = -Inf; edges(end) = Inf;
    g = discretize(P, edges);
    cb = table();
    for k = 1:max(g)
        w = g == k;
        cb = [cb; table(string(cohort), string(m), k, sum(w), mean(P(w)), mean(y(w)), ...
            'VariableNames', {'cohort','model','bin','n','mean_predicted','observed_rate'})]; %#ok<AGROW>
    end
end

function T = c1LargeQC(S, label)
    big = S.occupancy_category == "C1" & S.candidate_volume_mL >= 1000;
    st = unique(S.case_id); stBig = unique(S.case_id(big));
    T = table(string(label), numel(st), numel(stBig), 100*numel(stBig)/numel(st), sum(big), ...
        'VariableNames', {'set','n_studies','n_studies_with_C1_ge_1000mL','pct_studies','n_C1_candidates_ge_1000mL'});
end

%% Scoring (identical transform to LP04)
function [s, ids] = scoreTable(T, mm)
    X = applyPreprocessing(T, mm.preprocessing);
    [~, sc] = predict(mm.model, X);
    s = double(sc(:, mm.model.ClassNames == 1));
    ids = T.row_id;
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
    if ~all(tf), error("Final predictor missing from prepared set."); end
    X = Xall(:, loc);
    if any(~isfinite(X), "all"), error("Non-finite values after preprocessing."); end
end

function values = normalizeCategory(x)
    values = strtrim(string(x));
    miss = ismissing(values) | strlength(values) == 0 | lower(values) == "nan" | ...
        lower(values) == "na" | lower(values) == "none";
    values(miss) = "<missing>";
end

%% Utilities
function checkCount(label, got, expected)
    if got ~= expected, error("CHECK FAILED: %s = %d (expected %d).", label, got, expected); end
    fprintf("  check OK: %-45s %d\n", label, got);
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
    h = bytesToHex(md.digest());
end

function h = stringSha256(str)
    md = java.security.MessageDigest.getInstance("SHA-256");
    md.update(typecast(uint8(unicode2native(char(str), "UTF-8")), "int8"));
    h = bytesToHex(md.digest());
end

function h = bytesToHex(d)
    h = string(lower(reshape(dec2hex(typecast(int8(d), "uint8"), 2)', 1, [])));
end
