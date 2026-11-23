function LP05aD2_SUVmin_band_diagnostic()
%LP05AD2_SUVMIN_BAND_DIAGNOSTIC  Leakage paper - FINAL read-only mechanism diagnostic.
%
% Post hoc, descriptive, hypothesis-generating. Run once.
% Pre-run commitment: LeakagePaper_Evaluation\LP05aD2_PRE_RUN_COMMITMENT.md
%   (SHA256 3142c748f907a28988c5587cac94f5ef47184416c6c7ccdc586d731064ec9c0c), verified below.
%
% Fixed SUVmin bands (PET_original_firstorder_Minimum):
%   1  < 2.95      2  2.95 to < 3.00      3  3.00 to <= 3.05      4  > 3.05
% Band edges were chosen AFTER reviewing summary SUVmin distributions
% (LP05a-D Part 1) and BEFORE examining any band-specific acceptance.
%
% Populations: internal reference-derived positives (A-hist test; A-hist
% training shown separately, in-sample), external reference-derived
% positives, internal and external Pipeline B positives.
% Per population x band: n, %, median SUVmin, median volume, median SUVmean,
% and for B, A-std, A-hist, A-hist ablation: % accepted at locked threshold
% and median raw score.
%
% Primary descriptive contrast: band 3 vs band 4 WITHIN the same population.
% No inference, no bootstrap, no model or threshold change.
% STOPPING RULE: mechanistic exploration ends here, whatever the result.

clc;
t0 = tic;
V2 = string(fullfile(psma_root("PSMA_PROJECT_ROOT"),"PSMA_AI_Project_V2_PET_Derived"));
MD = fullfile(V2,"06_results","LeakagePaper_ModelData");
EV = fullfile(V2,"06_results","LeakagePaper_Evaluation");
lp03Mat = fullfile(MD,"LP03_20261109_203442","LP03_candidate_sets.mat");
lp04Mat = fullfile(MD,"LP04_20261109_204900","LP04_models_and_thresholds.mat");
commitFile = fullfile(EV,"LP05aD2_PRE_RUN_COMMITMENT.md");
shaCommit = "3142c748f907a28988c5587cac94f5ef47184416c6c7ccdc586d731064ec9c0c";
sha03 = "6ce98177fb7dffda4804d1a0717849be331571dfe252b842514b0f15583f3eb2";
sha04 = "bc87e01862d59086d3134aa73e61382276aaec7bce4dcd584e109da5a209a4d6";

if ~isempty(dir(fullfile(EV,"LP05aD2_2*")))
    error("An LP05a-D2 output folder already exists. This diagnostic is run once only.");
end
outDir = fullfile(EV, "LP05aD2_" + string(datetime("now","Format","yyyyMMdd_HHmmss")));
mkdir(outDir);
diary(char(fullfile(outDir,"LP05aD2_log.txt"))); cleanupDiary = onCleanup(@() diary("off"));
fprintf("LP05a-D2 (final read-only diagnostic) started %s\nMATLAB %s\nOutput: %s\n\n", string(datetime("now")), version, outDir);

if fileSha256(commitFile) ~= shaCommit, error("Pre-run commitment file missing or changed."); end
if fileSha256(lp03Mat) ~= sha03, error("LP03 hash mismatch."); end
if fileSha256(lp04Mat) ~= sha04, error("LP04 hash mismatch."); end
fprintf("Verified: pre-run commitment %s\n          LP03 %s\n          LP04 %s\n", shaCommit, sha03, sha04);
fprintf("Band edges fixed after reviewing summary SUVmin distributions, before examining band-specific acceptance.\n\n");

D = load(lp03Mat, "data"); data = D.data;
M = load(lp04Mat, "models"); models = M.models;
mods = ["B","Astd","Ahist","Ahist_ablation"];
FMIN = "PET_original_firstorder_Minimum";

P = struct();
P.internal_ref_pos_TEST     = data.internal_Ahist(data.internal_Ahist.split=="test" & data.internal_Ahist.y==1, :);
P.internal_ref_pos_TRAIN_insample = data.internal_Ahist(data.internal_Ahist.split=="train" & data.internal_Ahist.y==1, :);
P.external_ref_pos          = data.external_Ahist(data.external_Ahist.y==1, :);
P.internal_B_pos_TEST       = data.internal_B(data.internal_B.split=="test" & data.internal_B.y==1, :);
P.external_B_pos            = data.external_B(data.external_B.y==1, :);
pn = string(fieldnames(P))';

bandNames = ["1: <2.95","2: 2.95-<3.00","3: 3.00-<=3.05","4: >3.05"];
R = table();
for p = pn
    T = P.(p);
    x = double(T.(FMIN));
    band = zeros(height(T),1);
    band(x < 2.95) = 1;
    band(x >= 2.95 & x < 3.00) = 2;
    band(x >= 3.00 & x <= 3.05) = 3;
    band(x > 3.05) = 4;
    if any(band == 0), error("Unassigned SUVmin values in %s.", p); end
    S = struct();
    for m = mods, S.(m) = scoreTable(T, models.(m)); end
    for b = 1:4
        k = band == b;
        row = table(p, bandNames(b), sum(k), 100*mean(k), med(x(k)), med(T.candidate_volume_mL(k)), med(T.SUVmean(k)), ...
            'VariableNames', {'population','SUVmin_band','n','pct_of_population','SUVmin_median','volume_mL_median','SUVmean_median'});
        for m = mods
            s = S.(m)(k);
            row.("pct_accepted_" + m) = 100*mean(s >= models.(m).threshold_raw);
            row.("score_median_" + m) = med(s);
        end
        R = [R; row]; %#ok<AGROW>
    end
end

fprintf("SUVmin bands: counts, lesion descriptors and acceptance at locked thresholds\n");
fprintf("(acceptance is NaN where a band is empty; Pipeline B positives can only occupy bands 3-4)\n\n");
disp(R(:, {'population','SUVmin_band','n','pct_of_population','SUVmin_median','volume_mL_median','SUVmean_median'}));
disp(R(:, {'population','SUVmin_band','n','pct_accepted_B','pct_accepted_Astd','pct_accepted_Ahist','pct_accepted_Ahist_ablation'}));
disp(R(:, {'population','SUVmin_band','n','score_median_B','score_median_Astd','score_median_Ahist','score_median_Ahist_ablation'}));

fprintf("\nPrimary descriptive contrast (band 3 vs band 4 within the same population), A-hist acceptance:\n");
for p = pn
    a3 = R.pct_accepted_Ahist(R.population==p & R.SUVmin_band==bandNames(3));
    a4 = R.pct_accepted_Ahist(R.population==p & R.SUVmin_band==bandNames(4));
    n3 = R.n(R.population==p & R.SUVmin_band==bandNames(3));
    n4 = R.n(R.population==p & R.SUVmin_band==bandNames(4));
    fprintf("  %-34s band 3: %6.2f%% (n=%d)   band 4: %6.2f%% (n=%d)\n", p, a3, n3, a4, n4);
end
fprintf("\nInterpretation rules (fixed in advance): association only, not proof of a numerical cut-off;\n");
fprintf("volume and SUVmean differ between bands and may contribute; cross-population differences may reflect other shifts.\n");

writetable(R, fullfile(outDir,"LP05aD2_SUVmin_band_acceptance.csv"));
save(fullfile(outDir,"LP05aD2_results.mat"), "R", "-v7.3");
fprintf("\nLP05a-D2 COMPLETED (%.1f min). Read-only. STOPPING RULE: no further mechanistic exploration; next step LP05b.\nResults: %s\n", toc(t0)/60, outDir);
end

%% ======================================================================
function v = med(x)
    if isempty(x), v = NaN; else, v = median(x); end
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
    if ~isfile(f), h = "missing"; return; end
    md = java.security.MessageDigest.getInstance("SHA-256");
    fid = fopen(f, "r");
    c = onCleanup(@() fclose(fid));
    while true
        buf = fread(fid, 2^24, "*uint8");
        if isempty(buf), break; end
        md.update(typecast(buf, "int8"));
    end
    d = md.digest();
    h = string(lower(reshape(dec2hex(typecast(int8(d), "uint8"), 2)', 1, [])));
end
