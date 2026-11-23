function LP06_main_display_items()
%LP06_MAIN_DISPLAY_ITEMS  Leakage paper - seven main display items from FROZEN outputs only.
%
%   Fig 1  study design schematic (no data)
%   Table 1 cohorts: studies/patients, CANDIDATE populations, REFERENCE-LESION
%           populations (different units and denominators, kept separate)
%   Fig 2  primary result forest plot (+ prespecified interpretation sentence)
%   Table 2 candidate level: discrimination and calibration in separate blocks
%   Fig 3  EXPLORATORY (post hoc) mechanism panels
%   Table 3 burden (volume errors first), segmentation, detection; internal and
%           external separate; patient-clustered CIs
%   Fig 4  clinical centrepiece: accepted tumour / non-tumour volume by C0-C3,
%           and lesion detection outcomes
%
% No new statistical analysis: tables re-format frozen estimates; figures plot
% frozen estimates or frozen per-candidate / per-study values. Every input
% file is SHA256-verified against the freeze records.

clc;
V2 = string(fullfile(psma_root("PSMA_PROJECT_ROOT"),"PSMA_AI_Project_V2_PET_Derived"));
R = fullfile(V2,"06_results");
EV = fullfile(R,"LeakagePaper_Evaluation");
f = struct();
f.lp03   = fullfile(R,"LeakagePaper_ModelData","LP03_20261109_203442","LP03_candidate_sets.mat");
f.lp04s  = fullfile(R,"LeakagePaper_ModelData","LP04_20261109_204900","LP04_summary.csv");
f.aEst   = fullfile(EV,"LP05a_20261109_210731","LP05a_estimates_with_bootstrap_CI.csv");
f.aAcc   = fullfile(EV,"LP05a_20261109_210731","LP05a_acceptance_by_occupancy.csv");
f.aCal   = fullfile(EV,"LP05a_20261109_210731","LP05a_calibration_summary.csv");
f.aQC    = fullfile(EV,"LP05a_20261109_210731","LP05a_QC_C1_ge_1000mL_studies.csv");
f.d1     = fullfile(EV,"LP05aD_20261109_212046","LP05aD_part1_feature_distributions.csv");
f.d5     = fullfile(EV,"LP05aD_20261109_212046","LP05aD_part5_predictor_importance.csv");
f.d2     = fullfile(EV,"LP05aD2_20261109_214757","LP05aD2_SUVmin_band_acceptance.csv");
f.bEst   = fullfile(EV,"LP05b_20261109_220537","LP05b_estimates_with_bootstrap_CI.csv");
f.bFail  = fullfile(EV,"LP05b_20261109_220537","LP05b_failure_decomposition.csv");
f.bC1    = fullfile(EV,"LP05b_20261109_220537","LP05b_C1_only_detection.csv");
f.bNeg   = fullfile(EV,"LP05b_20261109_220537","LP05b_tumour_negative_studies.csv");
h = struct( ...
 "lp03","6ce98177fb7dffda4804d1a0717849be331571dfe252b842514b0f15583f3eb2", ...
 "lp04s","2b8a4ee6d9d2dc8310e2e0adf92e1bc446c5992f51ea2df92e6fe4d9f7e32127", ...
 "aEst","2a7c792ca3a96443e523c460eb2ce9569f09e0aa91585755942fc0cd67ab1e56", ...
 "aAcc","2f36134662846e11931ecfb5de550589dda7c9e4c4995d1b324d21acf901a4b8", ...
 "aCal","71b5731c1b12c79a0220b6d5999f2aaf26bed554b5f99b7e19aadeb50eac4337", ...
 "aQC","97224d356870162f0f434aac4851182d15e6cf1521254fb053427faa6d499596", ...
 "d1","91e6aa656ae207a16ca9467ea6772fa5082d6d38b3334ca46b5f27011383d851", ...
 "d5","70c01488b950d82c6f7cbcb85597c47d791a720ec8827533ffc95f4119f3091d", ...
 "d2","be4a1b35e17464a9ae86a58ee55ec1b1e1c1ce968cb2287d68d9e6b1c479073a", ...
 "bEst","88131b7d4f102c70a4d1acc9bf5085b6d6c6d471b7219ebe3df10b58e596fb89", ...
 "bFail","52980c24443a48bd0d9cd299de069ed2c80356b4b10e141bef5275e9a219c6e4", ...
 "bC1","4b02e63e9878d9ea2374e1365ff80e7317a989a53256be56ef396ad0e9259787", ...
 "bNeg","fd1d4dbace207c40c1813952a5de0e99115853778e86d02fdcf8462e8a76b913");
for k = string(fieldnames(f))'
    if fileSha256(f.(k)) ~= h.(k), error("Input hash mismatch: %s", f.(k)); end
end

outDir = fullfile(V2,"LeakagePaper","Manuscript","LP06_main_" + string(datetime("now","Format","yyyyMMdd_HHmmss")));
mkdir(outDir);
diary(char(fullfile(outDir,"LP06_main_log.txt"))); cleanupDiary = onCleanup(@() diary("off"));
fprintf("LP06 main display items %s\nMATLAB %s\nAll %d input files verified against frozen SHA256.\nOutput: %s\n\n", ...
    string(datetime("now")), version, numel(fieldnames(f)), outDir);

rd = @(p) readtable(p, "TextType","string", "VariableNamingRule","preserve");
D = load(f.lp03, "data", "cnt"); data = D.data; cnt = D.cnt;
L4 = rd(f.lp04s); AE = rd(f.aEst); AA = rd(f.aAcc); AC = rd(f.aCal); QC = rd(f.aQC);
D1 = rd(f.d1); D5 = rd(f.d5); D2 = rd(f.d2);
BE = rd(f.bEst); BF = rd(f.bFail); BC = rd(f.bC1); BN = rd(f.bNeg);
mods = ["B","Astd","Ahist"]; modLab = ["B (PET-only)","A-std","A-hist"];
coh = ["internal","external"];
cap = strings(0,1);

%% ===================== TABLE 1 =====================
T1 = table('Size',[0 7], 'VariableTypes',repmat("string",1,7), ...
    'VariableNames',{'Panel','Item','Unit','Internal_train','Internal_validation','Internal_test','External'});
add1 = @(T,p,it,u,a,b,c,d) [T; cellfun(@string, {p,it,u,a,b,c,d}, 'UniformOutput', false)];
gc = @(set,sp,col) cnt.(col)(cnt.set==set & cnt.split==sp);
nb = @(x) string(sprintf('%d', x));
T1 = add1(T1,"A. Studies and patients","Studies","study", nb(gc("internal_B","train","n_studies")), nb(gc("internal_B","validation","n_studies")), nb(gc("internal_B","test","n_studies")), nb(gc("external_B","external","n_studies")));
T1 = add1(T1,"A. Studies and patients","Patients","patient", nb(gc("internal_B","train","n_patients")), nb(gc("internal_B","validation","n_patients")), nb(gc("internal_B","test","n_patients")), nb(gc("external_B","external","n_patients")));
negI = BN.n_tumour_negative_studies(BN.cohort=="internal" & BN.model=="B");
T1 = add1(T1,"A. Studies and patients","Tumour-negative studies","study","—","—", nb(negI), "0");
setNames = ["B","Astd","Ahist"]; setLab = ["Pipeline B (PET-only)","Pipeline A-std","Pipeline A-hist"];
for s = 1:3
    for kind = ["n","n_pos","n_neg"]
        switch kind
            case "n", lab = setLab(s) + ": all candidates";
            case "n_pos", lab = setLab(s) + ": tumour-positive candidates";
            otherwise, lab = setLab(s) + ": tumour-negative candidates";
        end
        v = [gc("internal_"+setNames(s),"train",kind), gc("internal_"+setNames(s),"validation",kind), gc("internal_"+setNames(s),"test",kind)];
        if setNames(s) == "B", e = gc("external_B","external",kind);
        elseif setNames(s) == "Ahist", e = gc("external_Ahist","external",kind);
        else, e = NaN; end
        T1 = add1(T1,"B. Candidate populations (unit: candidate)",lab,"candidate", nb(v(1)), nb(v(2)), nb(v(3)), nbNA(e));
    end
end
for oc = ["C0","C1","C2","C3"]
    ni = AA.n(AA.cohort=="internal" & AA.model=="B" & AA.occupancy==oc); ti = AA.n(AA.cohort=="internal" & AA.model=="B" & AA.occupancy=="all");
    ne = AA.n(AA.cohort=="external" & AA.model=="B" & AA.occupancy==oc); te = AA.n(AA.cohort=="external" & AA.model=="B" & AA.occupancy=="all");
    T1 = add1(T1,"B. Candidate populations (unit: candidate)","Pipeline B occupancy " + oc + " " + occDef(oc),"candidate","—","—", ...
        sprintf('%d (%.1f%%)', ni, 100*ni/ti), sprintf('%d (%.1f%%)', ne, 100*ne/te));
end
qi = QC(QC.set=="internal",:); qe = QC(QC.set=="external",:);
T1 = add1(T1,"B. Candidate populations (unit: candidate)","Studies with >= 1 C1 candidate >= 1,000 mL","study","—","—", ...
    sprintf('%d/%d (%.1f%%)', qi.n_studies_with_C1_ge_1000mL, qi.n_studies, qi.pct_studies), sprintf('%d/%d (%.1f%%)', qe.n_studies_with_C1_ge_1000mL, qe.n_studies, qe.pct_studies));
T1 = add1(T1,"B. Candidate populations (unit: candidate)","Reference-derived positives with SUVmin < 3","candidate", ...
    sprintf('%.1f%%', D1.pct_SUVmin_lt3(D1.population=="Ahist_train_pos")), "—", sprintf('%.1f%%', D1.pct_SUVmin_lt3(D1.population=="Ahist_test_pos")), sprintf('%.1f%%', D1.pct_SUVmin_lt3(D1.population=="ext_A_pos")));
fb = @(c,col) BF.(col)(BF.cohort==c & BF.model=="B");
for it = ["all","ge1mL"]
    lab = "Reference lesions, all"; if it=="ge1mL", lab = "Reference lesions >= 1 mL"; end
    T1 = add1(T1,"C. Reference-lesion populations (unit: lesion; test and external only)",lab,"lesion","—","—", nb(fb("internal","lesions_"+it)), nb(fb("external","lesions_"+it)));
end
T1 = add1(T1,"C. Reference-lesion populations (unit: lesion; test and external only)","Lesions >= 1 mL reachable by any Pipeline B candidate","lesion","—","—", ...
    sprintf('%d (%.1f%%)', fb("internal","lesions_ge1mL")-fb("internal","candidate_generation_miss_ge1mL"), 100*(1-fb("internal","candidate_generation_miss_ge1mL")/fb("internal","lesions_ge1mL"))), ...
    sprintf('%d (%.1f%%)', fb("external","lesions_ge1mL")-fb("external","candidate_generation_miss_ge1mL"), 100*(1-fb("external","candidate_generation_miss_ge1mL")/fb("external","lesions_ge1mL"))));
writeT(T1, outDir, "Table1_cohorts_candidates_lesions");
cap(end+1) = "Table 1. Studies, candidate populations (unit: candidate) and reference-lesion populations (unit: lesion). Candidate and lesion counts have different units and denominators. Reference lesions = 26-connected components of the reference mask; reachable = overlapped by at least one Pipeline B candidate (SUV >= 3, >= 1 mL, 26-connectivity); reachability is specific to this frozen candidate set. A-std is not built for the external cohort (5 mm rule excluded 0 external negatives).";

%% ===================== TABLE 2 =====================
T2 = table('Size',[0 4], 'VariableTypes',repmat("string",1,4), 'VariableNames',{'Block','Item','Internal_test','External'});
ge = @(E,c,m) E(E.cohort==c & E.metric==m, :);
for j = 1:3
    m = mods(j);
    T2 = [T2; {"Discrimination (AUC, any-overlap labels)", modLab(j)+" -> B candidates", fe(ge(AE,"internal","AUC_"+m+"_any"),3), fe(ge(AE,"external","AUC_"+m+"_any"),3)}]; %#ok<AGROW>
end
dl = ["dAUC_B_minus_Ahist_any","PRIMARY: dAUC B - A-hist (any overlap)"; "dAUC_B_minus_Ahist_25pct","dAUC B - A-hist (>= 25% occupancy labels)"; ...
      "dAUC_B_minus_Ahist_highpurity","dAUC B - A-hist (high-purity labels)"; "dAUC_B_minus_Astd_any","Sequential, descriptive: dAUC B - A-std"; "dAUC_Astd_minus_Ahist_any","Sequential, descriptive: dAUC A-std - A-hist"];
for r = 1:size(dl,1)
    T2 = [T2; {"Differences (paired)", dl(r,2), fe(ge(AE,"internal",dl(r,1)),3,true), fe(ge(AE,"external",dl(r,1)),3,true)}]; %#ok<AGROW>
end
thrM = ["sens","Sensitivity"; "spec","Specificity"; "PPV","PPV"; "MCC","MCC"];
for j = 1:3
    for q = 1:size(thrM,1)
        T2 = [T2; {"At locked source threshold", modLab(j)+": "+thrM(q,2), fe(ge(AE,"internal",thrM(q,1)+"_"+mods(j)),3), fe(ge(AE,"external",thrM(q,1)+"_"+mods(j)),3)}]; %#ok<AGROW>
    end
end
for j = 1:3
    ci = AC(AC.cohort=="internal" & AC.model==mods(j),:); ce = AC(AC.cohort=="external" & AC.model==mods(j),:);
    T2 = [T2; {"Calibration (doublelogit probabilities)", modLab(j)+": Brier", fe(ge(AE,"internal","Brier_"+mods(j)),3), fe(ge(AE,"external","Brier_"+mods(j)),3)}]; %#ok<AGROW>
    T2 = [T2; {"Calibration (doublelogit probabilities)", modLab(j)+": intercept (CITL) / slope", string(sprintf('%.2f / %.2f', ci.calibration_intercept_CITL, ci.calibration_slope)), string(sprintf('%.2f / %.2f', ce.calibration_intercept_CITL, ce.calibration_slope))}]; %#ok<AGROW>
end
writeT(T2, outDir, "Table2_candidate_level");
cap(end+1) = "Table 2. Candidate-level performance on the identical Pipeline B candidates (internal test n = 2,673; external n = 4,173). Estimates with 95% patient-clustered bootstrap CIs (2,000 resamples). Discrimination and threshold/calibration are reported separately. Sequential contrasts are descriptive, not causal. p < 0.001 where no bootstrap replicate crossed zero. Calibration intercept/slope are point estimates.";

%% ===================== TABLE 3 =====================
T3 = table('Size',[0 4], 'VariableTypes',repmat("string",1,4), 'VariableNames',{'Block','Item','Internal_test','External'});
rows3 = { ...
 "1. Tumour burden: total tumour volume (all studies)", "Volume bias, mL (predicted - reference)", "bias_volume_", 0;
 "1. Tumour burden: total tumour volume (all studies)", "Lower limit of agreement, mL", "LoAlow_volume_", 0;
 "1. Tumour burden: total tumour volume (all studies)", "Upper limit of agreement, mL", "LoAhigh_volume_", 0;
 "1. Tumour burden: total tumour volume (all studies)", "Mean absolute error, mL", "MAE_volume_", 0;
 "1. Tumour burden: total tumour volume (all studies)", "Lin CCC (headline)", "CCC_volume_", 3;
 "1. Tumour burden: total tumour volume (all studies)", "Summed false-positive volume, mL", "FPvol_total_", 0;
 "2. Segmentation (tumour-positive studies)", "Dice, study-weighted mean (headline)", "Dice_", 3;
 "2. Segmentation (tumour-positive studies)", "Voxel precision", "precision_", 3;
 "2. Segmentation (tumour-positive studies)", "Voxel recall", "recall_", 3;
 "3. Lesion detection", "Recall, lesions >= 1 mL (headline)", "recall_ge1mL_", 3;
 "3. Lesion detection", "Recall, all lesions", "recall_all_", 3;
 "3. Lesion detection", "Recall, reachable lesions >= 1 mL", "recall_reachable_ge1mL_", 3};
for r = 1:size(rows3,1)
    for j = 1:3
        m = mods(j);
        T3 = [T3; {rows3{r,1}, modLab(j)+": "+rows3{r,2}, fe(ge(BE,"internal",rows3{r,3}+m),rows3{r,4}), fe(ge(BE,"external",rows3{r,3}+m),rows3{r,4})}]; %#ok<AGROW>
    end
end
for j = 1:3
    v = BN(BN.cohort=="internal" & BN.model==mods(j),:);
    T3 = [T3; {"4. Tumour-negative studies (internal, n = 9)", modLab(j)+": studies with any FP; FP volume, mL", string(sprintf('%d; %.0f', v.n_with_any_FP, v.total_FP_volume_mL)), "— (none)"}]; %#ok<AGROW>
end
dd = ["d_CCC_volume_B_minus_Ahist","dCCC volume B - A-hist",3; "d_Dice_B_minus_Ahist","dDice B - A-hist",3; "d_recall_ge1mL_B_minus_Ahist","dRecall >= 1 mL B - A-hist",3; ...
      "d_Dice_B_minus_Astd","Sequential, descriptive: dDice B - A-std",3; "d_Dice_Astd_minus_Ahist","Sequential, descriptive: dDice A-std - A-hist",3];
for r = 1:size(dd,1)
    T3 = [T3; {"5. Paired differences", dd(r,2), fe(ge(BE,"internal",dd(r,1)),str2double(dd(r,3)),true), fe(ge(BE,"external",dd(r,1)),str2double(dd(r,3)),true)}]; %#ok<AGROW>
end
writeT(T3, outDir, "Table3_burden_segmentation_detection");
cap(end+1) = "Table 3. Tumour burden, segmentation and lesion detection at the locked thresholds (internal test: 87 studies, 58 patients; external: 95 studies). Volume errors are listed first: the PET-only model overestimates total tumour volume by about 1.4 L per study internally. Estimates with 95% patient-clustered bootstrap CIs. Predicted mask = union of accepted whole candidates. Reachability is specific to the frozen Pipeline B candidate set.";

%% ===================== FIGURE 1 =====================
fig = figure("Color","w","Units","centimeters","Position",[2 2 24 13]);
ax = axes(fig,"Position",[0 0 1 1]); axis(ax,"off"); xlim(ax,[0 1]); ylim(ax,[0 1]);
box1 = @(x,y,w,hh,txt,col) drawBox(ax,x,y,w,hh,txt,col);
box1(0.30,0.86,0.40,0.10,"Same patients and studies (internal: 597 studies, 378 patients; patient-grouped split 423/87/87)",[0.96 0.96 0.96]);
box1(0.02,0.56,0.30,0.22,sprintf("Pipeline B (PET-only)\nCandidates: SUV >= 3, >= 1 mL, 26-conn.\nLabel: any reference overlap\nReference read only for labels"),[0.90 0.90 0.90]);
box1(0.35,0.56,0.30,0.22,sprintf("Pipeline A-std\nPositives: reference components >= 1 mL\nNegatives: B negatives with no voxel\nwithin 5 mm of reference"),[0.80 0.80 0.80]);
box1(0.68,0.56,0.30,0.22,sprintf("Pipeline A-hist (historical rule)\nPositives: reference components >= 1 mL\nNegatives: SUV >= 3 minus dilated\nreference, carved fragments >= 1 mL"),[0.70 0.70 0.70]);
box1(0.15,0.33,0.70,0.14,sprintf("Identical for all pipelines: features (routine + PET/CT radiomics), training-only preprocessing,\nLogitBoost (200 cycles), threshold on own validation set (sensitivity >= 0.90, max MCC)"),[0.96 0.96 0.96]);
box1(0.10,0.06,0.80,0.18,sprintf("Evaluation on the SAME Pipeline B candidates (deployment condition)\nInternal test: 2,673 candidates, 87 studies, 58 patients (primary)   |   External locked 95: 4,173 candidates\nCandidate level (AUC, primary dAUC B - A-hist) -> segmentation, burden, lesion detection"),[0.86 0.86 0.86]);
for x = [0.17 0.50 0.83], annotation(fig,"arrow",[x x],[0.86 0.79]); annotation(fig,"arrow",[x 0.5+(x-0.5)*0.6],[0.56 0.47]); end
annotation(fig,"arrow",[0.5 0.5],[0.33 0.24]);
saveFig(fig, outDir, "Fig1_study_design");
cap(end+1) = "Figure 1. Study design. Three training candidate sets are built from the same patients; everything after candidate construction is identical. All models are evaluated on the same PET-only (Pipeline B) candidates.";

%% ===================== FIGURE 2 =====================
fig = figure("Color","w","Units","centimeters","Position",[2 2 22 11]);
t = tiledlayout(fig,1,2,"TileSpacing","compact","Padding","compact");
nexttile(t); hold on;
lab = ["Internal: any overlap (PRIMARY)","Internal: >= 25% occupancy","Internal: high purity","External: any overlap","External: >= 25% occupancy","External: high purity"];
mt = ["dAUC_B_minus_Ahist_any","dAUC_B_minus_Ahist_25pct","dAUC_B_minus_Ahist_highpurity"];
k = 0; Y = [];
for c = coh
    for m = mt
        k = k + 1; r = ge(AE,c,m); Y(k,:) = [r.estimate r.ci_low r.ci_high]; %#ok<AGROW>
    end
end
yy = numel(lab):-1:1;
errorbar(Y(:,1), yy, Y(:,1)-Y(:,2), Y(:,3)-Y(:,1), "horizontal", "o", "Color","k", "MarkerFaceColor","w", "MarkerEdgeColor","k", "LineWidth",1.2);
plot(Y(1,1), yy(1), "s", "MarkerSize",10, "MarkerFaceColor","k", "MarkerEdgeColor","k");
xline(0,"k--"); yticks(fliplr(yy)); yticklabels(fliplr(lab)); xlabel("\DeltaAUC (B \rightarrow B minus A-hist \rightarrow B)");
title("Primary and label-convention analyses"); grid on; xlim([-0.02 0.16]);
nexttile(t); hold on;
lab2 = ["Internal: B - A-std","Internal: A-std - A-hist","External: B - A-std","External: A-std - A-hist"];
mt2 = ["dAUC_B_minus_Astd_any","dAUC_Astd_minus_Ahist_any"]; Y2 = []; k = 0;
for c = coh, for m = mt2, k = k+1; r = ge(AE,c,m); Y2(k,:) = [r.estimate r.ci_low r.ci_high]; end, end %#ok<AGROW>
yy2 = numel(lab2):-1:1;
errorbar(Y2(:,1), yy2, Y2(:,1)-Y2(:,2), Y2(:,3)-Y2(:,1), "horizontal", "o", "Color","k", "MarkerFaceColor",[0.6 0.6 0.6], "MarkerEdgeColor","k", "LineWidth",1.2);
xline(0,"k--"); yticks(fliplr(yy2)); yticklabels(fliplr(lab2)); xlabel("\DeltaAUC"); grid on; xlim([-0.02 0.16]);
title("Sequential contrasts (descriptive, not causal)");
annotation(fig,"textbox",[0.01 0.0 0.98 0.07],"String", ...
    "\DeltaAUC measures how well each model reproduces Pipeline B's candidate-level label convention under a common deployment condition; it does not by itself establish which model segments tumour better.", ...
    "EdgeColor","none","FontSize",8,"FontAngle","italic","HorizontalAlignment","center");
saveFig(fig, outDir, "Fig2_primary_forest");
cap(end+1) = "Figure 2. Candidate-level AUC differences on identical Pipeline B candidates, with 95% patient-clustered bootstrap CIs. Left: primary endpoint (filled black square) and label-convention sensitivity analyses. Right: sequential contrasts (descriptive, not a causal decomposition). Interpretation (prespecified): dAUC measures how well each model reproduces Pipeline B's candidate-level label convention under a common deployment condition; it does not by itself establish which model segments tumour better.";

%% ===================== FIGURE 3 (EXPLORATORY) =====================
fig = figure("Color","w","Units","centimeters","Position",[2 2 30 11]);
t = tiledlayout(fig,1,3,"TileSpacing","compact","Padding","compact");
title(t, "EXPLORATORY (post hoc): consistent with a construction-associated shortcut; does not establish causality", "FontSize",10,"FontWeight","bold");
nexttile(t);
P = {data.internal_Ahist(data.internal_Ahist.split=="train" & data.internal_Ahist.y==1,:), ...
     data.internal_Ahist(data.internal_Ahist.split=="train" & data.internal_Ahist.y==0,:), ...
     data.internal_B(data.internal_B.split=="test" & data.internal_B.y==1,:), ...
     data.external_Ahist(data.external_Ahist.y==1,:), data.external_B(data.external_B.y==1,:)};
plab = ["A-hist train pos","A-hist train neg","Internal B test pos","External reference pos","External B pos"];
xv = []; gv = [];
for q = 1:numel(P), x = double(P{q}.PET_original_firstorder_Minimum); xv = [xv; x]; gv = [gv; q*ones(numel(x),1)]; end %#ok<AGROW>
boxchart(categorical(gv, 1:numel(P), plab), xv, "MarkerStyle","none", "BoxFaceColor",[0.6 0.6 0.6], "WhiskerLineColor","k", "BoxEdgeColor","k");
yline(3,"k--","LineWidth",1); ylim([1 7]); ylabel("SUVmin of candidate"); title("(a) Minimum SUV by population"); xtickangle(30);
nexttile(t);
pops = ["internal_ref_pos_TEST","external_ref_pos","internal_B_pos_TEST","external_B_pos"];
popL = ["Internal reference pos (test)","External reference pos","Internal B pos (test)","External B pos"];
bands = unique(D2.SUVmin_band, "stable");
A = NaN(numel(pops), numel(bands));
for i = 1:numel(pops), for b = 1:numel(bands)
    r = D2(D2.population==pops(i) & D2.SUVmin_band==bands(b),:);
    if r.n > 0, A(i,b) = r.pct_accepted_Ahist; end
end, end
ax3 = gca; b = bar(ax3, 1:numel(bands), A', "EdgeColor","k");
sh3 = [0.20 0.50 0.80 1.00]; pt3 = ["none","/","x","."];
for q = 1:numel(b), b(q).FaceColor = sh3(q)*[1 1 1]; end
xticks(1:numel(bands)); xticklabels(bands); xlim([0.4 numel(bands)+0.6]); ylim([0 160]); yticks(0:20:100);
ylabel("A-hist acceptance at locked threshold (%)"); title("(b) A-hist acceptance by SUVmin band");
hatchBars(ax3, b, pt3, true);
hatchLegend(ax3, popL, sh3, pt3, [0.03 0.97], 0.55);
nexttile(t);
I = D5(D5.model=="Ahist",:); I = I(1:10,:);
barh(flipud(I.pct_of_total), "FaceColor",[0.55 0.55 0.55], "EdgeColor","k"); yticks(1:10); yticklabels(flipud(replace(I.predictor,"_","\_")));
xlabel("% of total predictor importance"); title("(c) A-hist predictor importance (top 10)");
set(gca,"FontSize",7);
saveFig(fig, outDir, "Fig3_exploratory_mechanism");
cap(end+1) = "Figure 3. EXPLORATORY, post hoc. (a) Candidate minimum SUV; threshold-derived candidates cluster just above SUV 3 (dashed line), manually contoured internal positives spread widely. (b) A-hist acceptance by prespecified SUVmin band (band edges set after viewing summary distributions, before band-specific acceptance); empty bands not shown. Bands contain different candidates, not perturbations of identical candidates. (c) Model-based predictor importance (importance is shared among correlated predictors). Together: strong evidence of a construction-associated shortcut; not proof of what the model learned.";

%% ===================== FIGURE 4 (CLINICAL CENTREPIECE) =====================
fig = figure("Color","w","Units","centimeters","Position",[2 2 30 13]);
t = tiledlayout(fig,2,2,"TileSpacing","compact","Padding","compact");
title(t, "Better lesion detection does not automatically produce reliable tumour-burden quantification", "FontSize",11,"FontWeight","bold");
shV = [0.20 0.97 0.70 0.50 0.85]; ptV = ["none","\","x","-","|"];
for ci = 1:2
    c = coh(ci);
    nexttile(t, ci);
    V = zeros(3,5);
    for j = 1:3
        a = AA(AA.cohort==c & AA.model==mods(j),:);
        V(j,1) = a.accepted_tumour_volume_mL(a.occupancy=="all");
        for o = 1:4, V(j,o+1) = a.accepted_nontumour_volume_mL(a.occupancy=="C"+(o-1)); end
    end
    axA = gca; b = bar(axA, 1:3, V/1000, "stacked", "EdgeColor","k");
    for q = 1:5, b(q).FaceColor = shV(q)*[1 1 1]; end
    xticks(1:3); xticklabels(modLab); xlim([0.4 3.6]); ylim([0 1.75*max(sum(V,2))/1000]);
    ylabel("Accepted volume, L (sum over studies)");
    hatchBars(axA, b, ptV, false);
    bi = ge(BE,c,"bias_volume_B");
    title(sprintf("(%s) %s: accepted volume by occupancy class\nB volume bias %+.0f mL/study (95%% CI %+.0f to %+.0f)", char('a'+ci-1), c, bi.estimate, bi.ci_low, bi.ci_high));
    hatchLegend(axA, ["Tumour (all classes)","Non-tumour, C0 (f = 0)","Non-tumour, C1 (0 < f < 0.25)","Non-tumour, C2","Non-tumour, C3 (f >= 0.5)"], shV, ptV, [0.50 0.97], 0.48);
    nexttile(t, ci+2);
    W = zeros(3,4);
    for j = 1:3
        fr = BF(BF.cohort==c & BF.model==mods(j),:); c1 = BC(BC.cohort==c & BC.model==mods(j),:);
        W(j,:) = [fr.detected_ge1mL - c1.detected_exclusively_via_C1_ge1mL, c1.detected_exclusively_via_C1_ge1mL, ...
                  fr.classifier_miss_ge1mL, fr.candidate_generation_miss_ge1mL] / fr.lesions_ge1mL * 100;
    end
    axB = gca; b = bar(axB, 1:3, W, "stacked", "EdgeColor","k");
    shW = [0.25 0.65 0.92 1.00]; ptW = ["none","x","/","."];
    for q = 1:4, b(q).FaceColor = shW(q)*[1 1 1]; end
    xticks(1:3); xticklabels(modLab); xlim([0.4 3.6]); ylim([0 150]); yticks(0:20:100);
    ylabel("% of reference lesions >= 1 mL");
    hatchBars(axB, b, ptW, false);
    title(sprintf("(%s) %s: outcome of reference lesions >= 1 mL (n = %d)", char('c'+ci-1), c, BF.lesions_ge1mL(BF.cohort==c & BF.model=="B")));
    hatchLegend(axB, ["Detected (not only via C1)","Detected only via accepted C1","Classifier miss","Candidate-generation miss"], shW, ptW, [0.03 0.97], 0.60);
end
saveFig(fig, outDir, "Fig4_clinical_burden_vs_detection");
cap(end+1) = "Figure 4. Clinical consequence. (a, b) Accepted candidate volume at each model's locked threshold, split into reference-tumour volume and non-tumour volume by occupancy class (C0 f = 0; C1 0 < f < 0.25; C2 0.25 <= f < 0.5; C3 f >= 0.5; f = reference-tumour fraction of the candidate); from the frozen candidate-level acceptance table. C1 contributes an important part, not necessarily all, of the non-tumour volume. (c, d) Outcome of reference lesions >= 1 mL: detected (>= 1 voxel overlap with the predicted mask), detected only through accepted C1 candidates, classifier miss, or candidate-generation miss (no Pipeline B candidate overlaps; specific to this frozen candidate set). Bias and CIs from the patient-clustered bootstrap.";

writelines(cap, fullfile(outDir,"LP06_main_captions.txt"));
fprintf("\nLP06 main display items COMPLETED. Frozen outputs only; no new statistical analysis.\nOutput: %s\n", outDir);
end

%% ======================================================================
function s = fe(r, dec, isDiff)
    if nargin < 3, isDiff = false; end
    if isempty(r) || ~isfinite(r.estimate(1)), s = "—"; return; end
    fmt = "%." + dec + "f";
    s = sprintf(fmt + " (" + fmt + " to " + fmt + ")", r.estimate(1), r.ci_low(1), r.ci_high(1));
    if isDiff && isfinite(r.bootstrap_p_two_sided(1))
        p = r.bootstrap_p_two_sided(1);
        if p < 0.001, s = s + "; p < 0.001"; else, s = s + sprintf("; p = %.3f", p); end
    end
    s = string(s);
end

function hatchBars(ax, b, pats, grouped)
% Black hatching inside bar segments (numeric x axis; linear y axis).
    ns = numel(b);
    for q = 1:ns
        if pats(q) == "none", continue; end
        xc = b(q).XEndPoints; yt = b(q).YEndPoints; hgt = b(q).YData;
        if grouped
            w = b(q).BarWidth * min(0.8, ns/(ns+1.5)) / ns;
        else
            w = b(q).BarWidth;
        end
        for k = 1:numel(xc)
            if ~isfinite(hgt(k)) || hgt(k) == 0, continue; end
            hatchRect(ax, xc(k)-w/2, xc(k)+w/2, yt(k)-hgt(k), yt(k), pats(q));
        end
    end
end

function hatchRect(ax, x0, x1, y0, y1, pat)
% Draw a hatch pattern inside a data-space rectangle (normalized spacing).
    xl = xlim(ax); yl = ylim(ax); dx = diff(xl); dy = diff(yl);
    u0 = (x0-xl(1))/dx; u1 = (x1-xl(1))/dx; v0 = (min(y0,y1)-yl(1))/dy; v1 = (max(y0,y1)-yl(1))/dy;
    g = 0.022; X = []; Y = [];
    toData = @(u,v) deal(xl(1)+u*dx, yl(1)+v*dy);
    if any(pat == ["/","x"])
        for c = (v0-u1):g:(v1-u0)
            a = max(u0, v0-c); bb = min(u1, v1-c);
            if a < bb, [xa,ya] = toData([a bb],[a bb]+c); X = [X, xa, NaN]; Y = [Y, ya, NaN]; end %#ok<AGROW>
        end
    end
    if any(pat == ["\","x"])
        for c = (v0+u0):g:(v1+u1)
            a = max(u0, c-v1); bb = min(u1, c-v0);
            if a < bb, [xa,ya] = toData([a bb],c-[a bb]); X = [X, xa, NaN]; Y = [Y, ya, NaN]; end %#ok<AGROW>
        end
    end
    if pat == "|"
        for u = (u0+g/2):g:u1, [xa,ya] = toData([u u],[v0 v1]); X = [X, xa, NaN]; Y = [Y, ya, NaN]; end %#ok<AGROW>
    end
    if pat == "-"
        for v = (v0+g/2):g:v1, [xa,ya] = toData([u0 u1],[v v]); X = [X, xa, NaN]; Y = [Y, ya, NaN]; end %#ok<AGROW>
    end
    if ~isempty(X)
        line(ax, X, Y, "Color","k", "LineWidth",0.4, "HandleVisibility","off");
    end
    if pat == "."
        [uu,vv] = meshgrid((u0+g/2):g:u1, (v0+g/2):g:v1);
        if ~isempty(uu)
            [xa,ya] = toData(uu(:),vv(:));
            line(ax, xa, ya, "LineStyle","none", "Marker",".", "MarkerSize",3, "Color","k", "HandleVisibility","off");
        end
    end
end

function hatchLegend(ax, labels, shades, pats, topLeft, width)
% Legend with grey swatches and the same hatch patterns, drawn in data space.
    hold(ax, "on");
    xl = xlim(ax); yl = ylim(ax); dx = diff(xl); dy = diff(yl);
    sw = 0.06*dx; sh = 0.045*dy; gap = 0.065*dy;
    x0 = xl(1) + topLeft(1)*dx; y0 = yl(1) + topLeft(2)*dy;
    rectangle(ax, "Position",[x0-0.01*dx, y0-numel(labels)*gap-0.01*dy, width*dx, numel(labels)*gap+0.02*dy], ...
        "FaceColor","w", "EdgeColor",[0.4 0.4 0.4]);
    for k = 1:numel(labels)
        yk = y0 - k*gap;
        rectangle(ax, "Position",[x0, yk, sw, sh], "FaceColor",shades(k)*[1 1 1], "EdgeColor","k");
        if pats(k) ~= "none", hatchRect(ax, x0, x0+sw, yk, yk+sh, pats(k)); end
        text(ax, x0+sw+0.015*dx, yk+sh/2, labels(k), "FontSize",7, "VerticalAlignment","middle", "Interpreter","none");
    end
end

function drawBox(ax, x, y, w, hh, txt, col)
    rectangle(ax, "Position",[x y w hh], "Curvature",0.08, "FaceColor",col, "EdgeColor",[0.3 0.3 0.3]);
    text(ax, x+w/2, y+hh/2, txt, "HorizontalAlignment","center", "VerticalAlignment","middle", "FontSize",8.5, "Interpreter","none");
end

function s = nbNA(x)
    if ~isfinite(x), s = "—"; else, s = string(sprintf('%d', x)); end
end

function s = occDef(oc)
    switch oc
        case "C0", s = "(f = 0)";
        case "C1", s = "(0 < f < 0.25)";
        case "C2", s = "(0.25 <= f < 0.50)";
        otherwise, s = "(f >= 0.50)";
    end
end

function writeT(T, outDir, name)
    writetable(T, fullfile(outDir, name + ".csv"));
    writetable(T, fullfile(outDir, name + ".xlsx"));
    fprintf("%s: %d rows written\n", name, height(T));
end

function saveFig(fig, outDir, name)
    exportgraphics(fig, fullfile(outDir, name + ".png"), "Resolution", 600);
    exportgraphics(fig, fullfile(outDir, name + ".tif"), "Resolution", 600);
    exportgraphics(fig, fullfile(outDir, name + ".pdf"), "ContentType", "vector");
    savefig(fig, fullfile(outDir, name + ".fig"));
    fprintf("%s saved (PNG 600 dpi, TIFF 600 dpi, PDF vector, FIG)\n", name);
    close(fig);
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
