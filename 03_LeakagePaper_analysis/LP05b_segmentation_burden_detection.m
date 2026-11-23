function LP05b_segmentation_burden_detection()
%LP05B_SEGMENTATION_BURDEN_DETECTION  Leakage paper - study and lesion level evaluation.
%
% Locked methods v1.0 (Sections 10.2 items 3-5, 10.3 C1 acceptance, 10.4,
% 10.6, 11) + Amendment 01. Frozen inputs only; no retraining, no
% preprocessing or threshold change, no mechanistic exploration.
%
% Models applied to Pipeline B candidates at their LOCKED thresholds:
%   B, A-std, A-hist (main), A-hist ablation (supplementary, LP04 model).
% Descriptive: A-hist -> A-hist "reference-informed reconstruction".
%
% Per study: predicted mask = union of accepted candidates (whole candidates).
%   Segmentation (tumour-positive studies): Dice (headline, study-weighted
%     mean), voxel precision, recall, Jaccard.
%   Burden (all studies, zeros included): Lin CCC for total tumour volume
%     (headline); bias, Bland-Altman limits, MAE, RMSE, Pearson, Spearman for
%     volume, SUVmean (studies with both masks non-empty) and TL-PSMA.
%   Tumour-negative studies: false-positive volume, studies with any FP.
% Per reference lesion (26-connected GT components):
%   detected = >= 1 voxel overlap with predicted mask; recall for lesions
%   >= 1 mL (headline), all lesions, lesions reachable by any B candidate;
%   misses = candidate-generation miss (no B candidate overlaps) or
%   classifier miss (overlapping candidates all rejected); lesions detected
%   exclusively through accepted C1 candidates.
% Candidate level (supplementary ablation): AUC, sensitivity, specificity,
%   paired differences vs A-hist and vs B.
% Uncertainty: patient-clustered paired bootstrap, 2,000 resamples, seed
%   20261109; one set of patient draws per cohort shared by all study-level,
%   lesion-level and candidate-level metrics (paired across models).
%
% SAFEGUARDS (stop with error on failure):
%   T1 frozen-input hashes (LP03, LP04, four frozen LP05a score files)
%   T2 recomputed scores and accept flags of B/A-std/A-hist equal the frozen
%      LP05a values; ablation scored with the frozen LP04 model
%   T3a spatial alignment: candidate label map, reference mask and
%      reconstruction map vs PET, all eight grid corners in world space
%      (tolerance 1% of the smallest voxel side); spacing differences logged
%   T3b exact voxel-count agreement: every candidate present; candidate,
%      reference-overlap and reconstruction voxel counts equal the frozen
%      volumes divided by the PET voxel volume (relative error < 1e-9)
%   T3c PET scale check (all B candidate voxels have SUV >= 3)
%   All volumes use the PET-grid voxel volume (grid of candidate generation).
%   T4 lesion counts: reference components >= 1 mL equal the frozen
%      reference-derived positive counts (internal test 1,577; external 3,250).
%      Count agreement with identical reference-mask files; not a voxelwise
%      lesion-by-lesion match to V1 components.
% Revision: first run (LP05b_20261109_220050) aborted at T3 because the
% voxel volume was taken from the reference-mask header (z spacing differs
% from PET by ~4e-5 mm in some internal studies). Corrected to PET grid.
%   T5 study counts: tumour-positive 78 internal / 95 external; tumour-negative 9
%   T6 bootstrap resamples patients (multiplicity verified every replicate)

clc;
t0 = tic;
V2 = string(fullfile(psma_root("PSMA_PROJECT_ROOT"),"PSMA_AI_Project_V2_PET_Derived"));
MD = fullfile(V2,"06_results","LeakagePaper_ModelData");
EV = fullfile(V2,"06_results","LeakagePaper_Evaluation");
R9 = fullfile(V2,"06_results","Objective9_PET_Derived_Candidate_Pipeline");
RA = fullfile(V2,"06_results","Objective9_LockedExternal_Addendum_train0096_0100");
L5 = fullfile(EV,"LP05a_20261109_210731");
files = struct( ...
  "lp03", fullfile(MD,"LP03_20261109_203442","LP03_candidate_sets.mat"), ...
  "lp04", fullfile(MD,"LP04_20261109_204900","LP04_models_and_thresholds.mat"), ...
  "scInt", fullfile(L5,"LP05a_scores_internal_test_B.csv"), ...
  "scExt", fullfile(L5,"LP05a_scores_external_B.csv"), ...
  "rcInt", fullfile(L5,"LP05a_scores_Ahist_on_Ahist_internal.csv"), ...
  "rcExt", fullfile(L5,"LP05a_scores_Ahist_on_Ahist_external.csv"));
hashes = struct( ...
  "lp03", "6ce98177fb7dffda4804d1a0717849be331571dfe252b842514b0f15583f3eb2", ...
  "lp04", "bc87e01862d59086d3134aa73e61382276aaec7bce4dcd584e109da5a209a4d6", ...
  "scInt", "7c8932007fa45826fbbb3847a339f43f3a22ff273d753310d4231448f9440470", ...
  "scExt", "b33bb7bc80e73c6621532ad7d4f4f57c34bdfbba39f70bda208c78dd44d5c121", ...
  "rcInt", "4fd366e8461a17225d295ab00a8f9d42d25720e1f1a27da9331d921efe7780fa", ...
  "rcExt", "42bd682ad24611dbbef4df7ee6705cd07f57ded659ed442c80bdf16a5ac17a7e");
petSrc = [fullfile(R9,"Objective9F_MASTER_Candidate_Dataset","internal","PSMA_objective9F_internal_MASTER_candidate_dataset.csv"), ...
          fullfile(R9,"Objective9B_GT_Overlap_and_Label_Assignment","external","external","PSMA_objective9B_candidate_overlap_table_20260716_090139.csv"), ...
          fullfile(RA,"Objective9B_GT_Overlap_and_Label_Assignment","external","external","PSMA_objective9B_candidate_overlap_table_20261109_150612.csv")];

done = dir(fullfile(EV,"LP05b_2*","LP05b_results.mat"));
if ~isempty(done)
    error("A completed LP05b run already exists (%s). LP05b is run once.", done(1).folder);
end
outDir = fullfile(EV, "LP05b_" + string(datetime("now","Format","yyyyMMdd_HHmmss")));
mkdir(outDir);
diary(char(fullfile(outDir,"LP05b_log.txt"))); cleanupDiary = onCleanup(@() diary("off"));
fprintf("LP05b started %s\nMATLAB %s\nOutput: %s\n\n", string(datetime("now")), version, outDir);

cfg.nBoot = 2000; cfg.seed = 20261109;
cfg.models = ["B","Astd","Ahist","Ahist_ablation"];
cfg.segModels = [cfg.models, "Ahist_recon"];
cfg.pairs = ["B","Ahist","key_secondary"; ...
             "B","Astd","mechanistic_sequential_descriptive_not_causal"; ...
             "Astd","Ahist","mechanistic_sequential_descriptive_not_causal"; ...
             "Ahist_ablation","Ahist","supplementary_ablation"; ...
             "B","Ahist_ablation","supplementary_ablation"];

%% T1 ----------------------------------------------------------------------
for f = string(fieldnames(files))'
    h = fileSha256(files.(f));
    if h ~= hashes.(f), error("T1 FAILED: %s hash mismatch (%s).", f, files.(f)); end
end
fprintf("T1  OK: LP03, LP04 and four frozen LP05a score files match recorded SHA256.\n");
D = load(files.lp03, "data"); data = D.data;
M = load(files.lp04, "models"); models = M.models;

% PET paths per study
petMap = containers.Map;
for f = petSrc
    o = detectImportOptions(f, "FileType","text","Delimiter",",","VariableNamingRule","preserve");
    nm = erase(string(o.VariableNames), [string(char(65279)), string(char([239 187 191]))]);
    o.VariableNames = cellstr(nm);
    o.SelectedVariableNames = {'case_id','PET_path'};
    o = setvartype(o, {'case_id','PET_path'}, "string");
    T = readtable(f, o);
    [u, ia] = unique(T.case_id);
    for i = 1:numel(u), petMap(char(u(i))) = char(T.PET_path(ia(i))); end
end

allEst = table(); geoAll = table(); stAll = table(); lesAll = table(); fdAll = table(); c1All = table(); tnAll = table(); sg = strings(0,1);
cohorts = ["internal","external"];
for c = cohorts
    fprintf("\n==================== %s cohort ====================\n", upper(c));
    %% candidates, frozen scores, T2 -------------------------------------
    if c == "internal"
        X = data.internal_B(data.internal_B.split=="test", :);
        F = readtable(files.scInt, "VariableNamingRule","preserve", "TextType","string");
        XA = data.internal_Ahist(data.internal_Ahist.split=="test", :);
        FA = readtable(files.rcInt, "VariableNamingRule","preserve", "TextType","string");
    else
        X = data.external_B;
        F = readtable(files.scExt, "VariableNamingRule","preserve", "TextType","string");
        XA = data.external_Ahist;
        FA = readtable(files.rcExt, "VariableNamingRule","preserve", "TextType","string");
    end
    if ~isequal(string(F.row_id), X.row_id), error("T2 FAILED: frozen score order differs from LP03 (%s).", c); end
    for m = ["B","Astd","Ahist"]
        s = scoreTable(X, models.(m));
        if max(abs(s - F.("s_" + m))) > 1e-9 || any((s >= models.(m).threshold_raw) ~= logical(F.("accept_" + m)))
            error("T2 FAILED: %s scores/accept flags differ from frozen LP05a (%s).", m, c);
        end
    end
    C = X(:, ["row_id","case_id","patient_id","candidate_label_value","candidate_labelmap_path","GT_path", ...
              "label_any","occupancy_category","candidate_volume_mL","overlap_volume_mL"]);
    for m = ["B","Astd","Ahist"]
        C.("s_" + m) = F.("s_" + m); C.("acc_" + m) = logical(F.("accept_" + m));
    end
    C.s_Ahist_ablation = scoreTable(X, models.Ahist_ablation);
    C.acc_Ahist_ablation = C.s_Ahist_ablation >= models.Ahist_ablation.threshold_raw;
    [tf, loc] = ismember(string(FA.row_id), XA.row_id);
    if ~all(tf) || height(FA) ~= height(XA), error("T2 FAILED: reconstruction rows (%s).", c); end
    CA = XA(loc, ["row_id","case_id","patient_id","candidate_label_value","candidate_labelmap_path","candidate_volume_mL"]);
    CA.acc = logical(FA.accept_Ahist);
    fprintf("T2  OK: B/A-std/A-hist scores and accept flags equal frozen LP05a; ablation scored with frozen LP04 model.\n");

    %% per-study loop ---------------------------------------------------
    studies = unique(C.case_id, "stable");
    nS = numel(studies); nM = numel(cfg.segModels);
    Z = @() NaN(nS, nM);
    st = struct("refVol",NaN(nS,1),"refSUVmean",NaN(nS,1),"refTLP",NaN(nS,1),"tumourPos",false(nS,1), ...
        "predVol",Z(),"inter",Z(),"dice",Z(),"prec",Z(),"rec",Z(),"jac",Z(),"predSUVmean",Z(),"predTLP",Z(),"fpVol",Z(), ...
        "nAll",zeros(nS,1),"nGe1",zeros(nS,1),"nReach",zeros(nS,1),"nReachGe1",zeros(nS,1), ...
        "detAll",zeros(nS,nM),"detGe1",zeros(nS,nM),"detReach",zeros(nS,nM),"detReachGe1",zeros(nS,nM), ...
        "c1only",zeros(nS,nM),"c1onlyGe1",zeros(nS,nM),"clfMiss",zeros(nS,nM),"clfMissGe1",zeros(nS,nM), ...
        "genMiss",zeros(nS,1),"genMissGe1",zeros(nS,1));
    patient = strings(nS,1);
    lesRows = cell(nS,1);
    geoRows = cell(nS,1);
    for i = 1:nS
        sid = studies(i);
        k = C.case_id == sid; Ci = C(k,:);
        patient(i) = Ci.patient_id(1);
        lp = Ci.candidate_labelmap_path(1); gp = Ci.GT_path(1);
        if ~isKey(petMap, char(sid)), error("PET path not found for %s", sid); end
        pp = string(petMap(char(sid)));
        iL = niftiinfo(lp); iG = niftiinfo(gp); iP = niftiinfo(pp);
        L = double(niftiread(iL)); G = niftiread(iG) > 0; PET = readScaled(iP);
        if ~isequal(size(L), size(G), size(PET)), error("Grid mismatch for %s", sid); end
        % Voxel volume from the PET grid (grid of candidate generation and of
        % every frozen volume in 9B/9C/LP03).
        vv = prod(double(iP.PixelDimensions(1:3))) / 1000;
        % T3a spatial alignment: label map and reference mask vs PET, all
        % eight grid corners in world space, tolerance 1% of the smallest voxel side
        [dL, sL] = geomCompare(iP, iL, sid, "candidate label map");
        [dG, sG] = geomCompare(iP, iG, sid, "reference mask");
        % T3b exact voxel-count agreement with frozen volumes (PET voxel volume)
        mx = max([L(:); Ci.candidate_label_value]);
        cnt = accumarray(L(L>0), 1, [mx 1]);
        ovl = accumarray(L(L>0 & G), 1, [mx 1]);
        lv = Ci.candidate_label_value;
        if any(cnt(lv) == 0), error("T3 FAILED: candidate label missing in map (%s).", sid); end
        if ~isequal(cnt(lv), round(Ci.candidate_volume_mL / vv)) || ...
                max(abs(cnt(lv)*vv - Ci.candidate_volume_mL) ./ Ci.candidate_volume_mL) > 1e-9
            error("T3 FAILED: candidate voxel counts/volumes differ from frozen values (%s).", sid);
        end
        if ~isequal(ovl(lv), round(Ci.overlap_volume_mL / vv)) || max(abs(ovl(lv)*vv - Ci.overlap_volume_mL)) > 1e-9
            error("T3 FAILED: reference-overlap voxel counts/volumes differ from frozen values (%s).", sid);
        end
        % T3c PET scale
        if min(PET(L>0)) < 3 - 1e-4
            error("T3 FAILED: PET scale inconsistent with SUV >= 3 candidates (%s).", sid);
        end
        % reference
        st.refVol(i) = nnz(G) * vv;
        st.tumourPos(i) = any(G(:));
        if st.tumourPos(i), st.refSUVmean(i) = mean(PET(G)); st.refTLP(i) = sum(PET(G)) * vv; else, st.refTLP(i) = 0; end
        % reconstruction map
        CAi = CA(CA.case_id == sid, :);
        dA = NaN; sA = NaN;
        if isempty(CAi)
            LA = zeros(size(G));          % no reference-derived candidates in this study
        else
            iA = niftiinfo(CAi.candidate_labelmap_path(1));
            LA = double(niftiread(iA));
            if ~isequal(size(LA), size(G)), error("Reconstruction map grid mismatch (%s).", sid); end
            [dA, sA] = geomCompare(iP, iA, sid, "reconstruction label map");
            mxA = max([LA(:); CAi.candidate_label_value]);
            cntA = accumarray(LA(LA>0), 1, [mxA 1]);
            if ~isequal(cntA(CAi.candidate_label_value), round(CAi.candidate_volume_mL / vv)) || ...
                    max(abs(cntA(CAi.candidate_label_value)*vv - CAi.candidate_volume_mL) ./ CAi.candidate_volume_mL) > 1e-9
                error("T3 FAILED: reconstruction candidate voxel counts/volumes differ (%s).", sid);
            end
        end
        geoRows{i} = table(c, sid, sL, sG, sA, dL, dG, dA, ...
            'VariableNames', {'cohort','case_id','max_spacing_diff_mm_labelmap_vs_PET','max_spacing_diff_mm_reference_vs_PET', ...
            'max_spacing_diff_mm_recon_vs_PET','max_corner_disp_vox_labelmap_vs_PET','max_corner_disp_vox_reference_vs_PET', ...
            'max_corner_disp_vox_recon_vs_PET'});
        % predicted masks
        Pm = cell(1, nM);
        for j = 1:nM
            mname = cfg.segModels(j);
            if mname == "Ahist_recon"
                Pm{j} = ismember(LA, CAi.candidate_label_value(CAi.acc));
            else
                Pm{j} = ismember(L, lv(Ci.("acc_" + mname)));
            end
            P = Pm{j};
            nP = nnz(P); nI = nnz(P & G); nR = nnz(G);
            st.predVol(i,j) = nP*vv; st.inter(i,j) = nI*vv; st.fpVol(i,j) = (nP - nI)*vv;
            st.predTLP(i,j) = sum(PET(P))*vv;
            if nP > 0, st.predSUVmean(i,j) = mean(PET(P)); st.prec(i,j) = nI/nP; end
            if nR > 0
                st.dice(i,j) = 2*nI/(nP + nR); st.rec(i,j) = nI/nR; st.jac(i,j) = nI/(nP + nR - nI);
            end
        end
        % lesions
        cc = bwconncomp(G, 26);
        nl = cc.NumObjects;
        lr = table();
        if nl > 0
            occOf = containers.Map(num2cell(lv), cellstr(Ci.occupancy_category));
            lvol = cellfun(@numel, cc.PixelIdxList)' * vv;
            reach = false(nl,1); det = false(nl,nM); c1o = false(nl,nM-1);
            for q = 1:nl
                vox = cc.PixelIdxList{q};
                labs = unique(L(vox)); labs(labs==0) = [];
                reach(q) = ~isempty(labs);
                for j = 1:nM
                    det(q,j) = any(Pm{j}(vox));
                    if j < nM && det(q,j)
                        accL = labs(ismember(labs, lv(Ci.("acc_" + cfg.segModels(j)))));
                        occ = string(values(occOf, num2cell(accL)));
                        c1o(q,j) = ~isempty(accL) && all(occ == "C1");
                    end
                end
            end
            ge1 = lvol >= 1.0;
            st.nAll(i) = nl; st.nGe1(i) = sum(ge1); st.nReach(i) = sum(reach); st.nReachGe1(i) = sum(reach & ge1);
            st.genMiss(i) = sum(~reach); st.genMissGe1(i) = sum(~reach & ge1);
            for j = 1:nM
                st.detAll(i,j) = sum(det(:,j)); st.detGe1(i,j) = sum(det(:,j) & ge1);
                st.detReach(i,j) = sum(det(:,j) & reach); st.detReachGe1(i,j) = sum(det(:,j) & reach & ge1);
                if j < nM
                    st.c1only(i,j) = sum(c1o(:,j)); st.c1onlyGe1(i,j) = sum(c1o(:,j) & ge1);
                    st.clfMiss(i,j) = sum(reach & ~det(:,j)); st.clfMissGe1(i,j) = sum(reach & ~det(:,j) & ge1);
                end
            end
            lr = table(repmat(sid,nl,1), (1:nl)', lvol, ge1, reach, ...
                'VariableNames', {'case_id','lesion_id','volume_mL','ge_1mL','reachable_by_B_candidate'});
            for j = 1:nM
                lr.("detected_" + cfg.segModels(j)) = det(:,j);
                if j < nM
                    lr.("C1only_" + cfg.segModels(j)) = c1o(:,j);
                    mt = repmat("detected", nl,1); mt(~det(:,j) & ~reach) = "candidate_generation_miss";
                    mt(~det(:,j) & reach) = "classifier_miss";
                    lr.("outcome_" + cfg.segModels(j)) = mt;
                end
            end
        end
        lesRows{i} = lr;
        if mod(i,10) == 0 || i == nS, fprintf("  studies processed: %d / %d\n", i, nS); end
    end
    lesRows = lesRows(~cellfun(@isempty, lesRows));
    les = vertcat(lesRows{:});
    geo = vertcat(geoRows{:});
    geoAll = [geoAll; geo]; %#ok<AGROW>
    fprintf("T3a OK: spatial alignment of label maps and reference masks with PET verified in %d studies (max corner displacement %.2g voxel; %d studies with a nonzero reference-vs-PET spacing difference, max %.2g mm)\n", ...
        nS, max([geo.max_corner_disp_vox_labelmap_vs_PET; geo.max_corner_disp_vox_reference_vs_PET; geo.max_corner_disp_vox_recon_vs_PET]), ...
        sum(geo.max_spacing_diff_mm_reference_vs_PET > 0), max(geo.max_spacing_diff_mm_reference_vs_PET));
    fprintf("T3b OK: exact voxel-count agreement with frozen candidate, overlap and reconstruction volumes (PET voxel volume)\nT3c OK: PET scale\n");

    %% T4 / T5 --------------------------------------------------------------
    expGe1 = struct("internal",1577,"external",3250);
    expPos = struct("internal",78,"external",95);
    expNeg = struct("internal",9,"external",0);
    chk(sprintf("T4 reference lesions >= 1 mL (%s)", c), sum(st.nGe1), expGe1.(c));
    chk(sprintf("T5 tumour-positive studies (%s)", c), sum(st.tumourPos), expPos.(c));
    chk(sprintf("T5 tumour-negative studies (%s)", c), sum(~st.tumourPos), expNeg.(c));

    %% study table output ---------------------------------------------------
    for j = 1:nM
        stAll = [stAll; table(repmat(c,nS,1), studies, patient, repmat(cfg.segModels(j),nS,1), st.tumourPos, ...
            st.refVol, st.predVol(:,j), st.inter(:,j), st.fpVol(:,j), st.dice(:,j), st.prec(:,j), st.rec(:,j), st.jac(:,j), ...
            st.refSUVmean, st.predSUVmean(:,j), st.refTLP, st.predTLP(:,j), st.nAll, st.detAll(:,j), st.nGe1, st.detGe1(:,j), ...
            'VariableNames', {'cohort','case_id','patient_id','model','tumour_positive','ref_volume_mL','pred_volume_mL', ...
            'intersection_mL','FP_volume_mL','Dice','voxel_precision','voxel_recall','Jaccard','ref_SUVmean','pred_SUVmean', ...
            'ref_TLPSMA','pred_TLPSMA','n_lesions','n_detected','n_lesions_ge1mL','n_detected_ge1mL'})]; %#ok<AGROW>
    end
    if ~isempty(les), les.cohort = repmat(c, height(les), 1); lesAll = [lesAll; les]; end %#ok<AGROW>

    %% failure decomposition, C1-only, tumour-negative (point) --------------
    for j = 1:nM-1
        m = cfg.segModels(j);
        fdAll = [fdAll; table(c, m, sum(st.nAll), sum(st.detAll(:,j)), sum(st.genMiss), sum(st.clfMiss(:,j)), ...
            sum(st.nGe1), sum(st.detGe1(:,j)), sum(st.genMissGe1), sum(st.clfMissGe1(:,j)), ...
            'VariableNames', {'cohort','model','lesions_all','detected_all','candidate_generation_miss_all','classifier_miss_all', ...
            'lesions_ge1mL','detected_ge1mL','candidate_generation_miss_ge1mL','classifier_miss_ge1mL'})]; %#ok<AGROW>
        c1All = [c1All; table(c, m, sum(st.detAll(:,j)), sum(st.c1only(:,j)), sum(st.detGe1(:,j)), sum(st.c1onlyGe1(:,j)), ...
            'VariableNames', {'cohort','model','lesions_detected_all','detected_exclusively_via_C1_all', ...
            'lesions_detected_ge1mL','detected_exclusively_via_C1_ge1mL'})]; %#ok<AGROW>
    end
    neg = ~st.tumourPos;
    if any(neg)
        for j = 1:nM
            tnAll = [tnAll; table(c, cfg.segModels(j), sum(neg), sum(st.predVol(neg,j) > 0), sum(st.predVol(neg,j)), median(st.predVol(neg,j)), ...
                'VariableNames', {'cohort','model','n_tumour_negative_studies','n_with_any_FP','total_FP_volume_mL','median_FP_volume_mL'})]; %#ok<AGROW>
        end
    end

    %% bootstrap (T6) -------------------------------------------------------
    [pList, ~, pidS] = unique(patient);
    [~, pidC] = ismember(C.patient_id, pList);
    if any(pidC == 0), error("Candidate patient not among study patients."); end
    fun = @(sIdx, cIdx) allMetrics(st, C, sIdx, cIdx, cfg);
    [pt, names] = fun((1:nS)', (1:height(C))');
    [B, nRej] = patientBootstrap(numel(pList), pidS, pidC, fun, st, C, cfg);
    est = summarise(c, names, pt, B, nRej, cfg);
    allEst = [allEst; est]; %#ok<AGROW>
    writetable(array2table(B, 'VariableNames', cellstr(names)), fullfile(outDir, "LP05b_bootstrap_replicates_" + c + ".csv"));
    printKey(est, c);
    sg(end+1,1) = sprintf("%s: T2 PASS; T3a alignment PASS, T3b exact voxel counts PASS, T3c PET scale PASS (%d studies); T4 PASS (%d lesions >= 1 mL; count agreement with identical reference-mask files, not a voxelwise lesion-by-lesion match to V1); T5 PASS; T6 PASS (%d patients)", ...
        c, nS, sum(st.nGe1), numel(pList)); %#ok<AGROW>
end

disp(fdAll); disp(c1All); disp(tnAll);
writetable(allEst, fullfile(outDir,"LP05b_estimates_with_bootstrap_CI.csv"));
writetable(stAll, fullfile(outDir,"LP05b_study_level_results.csv"));
writetable(geoAll, fullfile(outDir,"LP05b_geometry_alignment_check.csv"));
writetable(lesAll, fullfile(outDir,"LP05b_lesion_level_results.csv"));
writetable(fdAll, fullfile(outDir,"LP05b_failure_decomposition.csv"));
writetable(c1All, fullfile(outDir,"LP05b_C1_only_detection.csv"));
writetable(tnAll, fullfile(outDir,"LP05b_tumour_negative_studies.csv"));
writetable(table(["T1";sg], 'VariableNames', {'safeguard_summary'}), fullfile(outDir,"LP05b_safeguards.csv"));
save(fullfile(outDir,"LP05b_results.mat"), "allEst","geoAll","stAll","lesAll","fdAll","c1All","tnAll","cfg","files","hashes","-v7.3");
fprintf("\nLP05b COMPLETED (%.1f min). All safeguards passed. No model, preprocessing or threshold was changed.\nResults: %s\n", toc(t0)/60, outDir);
end

%% ======================================================================
function [v, names] = allMetrics(st, C, sIdx, cIdx, cfg)
    v = []; names = strings(1,0);
    segM = cfg.segModels; nM = numel(segM);
    tp = st.tumourPos(sIdx);
    rv = st.refVol(sIdx); rT = st.refTLP(sIdx); rS = st.refSUVmean(sIdx);
    val = struct();
    for j = 1:nM
        m = segM(j);
        d = st.dice(sIdx,j); pr = st.prec(sIdx,j); re = st.rec(sIdx,j); ja = st.jac(sIdx,j);
        pv = st.predVol(sIdx,j); pT = st.predTLP(sIdx,j); pS = st.predSUVmean(sIdx,j);
        x = struct();
        x.Dice = mean(d(tp)); x.precision = mean(pr(tp), "omitnan"); x.recall = mean(re(tp)); x.Jaccard = mean(ja(tp));
        x.CCC_volume = lin(rv, pv);
        [x.bias_volume, x.LoAlow_volume, x.LoAhigh_volume, x.MAE_volume, x.RMSE_volume, x.Pearson_volume, x.Spearman_volume] = agree(rv, pv);
        x.CCC_TLPSMA = lin(rT, pT);
        [x.bias_TLPSMA, x.LoAlow_TLPSMA, x.LoAhigh_TLPSMA, x.MAE_TLPSMA, x.RMSE_TLPSMA, x.Pearson_TLPSMA, x.Spearman_TLPSMA] = agree(rT, pT);
        ok = isfinite(rS) & isfinite(pS);
        [x.bias_SUVmean, x.LoAlow_SUVmean, x.LoAhigh_SUVmean, x.MAE_SUVmean, x.RMSE_SUVmean, x.Pearson_SUVmean, x.Spearman_SUVmean] = agree(rS(ok), pS(ok));
        x.n_SUVmean = sum(ok);
        x.recall_ge1mL = sum(st.detGe1(sIdx,j)) / sum(st.nGe1(sIdx));
        x.recall_all = sum(st.detAll(sIdx,j)) / sum(st.nAll(sIdx));
        if m ~= "Ahist_recon"
            x.recall_reachable = sum(st.detReach(sIdx,j)) / sum(st.nReach(sIdx));
            x.recall_reachable_ge1mL = sum(st.detReachGe1(sIdx,j)) / sum(st.nReachGe1(sIdx));
            x.C1only_lesions_ge1mL = sum(st.c1onlyGe1(sIdx,j));
            x.classifier_miss_ge1mL = sum(st.clfMissGe1(sIdx,j));
        end
        neg = ~tp;
        if any(neg), x.FPvol_tumourNeg_total = sum(pv(neg)); x.FPstudies_tumourNeg = sum(pv(neg) > 0);
        else, x.FPvol_tumourNeg_total = NaN; x.FPstudies_tumourNeg = NaN; end
        x.FPvol_total = sum(st.fpVol(sIdx,j));
        f = string(fieldnames(x))';
        for k = f
            v(end+1) = x.(k); names(end+1) = k + "_" + m; %#ok<AGROW>
        end
        val.(m) = x;
    end
    v(end+1) = sum(st.genMissGe1(sIdx)); names(end+1) = "candidate_generation_miss_ge1mL";
    % paired differences
    keyM = ["Dice","CCC_volume","recall_ge1mL","recall_all","precision","recall","Jaccard","CCC_TLPSMA","MAE_volume"];
    for p = 1:size(cfg.pairs,1)
        a = cfg.pairs(p,1); b = cfg.pairs(p,2);
        for km = keyM
            v(end+1) = val.(a).(km) - val.(b).(km); names(end+1) = "d_" + km + "_" + a + "_minus_" + b; %#ok<AGROW>
        end
    end
    % candidate-level ablation (paired, same patient draws)
    y = C.label_any(cIdx);
    auc = struct(); se = struct(); sp = struct();
    for m = cfg.models
        s = C.("s_" + m)(cIdx); acc = C.("acc_" + m)(cIdx);
        auc.(m) = aucMW(y, s);
        se.(m) = sum(acc & y==1)/sum(y==1); sp.(m) = sum(~acc & y==0)/sum(y==0);
    end
    v = [v, auc.Ahist_ablation, se.Ahist_ablation, sp.Ahist_ablation, ...
        auc.Ahist_ablation - auc.Ahist, auc.B - auc.Ahist_ablation, ...
        se.Ahist_ablation - se.Ahist, se.B - se.Ahist_ablation, sp.Ahist_ablation - sp.Ahist, sp.B - sp.Ahist_ablation];
    names = [names, "candAUC_Ahist_ablation", "candSens_Ahist_ablation", "candSpec_Ahist_ablation", ...
        "d_candAUC_Ahist_ablation_minus_Ahist", "d_candAUC_B_minus_Ahist_ablation", ...
        "d_candSens_Ahist_ablation_minus_Ahist", "d_candSens_B_minus_Ahist_ablation", ...
        "d_candSpec_Ahist_ablation_minus_Ahist", "d_candSpec_B_minus_Ahist_ablation"];
end

function c = lin(x, y)
    x = x(:); y = y(:);
    if numel(x) < 2, c = NaN; return; end
    mx = mean(x); my = mean(y);
    sxy = mean((x-mx).*(y-my)); sx = mean((x-mx).^2); sy = mean((y-my).^2);
    c = 2*sxy / (sx + sy + (mx-my)^2);
end

function [bias, lo, hi, mae, rmse, pr, spr] = agree(x, y)
    x = x(:); y = y(:); d = y - x;
    if numel(d) < 3, [bias,lo,hi,mae,rmse,pr,spr] = deal(NaN); return; end
    bias = mean(d); sd = std(d); lo = bias - 1.96*sd; hi = bias + 1.96*sd;
    mae = mean(abs(d)); rmse = sqrt(mean(d.^2));
    pr = corr(x, y, "Type","Pearson"); spr = corr(x, y, "Type","Spearman");
end

function a = aucMW(y, s)
    y = double(y(:)); s = double(s(:));
    n1 = sum(y==1); n0 = sum(y==0);
    if n1 == 0 || n0 == 0, a = NaN; return; end
    r = tiedrank(s);
    a = (sum(r(y==1)) - n1*(n1+1)/2) / (n1*n0);
end

function [B, nRej] = patientBootstrap(nP, pidS, pidC, fun, st, C, cfg)
    sByP = accumarray(pidS, (1:numel(pidS))', [nP 1], @(x) {x});
    cByP = accumarray(pidC, (1:numel(pidC))', [nP 1], @(x) {x});
    rng(cfg.seed, "twister");
    [~, names] = fun((1:numel(pidS))', (1:numel(pidC))');
    B = NaN(cfg.nBoot, numel(names)); nRej = 0; b = 0;
    while b < cfg.nBoot
        d = randi(nP, nP, 1);
        sIdx = vertcat(sByP{d}); cIdx = vertcat(cByP{d});
        pc = accumarray(d, 1, [nP 1]);
        if ~isequal(accumarray(sIdx,1,[numel(pidS) 1]), pc(pidS)) || ~isequal(accumarray(cIdx,1,[numel(pidC) 1]), pc(pidC))
            error("T6 FAILED: bootstrap replicate is not patient-clustered.");
        end
        y = C.label_any(cIdx);
        if ~(any(y==1) && any(y==0)) || ~any(st.tumourPos(sIdx)) || sum(st.nGe1(sIdx)) == 0
            nRej = nRej + 1; continue;
        end
        b = b + 1;
        B(b,:) = fun(sIdx, cIdx);
    end
    fprintf("T6  OK: %d replicates, resampling unit = PATIENT (%d patients); multiplicity of studies and candidates verified every replicate; %d redrawn\n", ...
        cfg.nBoot, nP, nRej);
end

function est = summarise(cohort, names, pt, B, nRej, cfg)
    k = numel(names);
    lo = NaN(k,1); hi = NaN(k,1); p = NaN(k,1); nv = zeros(k,1); fam = strings(k,1);
    for j = 1:k
        x = B(:,j); x = x(isfinite(x)); nv(j) = numel(x);
        if ~isempty(x), qq = quantile(x, [0.025 0.975]); lo(j) = qq(1); hi(j) = qq(2); end
        nm = names(j);
        if startsWith(nm, "d_") && ~isempty(x), p(j) = min(1, 2*min(mean(x <= 0), mean(x >= 0))); end
        fam(j) = "estimate";
        for q = 1:size(cfg.pairs,1)
            if contains(nm, "_" + cfg.pairs(q,1) + "_minus_" + cfg.pairs(q,2)), fam(j) = cfg.pairs(q,3); end
        end
        if contains(nm, "Ahist_recon"), fam(j) = "descriptive_reference_informed_reconstruction"; end
        if startsWith(nm, ["candAUC","candSens","candSpec","d_cand"]), fam(j) = "supplementary_ablation_candidate_level"; end
        if any(startsWith(nm, ["Dice_","CCC_volume_","recall_ge1mL_"])) && fam(j) == "estimate", fam(j) = "headline_estimate"; end
    end
    est = table(repmat(string(cohort),k,1), names(:), fam, pt(:), lo, hi, p, nv, repmat(nRej,k,1), ...
        'VariableNames', {'cohort','metric','family','estimate','ci_low','ci_high','bootstrap_p_two_sided','n_valid_replicates','n_redrawn_replicates'});
end

function printKey(est, c)
    show = ["Dice_B","Dice_Astd","Dice_Ahist","Dice_Ahist_ablation","Dice_Ahist_recon", ...
        "d_Dice_B_minus_Ahist","d_Dice_B_minus_Astd","d_Dice_Astd_minus_Ahist", ...
        "CCC_volume_B","CCC_volume_Astd","CCC_volume_Ahist","CCC_volume_Ahist_ablation","CCC_volume_Ahist_recon","d_CCC_volume_B_minus_Ahist", ...
        "bias_volume_B","bias_volume_Ahist", ...
        "recall_ge1mL_B","recall_ge1mL_Astd","recall_ge1mL_Ahist","recall_ge1mL_Ahist_ablation","recall_ge1mL_Ahist_recon","d_recall_ge1mL_B_minus_Ahist", ...
        "precision_B","precision_Ahist","FPvol_total_B","FPvol_total_Ahist", ...
        "candAUC_Ahist_ablation","d_candAUC_Ahist_ablation_minus_Ahist","d_candAUC_B_minus_Ahist_ablation", ...
        "d_Dice_Ahist_ablation_minus_Ahist","d_Dice_B_minus_Ahist_ablation"];
    fprintf("\n%-40s %12s %12s %12s %8s\n", "metric ("+c+")", "estimate", "CI low", "CI high", "p");
    for nm = show
        r = est(est.metric == nm, :);
        if isempty(r), continue; end
        fprintf("%-40s %12.4f %12.4f %12.4f %8.4f\n", nm, r.estimate, r.ci_low, r.ci_high, r.bootstrap_p_two_sided);
    end
    fprintf("(bootstrap p printed as 0.0000 is reported as p < 0.001)\n");
end

function chk(label, got, expected)
    if got ~= expected, error("CHECK FAILED: %s = %d (expected %d).", label, got, expected); end
    fprintf("%-50s OK: %d\n", label, got);
end

function [dispVox, spDiff] = geomCompare(iRef, iX, sid, label)
% World-space positions of all eight grid corners of iX vs the PET grid.
    if ~isequal(iRef.ImageSize(1:3), iX.ImageSize(1:3)), error("T3a FAILED: %s size differs from PET (%s).", label, sid); end
    d = double(iRef.ImageSize(1:3)) - 1;
    [a,b,c] = ndgrid([0 d(1)], [0 d(2)], [0 d(3)]);
    K = [a(:) b(:) c(:) ones(8,1)];
    wR = K * double(iRef.Transform.T); wX = K * double(iX.Transform.T);
    sp = double(iRef.PixelDimensions(1:3));
    dispVox = max(sqrt(sum((wR(:,1:3) - wX(:,1:3)).^2, 2))) / min(sp);
    spDiff = max(abs(double(iX.PixelDimensions(1:3)) - sp));
    if dispVox > 0.01 || spDiff > 1e-3
        error("T3a FAILED: %s not aligned with PET (%s): corner displacement %.4g voxel, spacing difference %.3g mm.", label, sid, dispVox, spDiff);
    end
end

function A = readScaled(info)
    A = double(niftiread(info));
    sl = double(info.MultiplicativeScaling); of = double(info.AdditiveOffset);
    if sl ~= 0 && ~(sl == 1 && of == 0), A = A*sl + of; end
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
