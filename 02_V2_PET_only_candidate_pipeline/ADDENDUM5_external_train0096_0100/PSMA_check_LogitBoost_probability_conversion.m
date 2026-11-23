function PSMA_check_LogitBoost_probability_conversion()
% READ-ONLY CHECK of how the locked V2 LogitBoost model's scores relate to
% probabilities. Nothing is written to disk; results print to the Command
% Window only. The locked model and threshold are not changed.
%
% Questions answered:
%   1. What ScoreTransform does the stored model use?
%   2. What do the two score columns returned by predict() look like?
%   3. Do 9H3's stored probabilities equal 1/(1+exp(-s)) of the raw score?
%   4. Does 'doublelogit' (1/(1+exp(-2s))) give columns that sum to 1?
%   5. What are the equivalent thresholds on the raw-score and doublelogit
%      scales for the locked threshold 0.445?

clc;
projectRoot = string(fullfile(psma_root("PSMA_PROJECT_ROOT"),"PSMA_AI_Project_V2_PET_Derived"));
resultsRoot = fullfile(projectRoot,"06_results","Objective9_PET_Derived_Candidate_Pipeline");
lockedMat   = fullfile(resultsRoot,"Objective9H_Model_Development", ...
    "Objective9H2_Locked_Validation_Threshold","PSMA_objective9H2_LOCKED_model_and_threshold.mat");
testCsv     = fullfile(resultsRoot,"Objective9G_Feature_Audit_and_Preparation","internal", ...
    "PSMA_objective9G_internal_test_model_ready.csv");
h3Csv       = fullfile(resultsRoot,"Objective9H_Model_Development", ...
    "Objective9H3_Single_Untouched_Test_Evaluation","PSMA_objective9H3_test_predictions_LOCKED.csv");

S   = load(lockedMat,"lockedArtifact");
A   = S.lockedArtifact;
mdl = A.model;
pred = string(A.predictors(:));
thr  = double(A.locked_threshold);

fprintf("Model class        : %s\n", class(mdl));
fprintf("Method             : %s\n", string(mdl.Method));
fprintf("ClassNames         : %s\n", strjoin(string(mdl.ClassNames),", "));
fprintf("ScoreTransform     : %s\n", string(mdl.ScoreTransform));
fprintf("Locked threshold   : %.6f\n\n", thr);

% Test predictors, read with a fixed comma delimiter (same as 9H3).
opts = detectImportOptions(testCsv,"FileType","text","Delimiter",",","VariableNamingRule","preserve");
T = readtable(testCsv,opts);
X = double(table2array(T(:,cellstr(pred))));
ids = string(T.candidate_global_id);

[~,s] = predict(mdl,X);                 % native scores
pos = find(mdl.ClassNames==1);
neg = find(mdl.ClassNames==0);
sp  = s(:,pos);
fprintf("Native scores: positive column range [%.4f, %.4f]\n", min(sp), max(sp));
fprintf("               max |s_pos + s_neg| = %.3g  (0 => symmetric margins)\n", max(abs(s(:,pos)+s(:,neg))));
fprintf("               max |s_pos + s_neg - 1| = %.3g  (0 => already probabilities)\n\n", max(abs(sum(s,2)-1)));

% Compare with stored 9H3 probabilities.
H = readtable(h3Csv,detectImportOptions(h3Csv,"FileType","text","Delimiter",",","VariableNamingRule","preserve"));
[tf,loc] = ismember(ids,string(H.candidate_global_id));
fprintf("Rows matched to 9H3 predictions: %d of %d\n", sum(tf), numel(ids));
p9H3 = H.predicted_probability(loc(tf));
pLogit = 1./(1+exp(-sp(tf)));
fprintf("max |9H3 prob - 1/(1+exp(-s))| = %.3g\n", max(abs(p9H3-pLogit)));

% doublelogit via the built-in transform, on a copy of the model.
mdl2 = mdl; mdl2.ScoreTransform = "doublelogit";
[~,s2] = predict(mdl2,X);
fprintf("doublelogit: max |row sum - 1| = %.3g ; range [%.4f, %.4f]\n", max(abs(sum(s2,2)-1)), min(s2(:,pos)), max(s2(:,pos)));
fprintf("doublelogit vs 1/(1+exp(-2s)): max diff = %.3g\n\n", max(abs(s2(:,pos)-1./(1+exp(-2*sp)))));

% Threshold on each scale (thresholds are monotone-equivalent).
rawThr = log(thr/(1-thr));
fprintf("Locked threshold 0.445 on logistic scale  = raw score %.6f\n", rawThr);
fprintf("Same cut on doublelogit probability scale = %.6f\n", 1/(1+exp(-2*rawThr)));
fprintf("Classes identical (logistic vs raw cut): %d\n", isequal(pLogit>=thr, sp(tf)>=rawThr));
fprintf("\nDONE. Nothing was written to disk.\n");
end
