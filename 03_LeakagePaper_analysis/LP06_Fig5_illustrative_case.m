function LP06_Fig5_illustrative_case()
%LP06_FIG5_ILLUSTRATIVE_CASE  Leakage paper - Figure 5, illustrative case (display only).
%
% Illustrative example selected by the authors: internal test study
% psma_ba6433757dfc9d21_2021-10-09, coronal slice 101 (0-based; MATLAB index
% 102; Mango Y = 98). Rows: tumour called by B, A-std and A-hist at their
% LOCKED thresholds (frozen LP05a accept flags), with a zoomed view of the
% liver margin. Values under each row are frozen LP05b per-study results.
% Display windows only; no image values, masks or results are changed or
% computed. Greyscale only.

clc;
V2 = string(fullfile(psma_root("PSMA_PROJECT_ROOT"),"PSMA_AI_Project_V2_PET_Derived"));
V1 = string(fullfile(psma_root("PSMA_PROJECT_ROOT"),"PSMA_AI_Project"));
cid = "psma_ba6433757dfc9d21_2021-10-09";
J = 101 + 1;                                % MATLAB index of raw coronal slice 101
f.ct  = fullfile(V1,"00_data_split_physical","test","imagesTr", cid + "_0000.nii.gz");
f.pet = fullfile(V1,"00_data_split_physical","test","imagesTr", cid + "_0001.nii.gz");
f.gt  = fullfile(V1,"00_data_split_physical","test","labelsTr", cid + ".nii.gz");
f.lm  = fullfile(V2,"06_results","Objective9_PET_Derived_Candidate_Pipeline","Objective9A_PET_Only_Candidate_Generation", ...
                 "internal","test","candidate_label_maps", cid + "_Objective9A_PET_candidate_labelmap.nii.gz");
f.sc  = fullfile(V2,"06_results","LeakagePaper_Evaluation","LP05a_20261109_210731","LP05a_scores_internal_test_B.csv");
f.st  = fullfile(V2,"06_results","LeakagePaper_Evaluation","LP05b_20261109_220537","LP05b_study_level_results.csv");
frozen = struct("sc","7c8932007fa45826fbbb3847a339f43f3a22ff273d753310d4231448f9440470", ...
                "st","75260b1e26cbf9c4589391a93130777452a642572c560d34a8ef00e5c6bcb5da");
for k = ["sc","st"]
    if fileSha256(f.(k)) ~= frozen.(k), error("Frozen input hash mismatch: %s", f.(k)); end
end
outDir = fullfile(V2,"LeakagePaper","Manuscript","LP06_Fig5_" + string(datetime("now","Format","yyyyMMdd_HHmmss")));
mkdir(outDir);
diary(char(fullfile(outDir,"LP06_Fig5_log.txt"))); cleanupDiary = onCleanup(@() diary("off"));
fprintf("Figure 5 (illustrative case) %s\nMATLAB %s\nFrozen LP05a/LP05b inputs verified.\n", string(datetime("now")), version);
for k = ["ct","pet","gt","lm"], fprintf("  %-3s %s  SHA256 %s\n", k, f.(k), fileSha256(f.(k))); end

iP = niftiinfo(f.pet);
CT = double(niftiread(f.ct)); PET = readScaled(iP);
GT = niftiread(f.gt) > 0; LM = double(niftiread(f.lm));
if ~isequal(size(CT), size(PET), size(GT), size(LM)), error("Grid mismatch."); end
px = double(iP.PixelDimensions(1:3));

sc = readtable(f.sc, "TextType","string", "VariableNamingRule","preserve"); sc = sc(sc.case_id == cid, :);
st = readtable(f.st, "TextType","string", "VariableNamingRule","preserve"); st = st(st.case_id == cid, :);
% documented facts used in the caption
c12 = sc(sc.candidate_label_value == 12, :);
fprintf("Candidate 12: %.0f mL, tumour %.2f mL, f = %.4f, class %s, accepted B/A-std/A-hist = %d/%d/%d\n", ...
    c12.candidate_volume_mL, c12.overlap_volume_mL, c12.f_overlap, c12.occupancy_category, c12.accept_B, c12.accept_Astd, c12.accept_Ahist);

% coronal slice -> display (superior up; orientation as in Mango)
cor = @(A) flipud(squeeze(A(:, J, :))');
ct = cor(CT); pet = cor(PET); gt = cor(GT); lm = cor(LM);
win = @(a, lo, hi) min(max((a - lo) ./ (hi - lo), 0), 1);
ctSoft = win(ct, 40 - 200, 40 + 200);
petInv = 1 - win(pet, 0, 10);
fused  = ctSoft .* (1 - 0.75 * win(pet, 0, 10));
ctComp = 0.50 + 0.38 * ctSoft;
xd = [0, (size(ct,2)-1)*px(1)]; yd = [0, (size(ct,1)-1)*px(3)];

% zoom box around the liver-margin lesions (display coordinates, voxels)
nz = size(CT,3);
colR = [45 105] + 1;                         % raw i 45-105
rowR = sort(nz - [225 295]);                 % raw k 225-295 -> display rows
rowR = [max(rowR(1),1), min(rowR(2),size(ct,1))];

mods = ["B","Astd","Ahist"]; modLab = ["B (PET-only)","A-std","A-hist"];
fig = figure("Color","w","Units","centimeters","Position",[1 1 19 30]);
t = tiledlayout(fig, 4, 3, "TileSpacing","compact", "Padding","compact");
showG(nexttile(t), petInv, xd, yd, "PET (inverse grey, SUV 0-10)");
showG(nexttile(t), ctSoft, xd, yd, "CT (W 400 / L 40)");
showG(nexttile(t), fused,  xd, yd, "Fused PET/CT");
for m = 1:3
    acc = sc.candidate_label_value(sc.("accept_" + mods(m)) == 1);
    pred = ismember(lm, acc);
    rgb = overlayRGB(ctComp, pred & gt, pred & ~gt, gt & ~pred);
    ax = nexttile(t);
    image(ax, xd, yd, rgb); axis(ax,"image"); axis(ax,"off");
    rectangle(ax, "Position",[(colR(1)-1)*px(1), (rowR(1)-1)*px(3), diff(colR)*px(1), diff(rowR)*px(3)], ...
        "EdgeColor","k", "LineStyle","--", "LineWidth",0.8);
    title(ax, modLab(m), "FontSize",9, "FontWeight","bold");
    az = nexttile(t);
    image(az, [0, diff(colR)*px(1)], [0, diff(rowR)*px(3)], rgb(rowR(1):rowR(2), colR(1):colR(2), :)); axis(az,"image"); axis(az,"off");
    title(az, "Liver margin (zoom)", "FontSize",8);
    at = nexttile(t); axis(at,"off");
    r = st(st.model == mods(m), :);
    text(at, 0.02, 0.5, sprintf("Dice %.2f\nPredicted %s mL\nReference %s mL\nFalse positive %s mL\nLesions >= 1 mL detected %d/%d", ...
        r.Dice, fmtN(r.pred_volume_mL), fmtN(r.ref_volume_mL), fmtN(r.FP_volume_mL), r.n_detected_ge1mL, r.n_lesions_ge1mL), ...
        "FontSize",8, "VerticalAlignment","middle", "Units","normalized");
end
lg = axes(fig, "Position",[0.05 0.0 0.9 0.03], "Visible","off"); hold(lg,"on");
h1 = patch(lg, NaN, NaN, "w", "EdgeColor","k", "LineWidth",1);
h2 = patch(lg, NaN, NaN, [0.25 0.25 0.25], "EdgeColor","w");
h3 = patch(lg, NaN, NaN, "w", "FaceColor","none", "EdgeColor","k", "LineWidth",1.5);
legend(lg, [h1 h2 h3], ["True positive (tumour called tumour)","False positive (non-tumour called tumour)","Missed tumour (outline)"], ...
    "Orientation","horizontal", "Location","south", "Box","off", "FontSize",7);

for ext = [".png", ".tif"], exportgraphics(fig, fullfile(outDir, "Fig5_illustrative_case" + ext), "Resolution", 600); end
exportgraphics(fig, fullfile(outDir, "Fig5_illustrative_case.pdf"), "ContentType", "vector");
savefig(fig, fullfile(outDir, "Fig5_illustrative_case.fig"));

cap = "Figure 5. Illustrative example selected by the authors (internal test study, coronal slice). Top row: PET (inverse grey, SUV 0-10), CT (W 400 / L 40) and fused PET/CT (uptake shown as darkening of CT). " + ...
      "Rows 2-4: tumour called by each model at its locked threshold on the same slice, with a zoomed view of the liver margin (dashed box). White fill with black ring: true positive; dark grey fill with white edge: false positive; black outline: missed reference tumour. " + ...
      "A 2.2-mL lesion at the liver margin connected the liver and bowel into one 4.5-L PET-only candidate (occupancy class C1; tumour fraction < 0.001), which the PET-only model accepted: it detected 15 of 16 lesions >= 1 mL but predicted 4,745 mL against a reference volume of 84 mL (Dice 0.03). " + ...
      "The reference-informed models rejected this candidate and missed most lesions (A-std 3/16; A-hist 0/16). Values are frozen per-study results. Display windows only; image values and masks are unchanged.";
writelines(cap, fullfile(outDir, "Fig5_caption.txt"));
fprintf("Figure 5 saved (PNG/TIFF 600 dpi, PDF, FIG) to %s\n", outDir);
end

%% ======================================================================
function showG(ax, img, xd, yd, ttl)
    image(ax, xd, yd, repmat(img, 1, 1, 3)); axis(ax,"image"); axis(ax,"off");
    title(ax, ttl, "FontSize",8);
end

function rgb = overlayRGB(base, tp, fp, miss)
    g = base;
    g(fp) = 0.25; g(bwperim(fp, 8)) = 1.0;              % FP: dark grey fill, white edge
    g(tp) = 1.0;  g(imdilate(tp, ones(3)) & ~tp) = 0.0; % TP: white fill, black ring outside
    if any(miss(:))
        mr = imdilate(miss, ones(3));
        g(bwperim(mr, 8)) = 0.0;                        % missed: black outline
    end
    rgb = repmat(g, 1, 1, 3);
end

function s = fmtN(x)
    s = string(sprintf('%.0f', x));
    if abs(x) >= 1000, s = string(sprintf('%d,%03d', floor(x/1000), round(mod(x,1000)))); end
end

function A = readScaled(info)
    A = double(niftiread(info));
    sl = double(info.MultiplicativeScaling); of = double(info.AdditiveOffset);
    if sl ~= 0 && ~(sl == 1 && of == 0), A = A*sl + of; end
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
