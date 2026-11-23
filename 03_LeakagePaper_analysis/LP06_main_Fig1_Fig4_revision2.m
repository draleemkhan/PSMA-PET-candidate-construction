function LP06_main_Fig1_Fig4_revision2()
%LP06_MAIN_FIG1_FIG4_REVISION2  Leakage paper - compact redesign of Figure 1 and Figure 4 (revision 2).
%
%   Revision 2 (after review of revision 1 output LP06_main_Fig1_Fig4_rev_*):
%   Fig 1  all rows share one width; the three pipeline boxes are narrowed to fit
%          exactly under the full-width rows (text wrapped, boxes taller);
%          bus connectors replace floating / diagonal arrows.
%   Fig 4  two-line tick labels (model / cohort) replace the group labels that
%          collided with the panel (b) title; legend column widened (no clipping).
%
%   DISPLAY-ONLY revision. No new statistical analysis; identical frozen inputs
%   (SHA256-verified) and identical numbers to LP06_main_display_items.m
%   (official run LP06_main_20261110_071007), which is kept unchanged.
%
%   Fig 1  study design: boxes sized to their text, double-column width,
%          arrows attached to box edges, reduced vertical spacing.
%   Fig 4  clinical consequence, two panels instead of four:
%          (a) accepted volume, grouped stacked bars (internal / external per model),
%              value labels on bars, legend outside the axes;
%          (b) outcome of reference lesions >= 1 mL as horizontal 100% stacked bars.
%          Interpretive headline removed (belongs in the caption).
%
%   Built-in checks: stacked totals reproduce the frozen accepted volumes;
%   lesion-outcome percentages sum to 100 for every bar.

clc;
V2 = string(fullfile(psma_root("PSMA_PROJECT_ROOT"),"PSMA_AI_Project_V2_PET_Derived"));
EV = fullfile(V2,"06_results","LeakagePaper_Evaluation");
f = struct();
f.aAcc  = fullfile(EV,"LP05a_20261109_210731","LP05a_acceptance_by_occupancy.csv");
f.bEst  = fullfile(EV,"LP05b_20261109_220537","LP05b_estimates_with_bootstrap_CI.csv");
f.bFail = fullfile(EV,"LP05b_20261109_220537","LP05b_failure_decomposition.csv");
f.bC1   = fullfile(EV,"LP05b_20261109_220537","LP05b_C1_only_detection.csv");
h = struct( ...
 "aAcc","2f36134662846e11931ecfb5de550589dda7c9e4c4995d1b324d21acf901a4b8", ...
 "bEst","88131b7d4f102c70a4d1acc9bf5085b6d6c6d471b7219ebe3df10b58e596fb89", ...
 "bFail","52980c24443a48bd0d9cd299de069ed2c80356b4b10e141bef5275e9a219c6e4", ...
 "bC1","4b02e63e9878d9ea2374e1365ff80e7317a989a53256be56ef396ad0e9259787");
for k = string(fieldnames(f))'
    if fileSha256(f.(k)) ~= h.(k), error("Input hash mismatch: %s", f.(k)); end
end

outDir = fullfile(V2,"LeakagePaper","Manuscript","LP06_main_Fig1_Fig4_rev2_" + string(datetime("now","Format","yyyyMMdd_HHmmss")));
mkdir(outDir);
diary(char(fullfile(outDir,"LP06_Fig1_Fig4_rev2_log.txt"))); cleanupDiary = onCleanup(@() diary("off"));
fprintf("LP06 Figure 1 and Figure 4 revision 2 (display only) %s\nMATLAB %s\nAll %d input files verified against frozen SHA256.\nOutput: %s\n\n", ...
    string(datetime("now")), version, numel(fieldnames(f)), outDir);

rd = @(p) readtable(p, "TextType","string", "VariableNamingRule","preserve");
AA = rd(f.aAcc); BE = rd(f.bEst); BF = rd(f.bFail); BC = rd(f.bC1);
mods = ["B","Astd","Ahist"]; modLab = ["B (PET-only)","A-std","A-hist"];
coh = ["internal","external"]; cohLab = ["internal","external"];
ge = @(E,c,m) E(E.cohort==c & E.metric==m, :);
FS = 8; FN = "Arial";
cap = strings(0,1);

%% ===================== FIGURE 1 (compact) =====================
T1 = "Same patients and studies (internal: 597 studies, 378 patients; patient-grouped split 423 / 87 / 87)";
TP = [ sprintf("Pipeline B (PET-only)\nCandidates: SUV ≥ 3, ≥ 1 mL, 26-connectivity\nLabel: any reference overlap\nReference read only for labels"), ...
       sprintf("Pipeline A-std\nPositives: reference components ≥ 1 mL\nNegatives: B negatives with no voxel\nwithin 5 mm of reference"), ...
       sprintf("Pipeline A-hist (historical rule)\nPositives: reference components ≥ 1 mL\nNegatives: SUV ≥ 3 minus dilated\nreference, carved fragments ≥ 1 mL") ];
T3 = sprintf("Identical for all pipelines: features (routine + PET/CT radiomics), training-only preprocessing,\nLogitBoost (200 cycles), threshold on own validation set (sensitivity ≥ 0.90, then maximum MCC)");
T4 = sprintf("Evaluation on the SAME Pipeline B candidates (deployment condition)\nInternal test: 2,673 candidates, 87 studies, 58 patients (primary)   |   External locked 95: 4,173 candidates\nCandidate level (AUC, primary ΔAUC B − A-hist) → segmentation, burden, lesion detection");
shP = [0.90 0.80 0.70];

% common width of all rows = widest full-width text + padding
pad = 0.22; gapV = 0.75; gapC = 0.30; marg = 0.12;
m1 = measureText(T1,FS,FN); m3 = measureText(T3,FS,FN); m4 = measureText(T4,FS,FN);
Wc = max([m1(1), m3(1), m4(1)]) + 2*pad;
wP = (Wc - 2*gapC)/3;                          % three pipeline boxes fill exactly the same width
TPw = strings(1,3); mP = zeros(3,2);
for j = 1:3, TPw(j) = wrapText(TP(j), wP - 2*pad, FS, FN); mP(j,:) = measureText(TPw(j),FS,FN); end
h1 = m1(2)+2*pad; hP = max(mP(:,2))+2*pad; h3 = m3(2)+2*pad; h4 = m4(2)+2*pad;
W = Wc + 2*marg;
H = marg + h1 + gapV + hP + gapV + h3 + gapV + h4 + marg;
fprintf("Figure 1 canvas: %.1f x %.1f cm (row width %.1f cm; pipeline boxes %.2f cm each)\n", W, H, Wc, wP);

fig = figure("Color","w","Units","centimeters","Position",[2 2 W H]);
ax = axes(fig,"Units","normalized","Position",[0 0 1 1]); axis(ax,"off"); hold(ax,"on");
xlim(ax,[0 W]); ylim(ax,[0 H]);
y1 = H - marg - h1;            % bottom of row 1
yP = y1 - gapV - hP;           % bottom of pipeline row
y3 = yP - gapV - h3;
y4 = y3 - gapV - h4;
xc = W/2; xL = xc - Wc/2;
xP = xL + wP/2 + (0:2)*(wP + gapC);
drawTB(ax, xc, y1, Wc, h1, T1, 0.96, FS, FN);
for j = 1:3, drawTB(ax, xP(j), yP, wP, hP, TPw(j), shP(j), FS, FN); end
drawTB(ax, xc, y3, Wc, h3, T3, 0.96, FS, FN);
drawTB(ax, xc, y4, Wc, h4, T4, 0.86, FS, FN);
N = @(x,y) deal(x/W, y/H);
lw = 0.6; lc = [0 0 0];
yb1 = y1 - gapV/2;                                   % bus 1: row 1 -> three pipelines
line(ax, [xc xc], [y1 yb1], "Color",lc, "LineWidth",lw);
line(ax, [xP(1) xP(3)], [yb1 yb1], "Color",lc, "LineWidth",lw);
for j = 1:3, [aX,aY] = N([xP(j) xP(j)], [yb1, yP+hP]); annotation(fig,"arrow",aX,aY,"HeadLength",5,"HeadWidth",5,"LineWidth",lw); end
yb2 = yP - gapV/2;                                   % bus 2: three pipelines -> identical box
for j = 1:3, line(ax, [xP(j) xP(j)], [yP yb2], "Color",lc, "LineWidth",lw); end
line(ax, [xP(1) xP(3)], [yb2 yb2], "Color",lc, "LineWidth",lw);
[aX,aY] = N([xc xc], [yb2, y3+h3]); annotation(fig,"arrow",aX,aY,"HeadLength",5,"HeadWidth",5,"LineWidth",lw);
[aX,aY] = N([xc xc], [y3, y4+h4]); annotation(fig,"arrow",aX,aY,"HeadLength",5,"HeadWidth",5,"LineWidth",lw);
saveFig(fig, outDir, "Fig1_study_design");
cap(end+1) = "Figure 1. Study design. Three training candidate sets are built from the same patients; everything after candidate construction is identical. All models are evaluated on the same PET-only (Pipeline B) candidates.";

%% ===================== FIGURE 4 (two panels) =====================
% ---- data (identical sources and arithmetic as LP06_main_display_items.m) ----
V = zeros(6,5); lab6 = strings(6,1); k = 0;     % rows: B int, B ext, Astd int, Astd ext, Ahist int, Ahist ext
Wp = zeros(6,4);
for j = 1:3
    for ci = 1:2
        k = k + 1; c = coh(ci);
        a = AA(AA.cohort==c & AA.model==mods(j),:);
        V(k,1) = a.accepted_tumour_volume_mL(a.occupancy=="all");
        for o = 1:4, V(k,o+1) = a.accepted_nontumour_volume_mL(a.occupancy=="C"+(o-1)); end
        tot = a.accepted_volume_mL(a.occupancy=="all");
        if abs(sum(V(k,:)) - tot) > 1e-6*max(1,tot), error("Volume check failed: %s %s", mods(j), c); end
        fr = BF(BF.cohort==c & BF.model==mods(j),:); c1 = BC(BC.cohort==c & BC.model==mods(j),:);
        Wp(k,:) = [fr.detected_ge1mL - c1.detected_exclusively_via_C1_ge1mL, c1.detected_exclusively_via_C1_ge1mL, ...
                   fr.classifier_miss_ge1mL, fr.candidate_generation_miss_ge1mL] / fr.lesions_ge1mL * 100;
        if abs(sum(Wp(k,:)) - 100) > 1e-9, error("Lesion-outcome check failed: %s %s", mods(j), c); end
        lab6(k) = modLab(j) + ", " + cohLab(ci);
        fprintf("Fig 4 check OK: %-22s accepted %9.1f mL | lesions: %5.1f / %4.1f / %4.1f / %4.1f %%\n", lab6(k), tot, Wp(k,:));
    end
end
nLes = [BF.lesions_ge1mL(BF.cohort=="internal" & BF.model=="B"), BF.lesions_ge1mL(BF.cohort=="external" & BF.model=="B")];

fig = figure("Color","w","Units","centimeters","Position",[2 2 17.5 13]);
t = tiledlayout(fig,2,6,"TileSpacing","loose","Padding","compact");
pos = [1 1.85 3.3 4.15 5.6 6.45];          % grouped positions (internal, external per model)
shV = [0.20 0.97 0.70 0.50 0.85]; ptV = ["none","\","x","-","|"];
% ---- (a) accepted volume ----
axA = nexttile(t,1,[1 4]); hold(axA,"on");
b = bar(axA, pos, V/1000, 0.75, "stacked", "EdgeColor","k");
for q = 1:5, b(q).FaceColor = shV(q)*[1 1 1]; end
xlim(axA,[0.4 7.05]); ymax = max(sum(V,2))/1000; ylim(axA,[0 1.18*ymax]);
hatchBars(axA, b, ptV);
for k = 1:6
    v = sum(V(k,:))/1000;
    if v >= 10, s = sprintf("%.0f L", v); elseif v >= 1, s = sprintf("%.1f L", v); else, s = sprintf("%.2f L", v); end
    text(axA, pos(k), v + 0.025*ymax, s, "HorizontalAlignment","center", "VerticalAlignment","bottom", "FontSize",FS-1, "FontName",FN);
end
tickLab = ["B\newlineint","B\newlineext","A-std\newlineint","A-std\newlineext","A-hist\newlineint","A-hist\newlineext"];
set(axA, "XTick", pos, "XTickLabel", tickLab, "TickLabelInterpreter","tex", "FontSize",FS, "FontName",FN, "TickLength",[0 0]);
ylabel(axA, "Accepted volume, L (sum over studies)", "FontSize",FS);
title(axA, "(a) Accepted volume at the locked threshold, by occupancy class", "FontSize",FS+1, "FontWeight","bold");
box(axA,"off");
axLA = nexttile(t,5,[1 2]); legendAxes(axLA, ["Tumour (all classes)","Non-tumour, C0 (f = 0)","Non-tumour, C1 (0 < f < 0.25)","Non-tumour, C2 (0.25 ≤ f < 0.5)","Non-tumour, C3 (f ≥ 0.5)"], shV, ptV, FS-1, FN);
% ---- (b) lesion outcomes, horizontal 100% bars ----
axB = nexttile(t,7,[1 4]); hold(axB,"on");
shW = [0.25 0.65 0.92 1.00]; ptW = ["none","x","/","."];
bb = barh(axB, pos, Wp, 0.75, "stacked", "EdgeColor","k");
for q = 1:4, bb(q).FaceColor = shW(q)*[1 1 1]; end
xlim(axB,[0 100]); ylim(axB,[0.4 7.05]); set(axB,"YDir","reverse");
hatchBarsH(axB, bb, ptW);
set(axB, "YTick", pos, "YTickLabel", lab6, "FontSize",FS, "FontName",FN, "TickLength",[0 0]);
xlabel(axB, "% of reference lesions ≥ 1 mL", "FontSize",FS);
title(axB, sprintf("(b) Outcome of reference lesions ≥ 1 mL (internal n = %s; external n = %s)", fmtInt(nLes(1)), fmtInt(nLes(2))), "FontSize",FS+1, "FontWeight","bold");
box(axB,"off");
axLB = nexttile(t,11,[1 2]); legendAxes(axLB, ["Detected (not only via C1)","Detected only via accepted C1","Classifier miss","Candidate-generation miss"], shW, ptW, FS-1, FN);
saveFig(fig, outDir, "Fig4_clinical_burden_vs_detection");
bi = ge(BE,"internal","bias_volume_B"); be = ge(BE,"external","bias_volume_B");
cap(end+1) = sprintf("Figure 4. Clinical consequence. (a) Accepted candidate volume at each model's locked threshold (B, PET-only; int, internal test; ext, external), split into reference-tumour volume and non-tumour volume by occupancy class (C0 f = 0; C1 0 < f < 0.25; C2 0.25 ≤ f < 0.5; C3 f ≥ 0.5; f = reference-tumour fraction of the candidate); totals above bars. For Pipeline B, the mean volume bias was %+.0f mL per study internally (95%% CI %+.0f to %+.0f) and %+.0f mL externally (%+.0f to %+.0f), patient-clustered bootstrap. C1 contributes an important part, not necessarily all, of the non-tumour volume. (b) Outcome of reference lesions ≥ 1 mL: detected (≥ 1 voxel overlap with the predicted mask), detected only through accepted C1 candidates, classifier miss, or candidate-generation miss (no Pipeline B candidate overlaps; specific to this frozen candidate set). Better lesion detection does not automatically produce reliable tumour-burden quantification.", ...
    bi.estimate, bi.ci_low, bi.ci_high, be.estimate, be.ci_low, be.ci_high);

writelines(cap, fullfile(outDir,"LP06_Fig1_Fig4_rev2_captions.txt"));
fprintf("\nFigure 1 and Figure 4 revision 2 COMPLETED. Display only; frozen inputs; no new statistical analysis.\nOutput: %s\n", outDir);
end

%% ======================================================================
function sz = measureText(txt, fs, fn)
% Width and height (cm) of a text block at the given font size (one hidden scratch figure, reused).
    persistent f0 a0
    if isempty(f0) || ~isvalid(f0)
        f0 = figure("Visible","off","Units","centimeters","Position",[0 0 40 30]);
        a0 = axes(f0,"Units","normalized","Position",[0 0 1 1]); axis(a0,"off");
    end
    t0 = text(a0, 0.5, 0.5, txt, "Units","centimeters", "FontSize",fs, "FontName",fn, ...
        "HorizontalAlignment","center", "VerticalAlignment","middle", "Interpreter","none");
    e = t0.Extent; sz = e(3:4); delete(t0);
end

function out = wrapText(txt, maxW, fs, fn)
% Greedy word wrap of each line so that no line is wider than maxW (cm).
    paras = split(string(txt), newline); lines = strings(0,1);
    for p = 1:numel(paras)
        words = split(strtrim(paras(p))); cur = "";
        for w = 1:numel(words)
            cand = strtrim(cur + " " + words(w));
            if cur ~= "" && measureText(cand, fs, fn)*[1;0] > maxW
                lines(end+1,1) = cur; cur = words(w); %#ok<AGROW>
            else
                cur = cand;
            end
        end
        lines(end+1,1) = cur; %#ok<AGROW>
    end
    out = strjoin(lines, newline);
end

function drawTB(ax, xc, y0, w, hh, txt, shade, fs, fn)
% Rounded box of width w and height hh with bottom edge at y0, centred at xc (data units = cm).
    rectangle(ax, "Position",[xc-w/2, y0, w, hh], "Curvature",[0.06 0.18], ...
        "FaceColor",shade*[1 1 1], "EdgeColor",[0.3 0.3 0.3], "LineWidth",0.6);
    text(ax, xc, y0+hh/2, txt, "HorizontalAlignment","center", "VerticalAlignment","middle", ...
        "FontSize",fs, "FontName",fn, "Interpreter","none");
end

function legendAxes(ax, labels, shades, pats, fs, fn)
% Stand-alone legend with grey swatches and hatch patterns, drawn in its own axes.
    axis(ax,"off"); hold(ax,"on"); xlim(ax,[0 1]); ylim(ax,[0 1]);
    n = numel(labels); sw = 0.11; sh = 0.075; gap = 0.13; y0 = 0.5 + n*gap/2;
    for k = 1:n
        yk = y0 - k*gap;
        rectangle(ax, "Position",[0.02, yk, sw, sh], "FaceColor",shades(k)*[1 1 1], "EdgeColor","k");
        if pats(k) ~= "none", hatchRect(ax, 0.02, 0.02+sw, yk, yk+sh, pats(k)); end
        text(ax, 0.02+sw+0.03, yk+sh/2, labels(k), "FontSize",fs, "FontName",fn, "VerticalAlignment","middle", "Interpreter","none");
    end
end

function hatchBars(ax, b, pats)
% Black hatching inside vertical stacked-bar segments.
    for q = 1:numel(b)
        if pats(q) == "none", continue; end
        xc = b(q).XEndPoints; yt = b(q).YEndPoints; hgt = b(q).YData; w = b(q).BarWidth * min(diff(b(q).XData));
        for k = 1:numel(xc)
            if ~isfinite(hgt(k)) || hgt(k) == 0, continue; end
            hatchRect(ax, xc(k)-w/2, xc(k)+w/2, yt(k)-hgt(k), yt(k), pats(q));
        end
    end
end

function hatchBarsH(ax, b, pats)
% Black hatching inside horizontal stacked-bar segments (barh).
    for q = 1:numel(b)
        if pats(q) == "none", continue; end
        yc = b(q).XEndPoints; xe = b(q).YEndPoints; len = b(q).YData; w = b(q).BarWidth * min(diff(b(q).XData));
        for k = 1:numel(yc)
            if ~isfinite(len(k)) || len(k) == 0, continue; end
            hatchRect(ax, xe(k)-len(k), xe(k), yc(k)-w/2, yc(k)+w/2, pats(q));
        end
    end
end

function hatchRect(ax, x0, x1, y0, y1, pat)
% Hatch pattern inside a data-space rectangle (normalized spacing; works with reversed axes).
    xl = xlim(ax); yl = ylim(ax); dx = diff(xl); dy = diff(yl);
    u0 = (min(x0,x1)-xl(1))/dx; u1 = (max(x0,x1)-xl(1))/dx; v0 = (min(y0,y1)-yl(1))/dy; v1 = (max(y0,y1)-yl(1))/dy;
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

function s = fmtInt(x)
    s = string(sprintf("%d", x));
    if x >= 1000, s = string(sprintf("%d,%03d", floor(x/1000), mod(x,1000))); end
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
