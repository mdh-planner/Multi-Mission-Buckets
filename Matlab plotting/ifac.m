% ===================== USER SETTINGS =====================

N = 20; % number of tasks
plot_route = 0; % 1 to plot, 0 no plot

approach = "OURS_" + N + "_0";

planID = 2;   % 1 = original plan, 2 = replan

if (planID == 1)
    close all; clearvars;
    planID = 1;
end

%%
path = "";  %******************* ENTER PATH TO THE FILES ************************

folder = path + approach + "\";
% =========================================================

%% --------------------- READ META ---------------------
metaPath = fullfile(folder,"sol_meta_" + planID + ".json");
meta = jsondecode(fileread(metaPath));
R     = meta.R;
M     = meta.M;
H     = meta.H;
Delta = meta.Delta;

% ---- FIXED MISSION COLOR MAP (consistent across all figures) ----
% We assume mission ids are 0..M-1. If your ids are sparse, adapt getColor().
Mtotal  = max([M, 8]);                 % stable size (>= #missions)
cmapFix = lines(Mtotal);               % fixed palette
getColor = @(m) cmapFix(m+1,:);        % mission id (0-based) -> RGB row

%% --------------------- ROUTES FIGURE ---------------------
if (plot_route)
    routesPath = fullfile(folder, "sol_routes_" + planID + ".csv");
    try
        T = readtable(routesPath,'VariableNamingRule','preserve');
    catch
        T = readtable(routesPath);
    end
    if isempty(T), error('sol_routes_%d.csv is empty or missing.', planID); end

    vn = T.Properties.VariableNames;
    robot       = pick(vn,'robot');
    seq         = pick(vn,'seq');
    fromX       = pick(vn,'from_x','fromX','fx');
    fromY       = pick(vn,'from_y','fromY','fy');
    toX         = pick(vn,'to_x','toX','tx');
    toY         = pick(vn,'to_y','toY','ty');
    fromTask    = pick(vn,'from_task','fromTask');
    toTask      = pick(vn,'to_task','toTask');
    fromMission = pick(vn,'from_mission','fromMission','mFrom');
    toMission   = pick(vn,'to_mission','toMission','mTo');

    figRoutes = figure('Color','w','Name','MMOP Routes','Units','pixels','Position',[100 100 1200 620]);
    clf; hold on; axis equal; grid on; xlabel('x'); ylabel('y'); title('Robot travel paths');

    isTaskFrom = T.(fromTask) >= 0;
    isTaskTo   = T.(toTask)   >= 0;

    tasksXY = unique([ ...
        T.(fromX)(isTaskFrom) T.(fromY)(isTaskFrom) T.(fromTask)(isTaskFrom) T.(fromMission)(isTaskFrom); ...
        T.(toX)(isTaskTo)     T.(toY)(isTaskTo)     T.(toTask)(isTaskTo)     T.(toMission)(isTaskTo)], ...
        'rows','stable');

    depXY = unique([ ...
        T.(fromX)(~isTaskFrom) T.(fromY)(~isTaskFrom) repmat(-1,sum(~isTaskFrom),1) repmat(-1,sum(~isTaskFrom),1); ...
        T.(toX)(~isTaskTo)     T.(toY)(~isTaskTo)     repmat(-1,sum(~isTaskTo),1)   repmat(-1,sum(~isTaskTo),1)], ...
        'rows','stable');

    % Mission-aware colors (use fixed map)
    if ~isempty(tasksXY)
        if istable(tasksXY), tasksXY = tasksXY{:,:}; end
        missionVals = round(tasksXY(:,4)); % mission ids from file
        ms = 50;
        for m0 = unique(missionVals,'stable').'
            idx = (missionVals == m0);
            scatter(tasksXY(idx,1), tasksXY(idx,2), ms, ...
                'MarkerFaceColor', getColor(m0), 'MarkerEdgeColor','k', ...
                'DisplayName', sprintf('M%d task', m0));
        end
        for k = 1:size(tasksXY,1)
            text(tasksXY(k,1), tasksXY(k,2), sprintf(' %d', tasksXY(k,3)), ...
                'FontSize',14, 'Color',[0.1 0.1 0.1], 'VerticalAlignment','bottom');
        end
    end

    if ~isempty(depXY)
        scatter(depXY(:,1), depXY(:,2), 80, 'w', 'filled', 'MarkerEdgeColor','k', ...
            'Marker','d', 'DisplayName','depot');
    end

    robots = unique(T.(robot),'stable');
    cmapR  = lines(max(numel(robots),3));
    LW = 1.6;

    for rix = 1:numel(robots)
        r  = robots(rix);
        RR = sortrows(T(T.(robot)==r, :), seq);
        for k = 1:height(RR)
            x1 = RR.(fromX)(k); y1 = RR.(fromY)(k);
            x2 = RR.(toX)(k);   y2 = RR.(toY)(k);
            plot([x1 x2],[y1 y2], '-', 'LineWidth', LW, ...
                'Color', cmapR(rix,:), 'HandleVisibility','off'); hold on;
        end
        plot(NaN,NaN,'-','LineWidth',LW,'Color',cmapR(rix,:), 'DisplayName', sprintf('r%d path', r));
    end
    legend('Location','eastoutside');
    hold off;
end
%% --------------------- READ SCHEDULE / BUCKETS ---------------------
schedPath   = fullfile(folder,"sol_schedule_" + planID + ".csv");
bucketsPath = fullfile(folder,"sol_buckets_" + planID + ".csv");

%% original schedule
% schedPath1   = fullfile(folder,"sol_schedule_1.csv");
% schedOriginal   = readtable(schedPath1,   'VariableNamingRule','preserve');
%%
try
    sched   = readtable(schedPath,   'VariableNamingRule','preserve');
    buckets = readtable(bucketsPath, 'VariableNamingRule','preserve');
catch
    sched   = readtable(schedPath);
    buckets = readtable(bucketsPath);
end
if isempty(sched)
    warning('sol_schedule_%d.csv is empty. Nothing to plot.', planID); return;
end

sched = sortrows(sched,"task");
% Which tasks are NEW in plan 2? (from C++ export)
newSet = [];
if planID == 2
    newPath = fullfile(folder, "sol_new_2.csv");
    if exist(newPath,'file')
        try
            NT = readtable(newPath,'VariableNamingRule','preserve');
        catch
            NT = readtable(newPath);
        end
        if ~isempty(NT)
            tcol = pickCol(NT.Properties.VariableNames,'task');
            newSet = unique(double(NT.(tcol)).');
        end
    else
        % fallback using mapping file (optional export from C++)
        mapPath = fullfile(folder, "sol_origmap_2.csv");
        if exist(mapPath,'file')
            try
                OM = readtable(mapPath,'VariableNamingRule','preserve');
            catch
                OM = readtable(mapPath);
            end
            jcol = pickCol(OM.Properties.VariableNames,'task');
            ocol = pickCol(OM.Properties.VariableNames,'orig_id','origin','old');
            newSet = unique(double(OM.(jcol)(OM.(ocol) < 0)).');
        end
    end
end
isNewTask = @(tid) ~isempty(newSet) && any(tid == newSet);

% Columns for the active plan
vnS = sched.Properties.VariableNames;
robotCol   = pickCol(vnS,'robot');
taskCol    = pickCol(vnS,'task');
missionCol = pickCol(vnS,'mission');
startCol   = pickCol(vnS,'start');
durCol     = pickCol(vnS,'duration','dur','service');
% sched = sortrows(sched, {robotCol, startCol});
sched = sortrows(sched,"task");

% Release map (for plan 2 absolute alignment)
relMap = containers.Map('KeyType','double','ValueType','double');
for r = 0:R-1, relMap(r) = 0.0; end

% Also read original plan for executed overlay & original task set
schedOriginal = [];
origTasks = [];
if planID == 2
    relPath = fullfile(folder, "sol_release_2.csv");
    if exist(relPath,'file')
        try
            RLS = readtable(relPath, 'VariableNamingRule','preserve');
        catch
            RLS = readtable(relPath);
        end
        rvn   = RLS.Properties.VariableNames;
        rcol2 = pickCol(rvn,'robot');
        relc  = pickCol(rvn,'release','rel','offset');
        for k = 1:height(RLS), relMap(RLS.(rcol2)(k)) = RLS.(relc)(k); end
    end

    schedPath1 = fullfile(folder,"sol_schedule_1.csv");
    if exist(schedPath1,'file')
        try
            schedOriginal = readtable(schedPath1,'VariableNamingRule','preserve');
        catch
            schedOriginal = readtable(schedPath1);
        end
        if ~isempty(schedOriginal)
            vnS1 = schedOriginal.Properties.VariableNames;
            robotCol1   = pickCol(vnS1,'robot');
            taskCol1    = pickCol(vnS1,'task');
            missionCol1 = pickCol(vnS1,'mission');
            startCol1   = pickCol(vnS1,'start');
            durCol1     = pickCol(vnS1,'duration','dur','service');
            %             schedOriginal = sortrows(schedOriginal, {robotCol1, startCol1});
            origTasks = unique(schedOriginal.(taskCol1),'stable');
        end
    end
end

%% ================== ID SHIFT (keep old numbering, append new) ==================
% % Build the set of fully executed task IDs (plan 1) up to each robot's release.
executed_ids = [];

if planID == 2 && ~isempty(schedOriginal)
    schedOriginal = sortrows(schedOriginal, "task");
    vnS1 = schedOriginal.Properties.VariableNames;
    robotCol1   = pickCol(vnS1,'robot');
    taskCol1    = pickCol(vnS1,'task');
    startCol1   = pickCol(vnS1,'start');
    durCol1     = pickCol(vnS1,'duration','dur','service');

    for i = 0:R-1
        rows = (schedOriginal.(robotCol1) == i);
        Ti1  = schedOriginal(rows,:);
        if isempty(Ti1), continue; end
        rel = 0; if isKey(relMap, i), rel = relMap(i); end

        maxRel = max(cell2mat(values(relMap)));
        % A task is considered "executed before replan" iff it finished before the release.
        finished = Ti1.(startCol1) + Ti1.(durCol1) <= maxRel + 1e-3;
        executed_ids = [executed_ids; double(Ti1.(taskCol1)(finished))]; %#ok<AGROW>
    end
    executed_ids = unique(executed_ids(:)');

    % If there were executed tasks, shift all plan-2 task ids by your rule:
    % new_id(k) = k + count(executed_ids <= k), for k >= 0 (do not touch depots = -1)
    if ~isempty(executed_ids)
        % --- Shift in sched (plan 2)
        v = double(sched.(taskCol));
        mask = v >= 0;
        v(mask) = shift_ids_vec(v(mask), executed_ids);
        sched.(taskCol) = v;

        % --- Shift in routes (from/to task columns) if present
        if exist('T','var') && ~isempty(T)
            vn = T.Properties.VariableNames;
            fromTask    = pick(vn,'from_task','fromTask');
            toTask      = pick(vn,'to_task','toTask');

            vf = double(T.(fromTask));  mf = vf >= 0;
            vt = double(T.(toTask));    mt = vt >= 0;

            vf(mf) = shift_ids_vec(vf(mf), executed_ids);
            vt(mt) = shift_ids_vec(vt(mt), executed_ids);

            T.(fromTask) = vf;
            T.(toTask)   = vt;
        end

        % --- Shift in "newSet" (so N/O markers line up with shifted labels)
        if ~isempty(newSet)
            newSet = sched.task(end-numel(newSet)+1:end);
        end

        % Optional: expose the mapping for debugging/comparison
        % old_ids = unique([v(mask); vf(mf); vt(mt)]);
        % idmap = table(old_ids, shift_ids_vec(old_ids, executed_ids), ...
        %     'VariableNames', {'old_id','new_id'});
        % disp(idmap)
    end

end




%% --------------------- GANTT FIGURE ---------------------
figG = figure('Color','w','Name','MMOP Gantt', ...
    'Units','pixels','Position',[60 80 1200 560]);
ax = axes('Parent',figG); hold(ax,'on'); box(ax,'on');

yticks(1:R);
yticklabels(arrayfun(@(r) sprintf('r%d', r), 0:R-1, 'uni',0));
ylim([0.5, R+0.5]);
xlim([0, H*Delta]);
ylabel('robot');

% bottom axis: minutes (edges) + bucket IDs (centers)
edges   = (0:H) * Delta;
centers = edges(1:end-1) + Delta/2;
buckIDs = arrayfun(@(h) sprintf('b%02d', h), 0:H-1, 'uni', 0);

ticksAll = [edges, centers];
[xticksSorted, order] = sort(ticksAll);
edgeLbls   = string(edges);
centerLbls = string(buckIDs);
labelsAll  = [edgeLbls, centerLbls];
set(ax, 'XTick', xticksSorted, 'XTickLabel', labelsAll(order), ...
    'TickDir','out', 'Layer','top');
xlabel('time (min)  |  buckets');

% alternating shading + vertical grid
for h = 0:H-1
    if mod(h,2)==1
        patch('Parent',ax, ...
            'XData',[h*Delta (h+1)*Delta (h+1)*Delta h*Delta], ...
            'YData',[0.5 0.5 R+0.5 R+0.5], ...
            'FaceColor',[0.95 0.95 0.95], 'EdgeColor','none', ...
            'HandleVisibility','off', 'HitTest','off');
    end
end
for h = 0:H
    xline(h*Delta, ':', 'Color', [0.7 0.7 0.7], 'HandleVisibility','off');
end

% legend holders for missions 0..M-1
legendHandlesAdded = false(1,M);
legH = gobjects(1,M);
missionNames = arrayfun(@(m) sprintf('M%d', m), 0:M-1, 'uni', 0);

% bar height
hBar = 0.6;

%% --- Overlay executed portion of plan 1 up to per-robot release (semi-transparent) ---
if planID == 2 && ~isempty(schedOriginal)
    vnS1 = schedOriginal.Properties.VariableNames;
    robotCol1   = pickCol(vnS1,'robot');
    taskCol1    = pickCol(vnS1,'task');
    missionCol1 = pickCol(vnS1,'mission');
    startCol1   = pickCol(vnS1,'start');
    durCol1     = pickCol(vnS1,'duration','dur','service');

    alphaExec = 0.1;
    for i = 0:R-1
        rows = (schedOriginal.(robotCol1) == i);
        Ti1  = schedOriginal(rows,:);
        if isempty(Ti1), continue; end
        rel = 0; if isKey(relMap, i), rel = relMap(i); end

        for k = 1:height(Ti1)
            m0 = Ti1.(missionCol1)(k);
            x0 = Ti1.(startCol1)(k);
            w0 = Ti1.(durCol1)(k);
            x1 = x0;
            x2 = min(x0 + w0, rel);     % executed portion only
            if x2 <= x1, continue; end

            y  = i + 0.92 - hBar/2;
            yT = y + hBar;

            patch('Parent', ax, ...
                'XData', [x1 x2 x2 x1], ...
                'YData', [y  y  yT yT], ...
                'FaceColor', getColor(m0), ...
                'FaceAlpha', alphaExec, ...
                'EdgeColor', 'k', ...
                'LineStyle', '-', ...
                'LineWidth', 0.1, ...
                'HandleVisibility','off', ...
                'HitTest','off');

            % faint "O" tag at end of executed chunk (optional)
            text(x2, yT - 0.02, 'O', ...
                'FontSize',6, 'HorizontalAlignment','right', ...
                'VerticalAlignment','top', 'Color',[0 0 0]);
            % task id label
            text(x1 + 0.02*w0, y + 0.5*hBar, sprintf('%d', Ti1.(taskCol)(k)), ...
                'FontSize',8, 'HorizontalAlignment','left', 'VerticalAlignment','middle','Color', [0 0 0]);

        end
    end
end

%% -------- read reservations and build resMap (before drawing bars) --------
yPath = fullfile(folder, "sol_y_" + planID + ".csv");
zPath = fullfile(folder, "sol_z_" + planID + ".csv");

resMap = NaN(R, H);      % default: no reservation
yz = [];                 % keep table for later overlay

if exist(yPath,'file') || exist(zPath,'file')
    readY = exist(yPath,'file') ~= 0;
    try
        yz = readtable(readY*yPath + ~readY*zPath, 'VariableNamingRule','preserve');
    catch
        if readY, yz = readtable(yPath); else, yz = readtable(zPath); end
    end

    vny  = yz.Properties.VariableNames;
    rcol = pickCol(vny,'robot');
    hcol = pickCol(vny,'bucket','h');
    mcol = pickCol(vny,'mission','m');
    vcol = pickCol(vny,'value','val');

    % one mission per (robot,bucket) in OURS (F1); if multiple, keep first
    for i = 0:R-1
        for h = 0:H-1
            rows = yz.(rcol)==i & yz.(hcol)==h & yz.(vcol)>0.5;
            if any(rows)
                mids = unique(yz.(mcol)(rows),'stable');
                resMap(i+1, h+1) = mids(1);
            end
        end
    end
end

%% -------- Build resMap from schedule(s) only (no b / y / z needed) --------
% resMap(i+1,h+1) = mission id reserved/used in bucket h for robot i,
% inferred from the mission of the task that *starts* in that bucket.

resMap = NaN(R, H);  % default = no info

% Helper: majority mission in a vector (mode that works for doubles)
mode_or_first = @(v) v(1) * (numel(v)==1) + (numel(v)>1) * mode(v);

% ---------- 1) Fill pre-cut buckets from original plan's executed starts ----------
if planID == 2 && ~isempty(schedOriginal)
    vnS1 = schedOriginal.Properties.VariableNames;
    robotCol1   = pickCol(vnS1,'robot');
    missionCol1 = pickCol(vnS1,'mission');
    startCol1   = pickCol(vnS1,'start');

    for i = 0:R-1
        rel_i = 0; if isKey(relMap, i), rel_i = relMap(i); end
        rows  = (schedOriginal.(robotCol1) == i);

        if ~any(rows), continue; end
        S1i = schedOriginal(rows,:);

        % Consider only tasks that START before the robot's release (pre-cut)
        starts  = double(S1i.(startCol1));
        maskPre = starts < rel_i - 1e-9;
        if ~any(maskPre), continue; end

        hStart = floor( starts(maskPre) / Delta );                    % absolute buckets in plan 1
        mids   = double(S1i.(missionCol1)(maskPre));

        % For each bucket, take the majority mission among starts in that bucket
        for h = 0:H-1
            idx = (hStart == h);
            if any(idx)
                resMap(i+1, h+1) = mode_or_first(mids(idx));
            end
        end
    end
end

% ---------- 2) Fill/overwrite post-cut buckets from active plan (planID 1 or 2) ----------
vnS = sched.Properties.VariableNames;
robotCol   = pickCol(vnS,'robot');
missionCol = pickCol(vnS,'mission');
startCol   = pickCol(vnS,'start');

for i = 0:R-1
    rows = (sched.(robotCol) == i);
    if ~any(rows), continue; end
    Si = sched(rows,:);

    % In plan 2, your bars are plotted at x = start + relMap(i)
    rel_i = 0; if isKey(relMap, i), rel_i = relMap(i); end
    startsAbs = double(Si.(startCol)) + rel_i;   % absolute start time on the same global axis

    hStart = floor( startsAbs / Delta );
    mids   = double(Si.(missionCol));

    for h = 0:H-1
        idx = (hStart == h);
        if any(idx)
            resMap(i+1, h+1) = mode_or_first(mids(idx));
        end
    end
end


%% -------- Draw main plan bars (plan 1 or plan 2) --------
for i = 0:R-1
    rows = (sched.(robotCol) == i);
    Ti   = sched(rows,:);
    for k = 1:height(Ti)
        m0 = Ti.(missionCol)(k);
        x  = Ti.(startCol)(k) + relMap(i);  % offset by release if plan 2
        w  = Ti.(durCol)(k);
        y  = i + 0.92 - hBar/2;

        % main bar (opaque)
%         rectangle('Parent',ax, 'Position',[x, y, w, hBar], ...
%             'FaceColor', getColor(m0), 'EdgeColor','k', ...
%             'HandleVisibility','off');
        rectPatch(ax, x, y, w, hBar, ...
          getColor(m0), 'k', 0.5);


        % task id label
        text(x + 0.02*w, y + 0.5*hBar, sprintf('%d', Ti.(taskCol)(k)), ...
            'FontSize',8, 'HorizontalAlignment','left', 'VerticalAlignment','middle','Color', [1 1 1]);

        % O/N marker (top-right)
        labelChar = 'O';
        if planID == 2
            if ismember(Ti.(taskCol)(k), newSet)
                labelChar = 'N';
            end
        end
        text(x + w - 0.02*w, y + hBar - 0.02, labelChar, ...
            'FontSize',6, 'HorizontalAlignment','right', 'VerticalAlignment','top', 'Color', [1 1 1], 'FontWeight', 'bold');

        % legend entry once per mission color
        if m0 >= 0 && m0 < M && ~legendHandlesAdded(m0+1)
            legH(m0+1) = plot(ax, NaN,NaN,'s','MarkerFaceColor',getColor(m0), ...
                'MarkerEdgeColor','k','DisplayName', missionNames{m0+1});
            legendHandlesAdded(m0+1) = true;
        end

        % Determine the start bucket of this task (reservation is for start bucket)
        hStartAbs = floor( (x + 1e-9) / Delta );   % x already includes relMap(i)

        if hStartAbs >= 0 && hStartAbs < H
            mRes = resMap(i+1, hStartAbs+1);   % reserved mission for (i, hStartAbs)
        else
            mRes = NaN;
        end

        % If reservation exists and doesn't match task mission, mark it
        if ~isnan(mRes) && mRes ~= m0
            % red outline to flag mismatch
%             rectangle('Parent',ax, 'Position',[x, y, w, hBar], ...
%                 'EdgeColor',[0.85 0.2 0.2], 'LineWidth', 2.0, ...
%                 'FaceColor','none', 'HandleVisibility','off');
            rectPatch(ax, x, y, w, hBar, ...
          'none', [0.85 0.2 0.2], 2.0);


            % add a small "!" tag
            text(x + w/2, y + hBar + 0.02, '!', 'Color',[0.85 0.2 0.2], ...
                'FontWeight','bold', 'HorizontalAlignment','center', ...
                'VerticalAlignment','bottom');
        end


    end
end

%% -------- Build bands from resMap (mission-inferred + gray if empty) --------
bandH  = 0.16;  
gapTop = 0.04;

for i = 0:R-1
    yBase = (i + 1) + hBar/2 + gapTop;
    for h = 0:H-1
        x1 = h * Delta;  
        x2 = (h + 1) * Delta;
        y1 = yBase;      
        y2 = y1 + bandH;

        m0 = resMap(i+1, h+1);

        if isnan(m0)
            % --- No reservation: draw gray rectangle ---
            patch('Parent', ax, ...
                'XData', [x1 x2 x2 x1], ...
                'YData', [y1 y1 y2 y2], ...
                'FaceColor', [1 1 1], ...
                'FaceAlpha', 1, ...
               'EdgeColor', [0.9 0.9 0.9], ...
                'LineWidth', 0.1, ...
                'HandleVisibility', 'off', ...
                'HitTest', 'off');
        else
             % --- No reservation: draw white rectangle ---
            patch('Parent', ax, ...
                'XData', [x1 x2 x2 x1], ...
                'YData', [y1 y1 y2 y2], ...
                'FaceColor', [1 1 1], ...
                'FaceAlpha', 1, ...
                'EdgeColor', 'none', ...  %'EdgeColor', [0 0 0], ...
                'LineWidth', 2, ...
                'HandleVisibility', 'off', ...
                'HitTest', 'off');
            % --- Reservation exists: draw mission-colored band ---
            patch('Parent', ax, ...
                'XData', [x1 x2 x2 x1], ...
                'YData', [y1 y1 y2 y2], ...
                'FaceColor', getColor(m0), ...
                'FaceAlpha', 0.25, ...
                'EdgeColor', getColor(m0), ...
                'LineWidth', 0.1, ...
                'HandleVisibility', 'off', ...
                'HitTest', 'off');
        end
    end
end

% Add legend sample
resSample = patch('Parent',ax, 'XData',nan, 'YData',nan, ...
    'FaceColor',[0.5 0.5 0.5], 'FaceAlpha',0.25, ...
    'EdgeColor','k', 'DisplayName','reservation (by mission)', ...
    'HandleVisibility','off');


%% -------- overlay reservations (y or z) as thin bands --------
% resSample = gobjects(1);
% if exist(yPath,'file') || exist(zPath,'file')
% 
%     rvals = unique(sched.(robotCol)).';
%     bandH  = 0.16;  gapTop = 0.04;
% 
%     for i = rvals
%         yBase = (i + 1) + hBar/2 + gapTop;
%         for h = 0:H-1
%             rows = yz.(rcol)==i & yz.(hcol)==h & yz.(vcol)>0.5;
%             if ~any(rows), continue; end
%             mids = unique(yz.(mcol)(rows),'stable');
%             for s = 1:numel(mids)
%                 m0 = mids(s);
%                 x1 = h * Delta;  x2 = (h + 1) * Delta;
%                 y1 = yBase + (s-1)*(bandH+0.02); y2 = y1 + bandH;
% 
%                 patch('Parent', ax, ...
%                     'XData', [x1 x2 x2 x1], ...
%                     'YData', [y1 y1 y2 y2], ...
%                     'FaceColor', getColor(m0), ...
%                     'FaceAlpha', 0.25, ...
%                     'EdgeColor', 'k', ...
%                     'LineWidth', 0.5, ...
%                     'HandleVisibility', 'off', ...
%                     'HitTest', 'off');
%             end
%         end
%     end
% 
%     resSample = patch('Parent',ax, 'XData',nan, 'YData',nan, ...
%         'FaceColor',[0.5 0.5 0.5], 'FaceAlpha',0.25, ...
%         'EdgeColor','k', 'DisplayName','reservation (by mission)', ...
%         'HandleVisibility','off');
% end
% 
% % legend (missions + optional reservation sample)
% legList  = legH(legendHandlesAdded);
% legNames = missionNames(legendHandlesAdded);
% legList  = legList(:).';   legNames = legNames(:).';
% 
% if ~isempty(resSample) && isgraphics(resSample)
%     legendHandles = [legList, resSample];
%     legendLabels  = [legNames, {'reservation (by mission)'}];
% else
%     legendHandles = legList;
%     legendLabels  = legNames;
% end
% if ~isempty(legendHandles), legend(legendHandles, legendLabels, 'Location','eastoutside'); end
% 
% title
if isfield(meta,'baseline') && isfield(meta,'objective')
    if (planID == 1)
        title(sprintf('Original Plan - Cost: %.3f',  meta.objective));
    else
        title(sprintf('Replan - Cost: %.3f',  meta.objective));
    end
end
hold(ax,'off');

%% ---------- bucket legend (cap vs quota) ----------
% Compute caps directly from the plotted reservations (resMap),
% so the table matches the figure.
caps = zeros(M, H);               % robots reserved per (mission, bucket)
for i = 1:R
    for h = 1:H
        m0 = resMap(i, h);
        if ~isnan(m0)
            caps(m0+1, h) = caps(m0+1, h) + 1;
        end
    end
end

% Get quotas (req) from CSV if present; otherwise assume zeros
quotas = zeros(M, H);
if exist('buckets','var') && ~isempty(buckets)
    bvn       = buckets.Properties.VariableNames;
    bucketCol = pickCol(bvn,'bucket','h');
    missionB  = pickCol(bvn,'mission','m');
    quotaCol  = pickCol(bvn,'quota','req','need');
    for h = 0:H-1
        for m0 = 0:M-1
            idx = buckets.(bucketCol)==h & buckets.(missionB)==m0;
            if any(idx)
                quotas(m0+1, h+1) = buckets.(quotaCol)(find(idx,1,'first'));
            end
        end
    end
end

% % Build the text box using caps (from resMap) and quotas (from CSV)
% rowsTxt = strings(H,1);
% for h = 0:H-1
%     pieces = strings(1,M);
%     for m0 = 0:M-1
%         got = caps(m0+1, h+1);
%         req = quotas(m0+1, h+1);
%         pieces(m0+1) = sprintf('M%d %d/%d', m0, got, req);
%     end
%     rowsTxt(h+1) = sprintf('b%02d  %s', h, strjoin(pieces,'   '));
% end

% axpos = get(ax,'Position');  % normalized
% x0 = min(axpos(1)+axpos(3)+0.12, 0.76);
% w  = 0.2;
% hBox = min(0.045*H + 0.04, 0.72);
% y0 = max(min(axpos(2)+0.55, 0.96 - hBox), 0.10);
% 
% annotation('textbox',[x0 y0 w hBox], ...
%     'String', sprintf('%s\n', rowsTxt), ...
%     'Interpreter','none', ...
%     'FitBoxToText','off', ...
%     'BackgroundColor',[1 1 1], ...
%     'EdgeColor',[0.7 0.7 0.7], ...
%     'FontName','Consolas', ...
%     'FontSize', 10, ...
%     'HorizontalAlignment','left', ...
%     'VerticalAlignment','top');

%% ================= helpers =================

    function name = pickCol(varNames, preferred, varargin)

        all = [{preferred} varargin];

        for t = 1:numel(all)
            if any(strcmp(varNames, all{t})), name = all{t}; return; end
        end

        idx = find(strcmpi(varNames, preferred),1,'first');

        if ~isempty(idx), name = varNames{idx}; return; end

        idx = find(contains(varNames, preferred,'IgnoreCase',true),1,'first');

        if ~isempty(idx), name = varNames{idx}; return; end

        error('Required column not found: %s', preferred);

    end

    function name = pick(varNames, preferred, varargin)

        all = [{preferred} varargin];

        for t = 1:numel(all)
            if any(strcmp(varNames, all{t})), name = all{t}; return; end
        end

        idx = find(strcmpi(varNames, preferred),1);

        if ~isempty(idx), name = varNames{idx}; return; end

        idx = find(contains(varNames, preferred,'IgnoreCase',true),1);

        if ~isempty(idx), name = varNames{idx}; return; end

        error('Required column not found: %s', preferred);

    end

    function new_ids = shift_ids_vec(old_ids, executed_ids)

        for i = 1 : numel(executed_ids)
            for j = 1 : numel(old_ids)

                if (executed_ids(i) == old_ids(j))
                    old_ids(j:end) = old_ids(j:end)+1;
                end

            end
        end
        new_ids = old_ids;
        return;
    end

    function h = rectPatch(ax, x, y, w, hBar, fc, ec, lw)
        % Draw a filled axis-aligned rectangle as a patch
        X = [x      x+w    x+w   x   ];
        Y = [y      y      y+hBar y+hBar];
        h = patch('Parent', ax, ...
                  'XData', X, 'YData', Y, ...
                  'FaceColor', fc, ...
                  'EdgeColor', ec, ...
                  'LineWidth', lw, ...
                  'HandleVisibility', 'off', ...
                  'HitTest', 'off');
    end



