%% ieee13_pq_flex_scan_lindisflow_yalmip_gurobi_v1_1_adaptive_scan.m
% IEEE-13 simplified balanced-feeder adaptive single-period P-Q flexibility scan.
%
% This version uses the IEEE 13-node distribution feeder topology as a
% simplified single-phase balanced radial feeder and solves a LinDistFlow
% optimization model with MATLAB + YALMIP + Gurobi. It is intended only to
% verify the first-stage per-period (P0,Q0) flexibility-domain scanning
% workflow. It does not model the original IEEE-13 three-phase unbalanced
% details such as voltage regulators, capacitors, phase-specific line codes, or
% phase-specific loads.
%
% Important simplification:
%   This version uses IEEE 13-node topology and simplified single-phase
%   equivalent line parameters. It is for validating the P-Q flexibility scan
%   process, not as a strict three-phase power-flow benchmark.
%
% Sign convention used in every CSV/figure:
%   The root feeder branch is 650 -> 632. P_root = P_650,632 and Q_root =
%   Q_650,632 are positive when power is imported from the upstream/root grid
%   into the distribution feeder. Final reported interface powers are
%       P0 = -P_root,  Q0 = -Q_root.
%   Therefore P0 > 0 means active-power export from the distribution feeder to
%   the upstream grid, and Q0 > 0 means reactive-power export to the upstream
%   grid.
%
% Adaptive scan note:
%   This version uses an adaptive support-function scan. It first applies a
%   small set of coarse directions to scan the P-Q flexibility boundary, then
%   inserts new directions based on normalized distances between adjacent
%   boundary points. This reduces repeated solves caused by fixed-direction
%   scans hitting the same boundary point and improves boundary-approximation
%   efficiency. It is still a first-stage per-period P-Q section scan and does
%   not represent the complete cross-time flexibility domain.
%
% Scope exclusions for this first version:
%   No external power-flow package dependency, no sensitivity-based AC validation,
%   no multi-period coupling, no generator ramping, no ESS SOC
%   recursion, no VG/VB decomposition, no bound shrinking, and no infeasible
%   point search. T=4 is used only to scan independent load/PV operating points.

clear; clc; close all;

%% --------------------------- User parameters ---------------------------
RESULT_DIR = 'results_ieee13_pq_flex_scan_lindisflow_v1_1_adaptive_scan';
SCAN_MODE = 'adaptive_scan';            % 'adaptive_scan' or 'fixed_angle'
FIXED_N_THETA = 36;
ADAPTIVE_INITIAL_N_THETA = 8;
ADAPTIVE_MAX_ROUNDS = 6;
ADAPTIVE_MAX_SOLVES = 80;
ADAPTIVE_MIN_NEW_ANGLES = 1;
ADAPTIVE_DIST_TOL = 0.08;
ADAPTIVE_Q_TOL = 0.08;
ADAPTIVE_P_TOL = 0.08;
DUP_POINT_TOL = 1e-5;
ANGLE_TOL = 1e-10;
T = 4;
LOAD_SCALE = [0.90, 1.00, 1.10, 0.95];
PV_AVAIL_PROFILE = [0.20, 0.65, 0.85, 0.40]; % MW available at node 675

baseMVA = 5;
V0 = 1.0;
VMIN = 0.95^2;                          % squared voltage lower bound
VMAX = 1.05^2;                          % squared voltage upper bound
USE_BRANCH_SOC = false;                 % false: linear P/Q box; true: optional SOCP branch limit
USE_PV_SOC = false;                     % false: PV inverter box; true: optional SOCP inverter limit
USE_ESS_SOC_CONE = false;               % false: ESS inverter box; true: optional SOCP inverter limit
USE_ESS_BINARY = false;                 % false: continuous relaxation allows simultaneous charge/discharge
YALMIP_VERBOSE = 0;
GUROBI_VERBOSE = 0;

%% ------------------------- Initialization ------------------------------
if exist(RESULT_DIR, 'dir')
    warning('Result directory already exists. New run outputs may overwrite generated files inside %s.', RESULT_DIR);
else
    mkdir(RESULT_DIR);
end
log_file = fullfile(RESULT_DIR, 'run_log.txt');
if exist(log_file, 'file')
    delete(log_file);
end
diary(log_file);
diary on;
cleanup_diary = onCleanup(@() diary('off')); %#ok<NASGU>

fprintf('IEEE-13 simplified balanced LinDistFlow P-Q scan started at %s\n', datestr(now)); drawnow;
fprintf('Result directory: %s\n', RESULT_DIR); drawnow;
fprintf('SCAN_MODE=%s, T=%d, FIXED_N_THETA=%d, ADAPTIVE_INITIAL_N_THETA=%d, baseMVA=%.3f\n', ...
    SCAN_MODE, T, FIXED_N_THETA, ADAPTIVE_INITIAL_N_THETA, baseMVA); drawnow;

check_dependencies();
sys = build_ieee13_simplified_system(baseMVA, V0, VMIN, VMAX);
ts = build_time_series_data(sys, T, LOAD_SCALE, PV_AVAIL_PROFILE);

params = struct();
params.RESULT_DIR = RESULT_DIR;
params.SCAN_MODE = SCAN_MODE;
params.FIXED_N_THETA = FIXED_N_THETA;
params.ADAPTIVE_INITIAL_N_THETA = ADAPTIVE_INITIAL_N_THETA;
params.ADAPTIVE_MAX_ROUNDS = ADAPTIVE_MAX_ROUNDS;
params.ADAPTIVE_MAX_SOLVES = ADAPTIVE_MAX_SOLVES;
params.ADAPTIVE_MIN_NEW_ANGLES = ADAPTIVE_MIN_NEW_ANGLES;
params.ADAPTIVE_DIST_TOL = ADAPTIVE_DIST_TOL;
params.ADAPTIVE_Q_TOL = ADAPTIVE_Q_TOL;
params.ADAPTIVE_P_TOL = ADAPTIVE_P_TOL;
params.DUP_POINT_TOL = DUP_POINT_TOL;
params.ANGLE_TOL = ANGLE_TOL;
params.T = T;
params.baseMVA = baseMVA;
params.V0 = V0;
params.VMIN = VMIN;
params.VMAX = VMAX;
params.USE_BRANCH_SOC = USE_BRANCH_SOC;
params.USE_PV_SOC = USE_PV_SOC;
params.USE_ESS_SOC_CONE = USE_ESS_SOC_CONE;
params.USE_ESS_BINARY = USE_ESS_BINARY;
params.YALMIP_VERBOSE = YALMIP_VERBOSE;
params.GUROBI_VERBOSE = GUROBI_VERBOSE;

all_rows = struct([]);
dispatch_solutions = cell(T, 1);
summary_rows = struct([]);

%% ------------------------- Main scan loop ------------------------------
for t = 1:T
    fprintf('\n--- IEEE-13 LinDistFlow single-period scan t=%d/%d, mode=%s ---\n', t, T, SCAN_MODE); drawnow;
    if strcmpi(SCAN_MODE, 'adaptive_scan')
        [period_rows, period_dispatch, period_summary] = run_adaptive_scan_for_period(t, sys, ts, params);
    elseif strcmpi(SCAN_MODE, 'fixed_angle')
        [period_rows, period_dispatch, period_summary] = run_fixed_scan_for_period(t, sys, ts, params);
    else
        error('Unknown SCAN_MODE: %s', SCAN_MODE);
    end
    all_rows = append_struct(all_rows, period_rows);
    summary_rows = append_struct(summary_rows, period_summary);
    dispatch_solutions{t} = period_dispatch;
    plot_pq_flex_region(period_rows, t, params, period_summary);
end

%% ------------------------- Save outputs --------------------------------
save_results_csv(all_rows, summary_rows, RESULT_DIR);
save(fullfile(RESULT_DIR, 'dispatch_solutions.mat'), 'dispatch_solutions', 'all_rows', 'summary_rows', 'sys', 'ts', 'params', '-v7.3');
fprintf('Saved MAT: %s\n', fullfile(RESULT_DIR, 'dispatch_solutions.mat')); drawnow;
fprintf('\nCompleted IEEE-13 LinDistFlow P-Q scan at %s\n', datestr(now)); drawnow;
fprintf('Outputs saved under %s\n', RESULT_DIR); drawnow;

%% ============================= Local functions =========================
function check_dependencies()
    missing = {};
    if exist('sdpvar', 'file') ~= 2 || exist('optimize', 'file') ~= 2 || exist('sdpsettings', 'file') ~= 2
        missing{end+1} = 'YALMIP (install YALMIP and add it to the MATLAB path)'; %#ok<AGROW>
    end
    if exist('gurobi', 'file') ~= 3 && exist('gurobi', 'file') ~= 2 && exist('gurobi_mex', 'file') ~= 3 && exist('gurobi_mex', 'file') ~= 2
        missing{end+1} = 'Gurobi MATLAB interface (install Gurobi, activate a license, and add its MATLAB folder to the path)'; %#ok<AGROW>
    end
    if ~isempty(missing)
        fprintf('Dependency check failed. Missing components:\n'); drawnow;
        for i = 1:numel(missing)
            fprintf('  - %s\n', missing{i}); drawnow;
        end
        error('Please install/add the missing dependencies before running this script. MATPOWER is not required for this IEEE-13 LinDistFlow version.');
    end
    fprintf('Dependency check passed: YALMIP and Gurobi are visible. MATPOWER is not required.\n'); drawnow;
end

function sys = build_ieee13_simplified_system(baseMVA, V0, VMIN, VMAX)
    % Bus mapping used by this single-phase balanced equivalent:
    %   id  name
    %    1  650  root/PCC/slack
    %    2  632
    %    3  633
    %    4  634
    %    5  645
    %    6  646
    %    7  671
    %    8  680
    %    9  684
    %   10  611
    %   11  652
    %   12  692
    %   13  675
    bus_name = {'650','632','633','634','645','646','671','680','684','611','652','692','675'};
    bus_id = (1:numel(bus_name)).';
    branch_name = {'650-632','632-633','633-634','632-645','645-646','632-671', ...
                   '671-680','671-684','684-611','684-652','671-692','692-675'};
    from = [1; 2; 3; 2; 5; 2; 7; 7; 9; 9; 7; 12];
    to   = [2; 3; 4; 5; 6; 7; 8; 9;10;11;12; 13];

    nb = numel(bus_name);
    nl = numel(from);
    if nl ~= nb - 1
        error('IEEE-13 simplified topology is not radial by count: Nbranch=%d, Nbus-1=%d.', nl, nb-1);
    end
    if ~is_connected_radial(nb, from, to)
        error('IEEE-13 simplified topology failed radial connectivity check.');
    end

    % Simplified positive-sequence-equivalent line parameters in per-unit on
    % baseMVA. They are deliberately approximate and chosen only to exercise
    % LinDistFlow constraints on the IEEE-13 radial topology.
    r = [0.010; 0.006; 0.006; 0.008; 0.006; 0.012; 0.004; 0.007; 0.005; 0.005; 0.006; 0.007];
    x = [0.030; 0.018; 0.018; 0.024; 0.018; 0.036; 0.012; 0.021; 0.015; 0.015; 0.018; 0.021];
    Smax = [4.0; 2.5; 1.8; 1.6; 1.3; 3.5; 1.0; 1.8; 0.8; 0.8; 1.8; 1.6]; % MW/MVAr box or MVA SOC radius

    sys = struct();
    sys.baseMVA = baseMVA;
    sys.V0 = V0;
    sys.VMIN = VMIN;
    sys.VMAX = VMAX;
    sys.bus_name = bus_name;
    sys.bus_id = bus_id;
    sys.nb = nb;
    sys.nl = nl;
    sys.root_bus = 1;
    sys.root_branch = 1;
    sys.branch_name = branch_name;
    sys.from = from;
    sys.to = to;
    sys.r = r;
    sys.x = x;
    sys.Smax = Smax;
    sys.children_branches = cell(nb, 1);
    for b = 1:nb
        sys.children_branches{b} = find(from == b);
    end

    sys.gen_bus = find_bus(sys, '671');
    sys.pv_bus = find_bus(sys, '675');
    sys.ess_bus = find_bus(sys, '680');
    sys.PG_min = 0.00;
    sys.PG_max = 1.20;
    sys.QG_min = -0.60;
    sys.QG_max = 0.60;
    sys.SPV = 1.00;
    sys.QPV_max = 0.50;
    sys.PCH_max = 0.50;
    sys.PDIS_max = 0.50;
    sys.QESS_max = 0.35;
    sys.SESS = 0.60;
    % ESS energy/SOC parameters are intentionally omitted in this single-period
    % version; no E_t variable or SOC recursion is used.

    fprintf('Built IEEE-13 simplified balanced feeder with %d buses and %d radial branches.\n', nb, nl); drawnow;
    fprintf('Root/PCC bus=%s, root branch=%s, G1 bus=%s, PV bus=%s, ESS bus=%s.\n', ...
        sys.bus_name{sys.root_bus}, sys.branch_name{sys.root_branch}, sys.bus_name{sys.gen_bus}, sys.bus_name{sys.pv_bus}, sys.bus_name{sys.ess_bus}); drawnow;
end

function ok = is_connected_radial(nb, from, to)
    adjacency = cell(nb, 1);
    for e = 1:numel(from)
        adjacency{from(e)} = [adjacency{from(e)}, to(e)]; %#ok<AGROW>
        adjacency{to(e)} = [adjacency{to(e)}, from(e)]; %#ok<AGROW>
    end
    visited = false(nb, 1);
    queue = 1;
    visited(1) = true;
    while ~isempty(queue)
        b = queue(1);
        queue(1) = [];
        nbrs = adjacency{b};
        for k = 1:numel(nbrs)
            if ~visited(nbrs(k))
                visited(nbrs(k)) = true;
                queue(end+1) = nbrs(k); %#ok<AGROW>
            end
        end
    end
    ok = all(visited) && numel(from) == nb - 1;
end

function idx = find_bus(sys, name)
    idx = find(strcmp(sys.bus_name, name), 1);
    if isempty(idx)
        error('Bus name %s not found.', name);
    end
end

function ts = build_time_series_data(sys, T, load_scale, pv_avail_profile)
    if numel(load_scale) ~= T
        error('LOAD_SCALE length (%d) must equal T (%d).', numel(load_scale), T);
    end
    if numel(pv_avail_profile) ~= T
        error('PV_AVAIL_PROFILE length (%d) must equal T (%d).', numel(pv_avail_profile), T);
    end
    Pload_base = zeros(sys.nb, 1);
    Qload_base = zeros(sys.nb, 1);
    Pload_base(find_bus(sys, '634')) = 0.40; Qload_base(find_bus(sys, '634')) = 0.18;
    Pload_base(find_bus(sys, '646')) = 0.25; Qload_base(find_bus(sys, '646')) = 0.10;
    Pload_base(find_bus(sys, '652')) = 0.18; Qload_base(find_bus(sys, '652')) = 0.08;
    Pload_base(find_bus(sys, '671')) = 1.15; Qload_base(find_bus(sys, '671')) = 0.55;
    Pload_base(find_bus(sys, '675')) = 0.85; Qload_base(find_bus(sys, '675')) = 0.40;
    Pload_base(find_bus(sys, '692')) = 0.17; Qload_base(find_bus(sys, '692')) = 0.08;

    ts = struct();
    ts.T = T;
    ts.LOAD_SCALE = load_scale(:);
    ts.PV_AVAIL_PROFILE = pv_avail_profile(:);
    ts.Pload_base = Pload_base;
    ts.Qload_base = Qload_base;
    ts.Pload = zeros(sys.nb, T);
    ts.Qload = zeros(sys.nb, T);
    for t = 1:T
        ts.Pload(:, t) = Pload_base * load_scale(t);
        ts.Qload(:, t) = Qload_base * load_scale(t);
        fprintf('Time data t=%d: load scale=%.3f, total Pload=%.3f MW, total Qload=%.3f MVAr, PV_avail=%.3f MW.\n', ...
            t, load_scale(t), sum(ts.Pload(:, t)), sum(ts.Qload(:, t)), pv_avail_profile(t)); drawnow;
    end
end

function theta_list = get_fixed_angles(Ntheta)
    theta_list = linspace(0, 2*pi, Ntheta + 1);
    theta_list(end) = [];
end

function [period_rows, period_dispatch, period_summary] = run_fixed_scan_for_period(t, sys, ts, params)
    theta_list = get_fixed_angles(params.FIXED_N_THETA);
    period_rows = struct([]);
    period_dispatch = cell(0, 1);
    for k = 1:numel(theta_list)
        theta = theta_list(k);
        fprintf('Solving IEEE-13 LinDistFlow fixed_angle, t=%d/%d, theta=%.2f deg, idx=%d/%d ...\n', ...
            t, params.T, rad2deg(theta), k, numel(theta_list)); drawnow;
        soln = solve_pq_boundary_point_lindisflow(t, theta, sys, ts, params);
        row = build_result_row(params.SCAN_MODE, t, theta, soln, 0, true, NaN, NaN);
        period_rows = append_struct(period_rows, row);
        period_rows = mark_duplicate_points(period_rows, params.DUP_POINT_TOL);
        period_dispatch{end+1, 1} = collect_dispatch_solution(soln, sys, ts, t); %#ok<AGROW>
    end
    period_summary = summarize_period(t, period_rows, params, 0);
end

function [period_rows, period_dispatch, period_summary] = run_adaptive_scan_for_period(t, sys, ts, params)
    fprintf('Adaptive scan initial directions: %d, max_rounds=%d, max_solves=%d\n', ...
        params.ADAPTIVE_INITIAL_N_THETA, params.ADAPTIVE_MAX_ROUNDS, params.ADAPTIVE_MAX_SOLVES); drawnow;
    period_rows = struct([]);
    period_dispatch = cell(0, 1);
    solved_theta = [];
    pending_angles = get_fixed_angles(params.ADAPTIVE_INITIAL_N_THETA);
    pending_left = NaN(size(pending_angles));
    pending_right = NaN(size(pending_angles));
    pending_initial = true(size(pending_angles));
    rounds_used = 0;
    stop_reason = 'maximum rounds reached';

    for round_id = 1:params.ADAPTIVE_MAX_ROUNDS
        rounds_used = round_id;
        solve_count_before = numel(solved_theta);
        for k = 1:numel(pending_angles)
            if numel(solved_theta) >= params.ADAPTIVE_MAX_SOLVES
                stop_reason = 'reached ADAPTIVE_MAX_SOLVES before solving all pending angles';
                break;
            end
            theta = mod(pending_angles(k), 2*pi);
            if is_angle_already_solved(theta, solved_theta, params.ANGLE_TOL)
                fprintf('Skipping duplicate angle theta=%.8f deg in adaptive round %d.\n', rad2deg(theta), round_id); drawnow;
                continue;
            end
            fprintf('Solving IEEE-13 LinDistFlow adaptive_scan, t=%d/%d, round=%d/%d, theta=%.2f deg, idx=%d/%d, solved=%d/%d ...\n', ...
                t, params.T, round_id, params.ADAPTIVE_MAX_ROUNDS, rad2deg(theta), k, numel(pending_angles), numel(solved_theta)+1, params.ADAPTIVE_MAX_SOLVES); drawnow;
            soln = solve_pq_boundary_point_lindisflow(t, theta, sys, ts, params);
            row = build_result_row(params.SCAN_MODE, t, theta, soln, round_id, pending_initial(k), pending_left(k), pending_right(k));
            period_rows = append_struct(period_rows, row);
            period_dispatch{end+1, 1} = collect_dispatch_solution(soln, sys, ts, t); %#ok<AGROW>
            solved_theta(end+1) = theta; %#ok<AGROW>
        end
        period_rows = mark_duplicate_points(period_rows, params.DUP_POINT_TOL);
        temp_summary = summarize_period(t, period_rows, params, rounds_used);
        if numel(solved_theta) >= params.ADAPTIVE_MAX_SOLVES
            stop_reason = 'reached ADAPTIVE_MAX_SOLVES';
            fprintf('Adaptive round %d finished: total_solved=%d, unique_points=%d, new_angles=0, area=%.6f\n', ...
                round_id, numel(solved_theta), temp_summary.num_unique_boundary_points, temp_summary.polygon_area); drawnow;
            break;
        end
        [new_angles, parent_left, parent_right] = propose_new_angles_adaptive(period_rows, params);
        [new_angles, keep] = filter_unsolved_angles(new_angles, solved_theta, params.ANGLE_TOL);
        parent_left = parent_left(keep);
        parent_right = parent_right(keep);
        remaining_slots = params.ADAPTIVE_MAX_SOLVES - numel(solved_theta);
        if numel(new_angles) > remaining_slots
            new_angles = new_angles(1:remaining_slots);
            parent_left = parent_left(1:remaining_slots);
            parent_right = parent_right(1:remaining_slots);
        end
        fprintf('Adaptive round %d finished: total_solved=%d, unique_points=%d, new_angles=%d, area=%.6f\n', ...
            round_id, numel(solved_theta), temp_summary.num_unique_boundary_points, numel(new_angles), temp_summary.polygon_area); drawnow;
        if isempty(new_angles)
            stop_reason = 'all adjacent boundary segments satisfy adaptive tolerances or only duplicate-point segments remain';
            break;
        end
        if numel(new_angles) < params.ADAPTIVE_MIN_NEW_ANGLES
            stop_reason = 'new angle count below ADAPTIVE_MIN_NEW_ANGLES';
            break;
        end
        pending_angles = new_angles;
        pending_left = parent_left;
        pending_right = parent_right;
        pending_initial = false(size(pending_angles));
        if numel(solved_theta) == solve_count_before && isempty(pending_angles)
            stop_reason = 'no additional angle was solved in this round';
            break;
        end
    end
    period_rows = mark_duplicate_points(period_rows, params.DUP_POINT_TOL);
    period_summary = summarize_period(t, period_rows, params, rounds_used);
    fprintf('Adaptive scan stopped for t=%d: %s. total_solved=%d, unique_points=%d, area=%.6f\n', ...
        t, stop_reason, period_summary.num_total_solves, period_summary.num_unique_boundary_points, period_summary.polygon_area); drawnow;
end

function [new_angles, parent_left_deg, parent_right_deg] = propose_new_angles_adaptive(rows, params)
    new_angles = [];
    parent_left_deg = [];
    parent_right_deg = [];
    if isempty(rows), return; end
    tbl = struct2table(rows);
    ok = tbl.solver_status == 0 & ~isnan(tbl.P0) & ~isnan(tbl.Q0);
    tbl = tbl(ok, :);
    if height(tbl) < 2, return; end
    tbl = sortrows(tbl, 'theta_rad');
    P = tbl.P0;
    Q = tbl.Q0;
    theta = tbl.theta_rad;
    P_range = max(max(P) - min(P), 1e-6);
    Q_range = max(max(Q) - min(Q), 1e-6);
    n = height(tbl);
    for i = 1:n
        j = i + 1;
        if j > n, j = 1; end
        if abs(P(i) - P(j)) <= params.DUP_POINT_TOL && abs(Q(i) - Q(j)) <= params.DUP_POINT_TOL
            continue;
        end
        dP = abs((P(i) - P(j)) / P_range);
        dQ = abs((Q(i) - Q(j)) / Q_range);
        d = hypot(dP, dQ);
        if d > params.ADAPTIVE_DIST_TOL || dP > params.ADAPTIVE_P_TOL || dQ > params.ADAPTIVE_Q_TOL
            th_mid = angle_mid_circular(theta(i), theta(j));
            new_angles(end+1) = th_mid; %#ok<AGROW>
            parent_left_deg(end+1) = rad2deg(theta(i)); %#ok<AGROW>
            parent_right_deg(end+1) = rad2deg(theta(j)); %#ok<AGROW>
        end
    end
    [new_angles, unique_idx] = unique_angles_keep_order(new_angles, params.ANGLE_TOL);
    parent_left_deg = parent_left_deg(unique_idx);
    parent_right_deg = parent_right_deg(unique_idx);
end

function theta_mid = angle_mid_circular(theta_left, theta_right)
    theta_left = mod(theta_left, 2*pi);
    theta_right = mod(theta_right, 2*pi);
    if theta_right <= theta_left
        theta_right = theta_right + 2*pi;
    end
    theta_mid = mod((theta_left + theta_right) / 2, 2*pi);
end

function [P_unique, Q_unique, unique_idx, is_duplicate] = deduplicate_boundary_points(P, Q, tol)
    P = P(:);
    Q = Q(:);
    unique_idx = [];
    is_duplicate = false(numel(P), 1);
    for i = 1:numel(P)
        if isnan(P(i)) || isnan(Q(i))
            is_duplicate(i) = true;
            continue;
        end
        duplicate_found = false;
        for k = 1:numel(unique_idx)
            j = unique_idx(k);
            if abs(P(i) - P(j)) <= tol && abs(Q(i) - Q(j)) <= tol
                duplicate_found = true;
                break;
            end
        end
        if duplicate_found
            is_duplicate(i) = true;
        else
            unique_idx(end+1, 1) = i; %#ok<AGROW>
        end
    end
    P_unique = P(unique_idx);
    Q_unique = Q(unique_idx);
end

function rows = mark_duplicate_points(rows, tol)
    if isempty(rows), return; end
    tbl = struct2table(rows);
    ok = tbl.solver_status == 0 & ~isnan(tbl.P0) & ~isnan(tbl.Q0);
    dup_all = false(height(tbl), 1);
    if any(ok)
        ok_idx = find(ok);
        [~, ~, ~, dup_ok] = deduplicate_boundary_points(tbl.P0(ok), tbl.Q0(ok), tol);
        dup_all(ok_idx) = dup_ok;
    end
    for i = 1:numel(rows)
        rows(i).is_duplicate_point = dup_all(i);
    end
end

function tf = is_angle_already_solved(theta, solved_theta, tol)
    tf = false;
    for i = 1:numel(solved_theta)
        diff_val = abs(angle_diff_circular(theta, solved_theta(i)));
        if diff_val < tol
            tf = true;
            return;
        end
    end
end

function d = angle_diff_circular(a, b)
    d = mod(a - b + pi, 2*pi) - pi;
end

function [angles_out, keep_idx] = filter_unsolved_angles(angles_in, solved_theta, tol)
    keep_idx = [];
    angles_out = [];
    [angles_unique, unique_idx] = unique_angles_keep_order(angles_in, tol);
    for i = 1:numel(angles_unique)
        if ~is_angle_already_solved(angles_unique(i), solved_theta, tol)
            angles_out(end+1) = angles_unique(i); %#ok<AGROW>
            keep_idx(end+1) = unique_idx(i); %#ok<AGROW>
        end
    end
end

function [angles_unique, unique_idx] = unique_angles_keep_order(angles_in, tol)
    angles_unique = [];
    unique_idx = [];
    for i = 1:numel(angles_in)
        th = mod(angles_in(i), 2*pi);
        if ~is_angle_already_solved(th, angles_unique, tol)
            angles_unique(end+1) = th; %#ok<AGROW>
            unique_idx(end+1) = i; %#ok<AGROW>
        end
    end
end

function soln = solve_pq_boundary_point_lindisflow(t, theta, sys, ts, params)
    alpha = cos(theta);
    beta = sin(theta);
    nb = sys.nb;
    nl = sys.nl;

    Pbr = sdpvar(nl, 1, 'full');
    Qbr = sdpvar(nl, 1, 'full');
    v = sdpvar(nb, 1, 'full');
    PG = sdpvar(1, 1);
    QG = sdpvar(1, 1);
    PPV = sdpvar(1, 1);
    QPV = sdpvar(1, 1);
    Pch = sdpvar(1, 1);
    Pdis = sdpvar(1, 1);
    QESS = sdpvar(1, 1);
    if params.USE_ESS_BINARY
        zch = binvar(1, 1);
        zdis = binvar(1, 1);
    end
    PESS = Pdis - Pch;
    Proot = Pbr(sys.root_branch);
    Qroot = Qbr(sys.root_branch);
    P0 = -Proot;
    Q0 = -Qroot;

    cons = [];
    cons = [cons, v(sys.root_bus) == sys.V0^2];
    cons = [cons, params.VMIN <= v, v <= params.VMAX];

    % Device constraints.
    cons = [cons, sys.PG_min <= PG, PG <= sys.PG_max, sys.QG_min <= QG, QG <= sys.QG_max];
    cons = [cons, 0 <= PPV, PPV <= ts.PV_AVAIL_PROFILE(t), -sys.QPV_max <= QPV, QPV <= sys.QPV_max];
    if params.USE_PV_SOC
        cons = [cons, cone([PPV; QPV], sys.SPV)]; %#ok<AGROW>
    end
    cons = [cons, 0 <= Pch, Pch <= sys.PCH_max, 0 <= Pdis, Pdis <= sys.PDIS_max];
    cons = [cons, -sys.QESS_max <= QESS, QESS <= sys.QESS_max];
    cons = [cons, -sys.SESS <= PESS, PESS <= sys.SESS];
    if params.USE_ESS_SOC_CONE
        cons = [cons, cone([PESS; QESS], sys.SESS)]; %#ok<AGROW>
    end
    if params.USE_ESS_BINARY
        cons = [cons, Pch <= sys.PCH_max * zch, Pdis <= sys.PDIS_max * zdis, zch + zdis <= 1]; %#ok<AGROW>
    end

    % LinDistFlow constraints on each radial branch i -> j.
    for e = 1:nl
        i = sys.from(e);
        j = sys.to(e);
        child_edges = sys.children_branches{j};
        [p_gen_j, q_gen_j, p_pv_j, q_pv_j, p_ess_j, q_ess_j] = device_injection_at_bus(j, sys, PG, QG, PPV, QPV, PESS, QESS);
        cons = [cons, Pbr(e) == sum(Pbr(child_edges)) + ts.Pload(j, t) - p_gen_j - p_pv_j - p_ess_j]; %#ok<AGROW>
        cons = [cons, Qbr(e) == sum(Qbr(child_edges)) + ts.Qload(j, t) - q_gen_j - q_pv_j - q_ess_j]; %#ok<AGROW>
        cons = [cons, v(j) == v(i) - 2 * (sys.r(e) * Pbr(e) / params.baseMVA + sys.x(e) * Qbr(e) / params.baseMVA)]; %#ok<AGROW>
        cons = [cons, -sys.Smax(e) <= Pbr(e), Pbr(e) <= sys.Smax(e)]; %#ok<AGROW>
        cons = [cons, -sys.Smax(e) <= Qbr(e), Qbr(e) <= sys.Smax(e)]; %#ok<AGROW>
        if params.USE_BRANCH_SOC
            cons = [cons, cone([Pbr(e); Qbr(e)], sys.Smax(e))]; %#ok<AGROW>
        end
    end

    obj_expr = alpha * P0 + beta * Q0;
    ops = sdpsettings('solver', 'gurobi', 'verbose', params.YALMIP_VERBOSE, 'gurobi.OutputFlag', params.GUROBI_VERBOSE);

    soln = init_solution_struct(t, theta, alpha, beta, nb, nl);
    try
        sol = optimize(cons, -obj_expr, ops);
        soln.solver_status = sol.problem;
        soln.solver_info = sol.info;
    catch ME
        soln.solver_status = 100;
        soln.solver_info = sprintf('optimize exception: %s', ME.message);
        fprintf('  Optimization exception: %s\n', ME.message); drawnow;
        return;
    end

    if soln.solver_status ~= 0
        fprintf('  Optimization failed: problem=%d, info=%s\n', soln.solver_status, soln.solver_info); drawnow;
        return;
    end

    vals = [value(P0); value(Q0); value(obj_expr); value(PG); value(QG); value(PPV); value(QPV); value(Pch); value(Pdis); value(PESS); value(QESS); value(Pbr); value(Qbr); value(v)];
    if any(isnan(vals))
        soln.solver_status = 99;
        soln.solver_info = 'YALMIP returned NaN values';
        fprintf('  Optimization returned NaN values; marking as failed.\n'); drawnow;
        return;
    end

    soln.P0 = value(P0);
    soln.Q0 = value(Q0);
    soln.objective = value(obj_expr);
    soln.PG = value(PG);
    soln.QG = value(QG);
    soln.PPV = value(PPV);
    soln.QPV = value(QPV);
    soln.Pch = value(Pch);
    soln.Pdis = value(Pdis);
    soln.PESS = value(PESS);
    soln.QESS = value(QESS);
    soln.Pbr = value(Pbr);
    soln.Qbr = value(Qbr);
    soln.v = value(v);
    [soln.max_voltage_violation, soln.max_branch_p_violation, soln.max_branch_q_violation] = calc_model_violations(soln, sys, params);
    fprintf('  Success: P0=%.4f MW, Q0=%.4f MVAr, objective=%.4f\n', soln.P0, soln.Q0, soln.objective); drawnow;
end


function [p_gen, q_gen, p_pv, q_pv, p_ess, q_ess] = device_injection_at_bus(bus_idx, sys, PG, QG, PPV, QPV, PESS, QESS)
    p_gen = 0; q_gen = 0; p_pv = 0; q_pv = 0; p_ess = 0; q_ess = 0;
    if bus_idx == sys.gen_bus
        p_gen = PG;
        q_gen = QG;
    end
    if bus_idx == sys.pv_bus
        p_pv = PPV;
        q_pv = QPV;
    end
    if bus_idx == sys.ess_bus
        p_ess = PESS;
        q_ess = QESS;
    end
end

function soln = init_solution_struct(t, theta, alpha, beta, nb, nl)
    soln = struct();
    soln.t = t;
    soln.theta_deg = rad2deg(theta);
    soln.alpha = alpha;
    soln.beta = beta;
    soln.P0 = NaN;
    soln.Q0 = NaN;
    soln.objective = NaN;
    soln.solver_status = NaN;
    soln.solver_info = '';
    soln.PG = NaN;
    soln.QG = NaN;
    soln.PPV = NaN;
    soln.QPV = NaN;
    soln.Pch = NaN;
    soln.Pdis = NaN;
    soln.PESS = NaN;
    soln.QESS = NaN;
    soln.Pbr = NaN(nl, 1);
    soln.Qbr = NaN(nl, 1);
    soln.v = NaN(nb, 1);
    soln.max_voltage_violation = NaN;
    soln.max_branch_p_violation = NaN;
    soln.max_branch_q_violation = NaN;
end

function [max_v_viol, max_p_viol, max_q_viol] = calc_model_violations(soln, sys, params)
    v_low = max(params.VMIN - soln.v, 0);
    v_high = max(soln.v - params.VMAX, 0);
    max_v_viol = max([v_low; v_high]);
    max_p_viol = max(max(abs(soln.Pbr) - sys.Smax, 0));
    max_q_viol = max(max(abs(soln.Qbr) - sys.Smax, 0));
end

function row = build_result_row(scan_mode, t, theta, soln, adaptive_round, is_initial_angle, parent_theta_left_deg, parent_theta_right_deg)
    row = struct();
    row.scan_mode = scan_mode;
    row.t = t;
    row.adaptive_round = adaptive_round;
    row.theta_rad = mod(theta, 2*pi);
    row.theta_deg = rad2deg(mod(theta, 2*pi));
    row.alpha = cos(theta);
    row.beta = sin(theta);
    row.P0 = soln.P0;
    row.Q0 = soln.Q0;
    row.objective = soln.objective;
    row.solver_status = soln.solver_status;
    row.solver_info = soln.solver_info;
    row.is_initial_angle = logical(is_initial_angle);
    row.parent_theta_left_deg = parent_theta_left_deg;
    row.parent_theta_right_deg = parent_theta_right_deg;
    row.is_duplicate_point = false;
    row.max_voltage_violation = soln.max_voltage_violation;
    row.max_branch_p_violation = soln.max_branch_p_violation;
    row.max_branch_q_violation = soln.max_branch_q_violation;
end

function dispatch = collect_dispatch_solution(soln, sys, ts, t)
    dispatch = struct();
    dispatch.t = t;
    dispatch.theta_deg = soln.theta_deg;
    dispatch.bus_name = sys.bus_name;
    dispatch.branch_name = sys.branch_name;
    dispatch.Pload = ts.Pload(:, t);
    dispatch.Qload = ts.Qload(:, t);
    dispatch.PG = soln.PG;
    dispatch.QG = soln.QG;
    dispatch.PPV = soln.PPV;
    dispatch.QPV = soln.QPV;
    dispatch.Pch = soln.Pch;
    dispatch.Pdis = soln.Pdis;
    dispatch.PESS = soln.PESS;
    dispatch.QESS = soln.QESS;
    dispatch.Pbr = soln.Pbr;
    dispatch.Qbr = soln.Qbr;
    dispatch.v = soln.v;
    dispatch.P0 = soln.P0;
    dispatch.Q0 = soln.Q0;
    dispatch.objective = soln.objective;
    dispatch.solver_status = soln.solver_status;
    dispatch.solver_info = soln.solver_info;
end

function summary = summarize_period(t, rows, params, adaptive_rounds_used)
    summary = struct('scan_mode', params.SCAN_MODE, 't', t, 'num_total_solves', 0, ...
        'num_successful_solves', 0, 'num_unique_boundary_points', 0, 'num_duplicate_points', 0, ...
        'num_failed_optimizations', 0, 'adaptive_rounds_used', adaptive_rounds_used, 'polygon_area', NaN);
    if isempty(rows)
        return;
    end
    rows = mark_duplicate_points(rows, params.DUP_POINT_TOL);
    tbl = struct2table(rows);
    ok = tbl.solver_status == 0 & ~isnan(tbl.P0) & ~isnan(tbl.Q0);
    summary.num_total_solves = height(tbl);
    summary.num_successful_solves = sum(ok);
    summary.num_failed_optimizations = sum(tbl.solver_status ~= 0);
    summary.num_duplicate_points = sum(tbl.is_duplicate_point & ok);
    if any(ok)
        [P_unique, Q_unique] = deduplicate_boundary_points(tbl.P0(ok), tbl.Q0(ok), params.DUP_POINT_TOL);
        summary.num_unique_boundary_points = numel(P_unique);
        if summary.num_unique_boundary_points >= 3
            [P, Q] = sort_boundary_points(P_unique, Q_unique);
            summary.polygon_area = polygon_area_pq(P, Q);
        end
    end
    fprintf('Summary mode=%s, t=%d: total=%d, success=%d, unique=%d, duplicate=%d, failed=%d, rounds=%d, polygon_area=%.6f\n', ...
        params.SCAN_MODE, t, summary.num_total_solves, summary.num_successful_solves, ...
        summary.num_unique_boundary_points, summary.num_duplicate_points, summary.num_failed_optimizations, ...
        summary.adaptive_rounds_used, summary.polygon_area); drawnow;
end

function [Psort, Qsort, order] = sort_boundary_points(P, Q)
    P = P(:);
    Q = Q(:);
    cP = mean(P, 'omitnan');
    cQ = mean(Q, 'omitnan');
    ang = atan2(Q - cQ, P - cP);
    [~, order] = sort(ang);
    Psort = P(order);
    Qsort = Q(order);
end

function area = polygon_area_pq(P, Q)
    P = P(:);
    Q = Q(:);
    if numel(P) < 3 || any(isnan(P)) || any(isnan(Q))
        area = NaN;
        return;
    end
    area = polyarea(P, Q);
end

function plot_pq_flex_region(rows, t, params, summary)
    if isempty(rows), return; end
    rows = mark_duplicate_points(rows, params.DUP_POINT_TOL);
    tbl = struct2table(rows);
    ok = tbl.solver_status == 0 & ~isnan(tbl.P0) & ~isnan(tbl.Q0);
    if sum(ok) < 1
        fprintf('No successful points for t=%d; skip plot.\n', t); drawnow;
        return;
    end
    P_all = tbl.P0(ok);
    Q_all = tbl.Q0(ok);
    [P_unique, Q_unique] = deduplicate_boundary_points(P_all, Q_all, params.DUP_POINT_TOL);
    area = NaN;
    if numel(P_unique) >= 3
        [Ps, Qs] = sort_boundary_points(P_unique, Q_unique);
        area = polygon_area_pq(Ps, Qs);
    else
        Ps = P_unique;
        Qs = Q_unique;
    end

    fig = figure('Visible', 'off', 'Color', 'w');
    hold on; grid on; box on;
    if numel(Ps) >= 3
        fill([Ps; Ps(1)], [Qs; Qs(1)], [0.84 0.92 1.00], 'EdgeColor', [0.00 0.25 0.80], ...
            'LineWidth', 1.5, 'FaceAlpha', 0.60, 'DisplayName', 'P-Q polygon');
    end
    plot(P_all, Q_all, '.', 'Color', [0.55 0.55 0.55], 'MarkerSize', 12, 'DisplayName', 'all solved points');
    plot(P_unique, Q_unique, 'o', 'Color', [0.00 0.20 0.75], 'MarkerFaceColor', [0.00 0.45 0.95], ...
        'MarkerSize', 6, 'DisplayName', 'unique boundary points');
    xlabel('P_0 export to upstream grid (MW)');
    ylabel('Q_0 export to upstream grid (MVAr)');
    if strcmpi(params.SCAN_MODE, 'adaptive_scan')
        title(sprintf('IEEE-13 simplified balanced feeder | LinDistFlow | single_period | adaptive_scan | initial_N=%d | total_solves=%d | unique_points=%d | area=%.6f', ...
            params.ADAPTIVE_INITIAL_N_THETA, summary.num_total_solves, summary.num_unique_boundary_points, area));
        fname = fullfile(params.RESULT_DIR, sprintf('pq_flex_t%02d_ieee13_lindisflow_adaptive.png', t));
    else
        title(sprintf('IEEE-13 simplified balanced feeder | LinDistFlow | single_period | fixed_angle | Ntheta=%d | total_solves=%d | unique_points=%d | area=%.6f', ...
            params.FIXED_N_THETA, summary.num_total_solves, summary.num_unique_boundary_points, area));
        fname = fullfile(params.RESULT_DIR, sprintf('pq_flex_t%02d_ieee13_lindisflow_fixed.png', t));
    end
    legend('Location', 'best');
    save_png_compat(fig, fname);
    close(fig);
    fprintf('Saved figure: %s\n', fname); drawnow;
end

function save_png_compat(fig, fname)
    try
        exportgraphics(fig, fname, 'Resolution', 200);
    catch
        set(fig, 'PaperPositionMode', 'auto');
        print(fig, fname, '-dpng', '-r200');
    end
end

function save_results_csv(all_rows, summary_rows, result_dir)
    if isempty(all_rows)
        warning('No boundary rows to save.');
    else
        writetable(struct2table(all_rows), fullfile(result_dir, 'boundary_points_all.csv'));
        fprintf('Saved CSV: %s\n', fullfile(result_dir, 'boundary_points_all.csv')); drawnow;
    end
    if isempty(summary_rows)
        warning('No summary rows to save.');
    else
        writetable(struct2table(summary_rows), fullfile(result_dir, 'summary_area.csv'));
        fprintf('Saved CSV: %s\n', fullfile(result_dir, 'summary_area.csv')); drawnow;
    end
end

function out = append_struct(a, b)
    if isempty(b)
        out = a;
        return;
    end
    if isempty(a)
        out = b;
    else
        out = [a(:); b(:)];
    end
end
