%% ieee13_temp_p0_v1_2_decoup.m
% IEEE-13 simplified balanced-feeder temporal P0 trajectory-boundary Step 1.
%
% This v1_2_decoup version keeps the Step 1 temporal P0 boundary workflow,
% adds a matrix form of the LinDistFlow network constraints, validates that
% matrix form against original dispatch solutions, and then builds conservative
% VG/VB robust decoupled network limits. Coupled VG/VB projection bounds are
% full-network projection bounds and need not be independently composable;
% decoupled VG/VB robust bounds subtract worst-case network occupation by the
% other resource class, so they are conservative and more suitable as a first
% independent VG/VB equivalent. The robust margins here use only single-period
% resource boxes, not the other side's ramp/SOC coupling, and this is not bound
% shrinking or infeasible-trajectory checking.
%
% This script uses a simplified single-phase balanced equivalent of the IEEE
% 13-node distribution feeder with a multi-period LinDistFlow model. It solves
% the cross-time Step 1 boundary problems for total PCC active-power trajectory
% limits, namely per-period P0 power bounds, adjacent-period P0 ramp bounds, and
% cumulative P0 exchange-energy bounds.
%
% Scope and sign convention:
%   - The feeder is a simplified single-phase balanced LinDistFlow model, not
%     the original three-phase unbalanced IEEE-13 feeder.
%   - This version considers G1 active-power ramping and ESS SOC recursion.
%   - This version does NOT perform single-period P-Q boundary scanning.
%   - This version implements a first conservative VG/VB robust decoupling
%     based on matrixized network constraints and box-based worst-case margins.
%   - This version does NOT implement bound shrinking, Farkas/KKT infeasible-point
%     search, undecomposable trajectory search, or neural-network training.
%   - This is only a matrixization/robust-decoupling prototype and is not the
%     complete final stage of the paper workflow.
%   - Root branch is 650 -> 632. P_root > 0 means importing active power from
%     the upstream grid into the feeder. Reported P0 = -P_root, so P0 > 0 means
%     active-power export from the feeder to the upstream grid, and P0 < 0 means
%     feeder import/purchase from the upstream grid. Q0 is defined similarly.
%
% ESS convention:
%   PESS(t) > 0 means ESS discharge/injection to the feeder; PESS(t) < 0 means
%   charging/absorbing from the feeder. This LP version uses net PESS directly,
%   without Pch/Pdis, efficiency, or charge/discharge binary variables.

clear; clc; close all;

%% --------------------------- User parameters ---------------------------
RESULT_DIR = 'results_ieee13_temp_p0_v1_2_decoup';

T = 4;
DT = 1.0;
LOAD_SCALE = [0.90, 1.00, 1.10, 0.95];
PV_AVAIL_PROFILE = [1.50, 1.80, 2.10, 1.60];  % MW, enhanced to allow P0 to cross zero

baseMVA = 5;
V0 = 1.0;
VMIN = 0.95^2;
VMAX = 1.05^2;

% First temporal version is an LP by default. SOCP/binary switches are kept for
% future experiments but disabled here.
USE_BRANCH_SOC = false;
USE_PV_SOC = false;
USE_ESS_SOC_CONE = false;
USE_ESS_BINARY = false; %#ok<NASGU>

YALMIP_VERBOSE = 0;
GUROBI_VERBOSE = 0;

% Device parameters.
PG_min = 0.00;
PG_max = 1.20;
QG_min = -0.60;
QG_max = 0.60;
RG_up = 0.40;          % MW per time step
RG_down = 0.40;        % MW per time step
PG0 = 0.60;            % initial G1 active power for t=1 ramp constraint

SPV = 2.30;
QPV_max = 1.00;

PESS_ch_max = 0.80;    % MW, PESS >= -PESS_ch_max
PESS_dis_max = 0.80;   % MW, PESS <= PESS_dis_max
QESS_max = 0.50;
SESS = 1.00;
E_min = 0.20;          % MWh
E_max = 2.00;          % MWh
E0 = 1.00;             % MWh
ENFORCE_TERMINAL_SOC = true;

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

fprintf('Running file: ieee13_temp_p0_v1_2_decoup.m\n'); drawnow;
fprintf('IEEE-13 temporal P0 boundary LinDistFlow Step 1 + VG/VB decoupling started at %s\n', datestr(now)); drawnow;
fprintf('Result directory: %s\n', RESULT_DIR); drawnow;
fprintf('T=%d, DT=%.3f h, baseMVA=%.3f, terminal_SOC=%d\n', T, DT, baseMVA, ENFORCE_TERMINAL_SOC); drawnow;

check_dependencies();
sys = build_ieee13_simplified_system(baseMVA, V0, VMIN, VMAX);
ts = build_time_series_data(sys, T, LOAD_SCALE, PV_AVAIL_PROFILE);

params = struct();
params.RESULT_DIR = RESULT_DIR;
params.T = T;
params.DT = DT;
params.baseMVA = baseMVA;
params.V0 = V0;
params.VMIN = VMIN;
params.VMAX = VMAX;
params.USE_BRANCH_SOC = USE_BRANCH_SOC;
params.USE_PV_SOC = USE_PV_SOC;
params.USE_ESS_SOC_CONE = USE_ESS_SOC_CONE;
params.YALMIP_VERBOSE = YALMIP_VERBOSE;
params.GUROBI_VERBOSE = GUROBI_VERBOSE;
params.PG_min = PG_min;
params.PG_max = PG_max;
params.QG_min = QG_min;
params.QG_max = QG_max;
params.RG_up = RG_up;
params.RG_down = RG_down;
params.PG0 = PG0;
params.SPV = SPV;
params.QPV_max = QPV_max;
params.PESS_ch_max = PESS_ch_max;
params.PESS_dis_max = PESS_dis_max;
params.QESS_max = QESS_max;
params.SESS = SESS;
params.E_min = E_min;
params.E_max = E_max;
params.E0 = E0;
params.ENFORCE_TERMINAL_SOC = ENFORCE_TERMINAL_SOC;

%% ------------------------- Boundary solves -----------------------------
all_solutions = struct([]);

fprintf('\n=== Solving per-period P0 power bounds ===\n'); drawnow;
for t = 1:T
    for isense = 1:2
        sense = pick_sense(isense);
        fprintf('Solving temporal boundary: type=P0_power, sense=%s, target_t=%d/%d ...\n', sense, t, T); drawnow;
        soln = solve_temporal_boundary_problem('P0_power', t, sense, sys, ts, params);
        all_solutions = append_struct(all_solutions, soln);
    end
end

fprintf('\n=== Solving adjacent-period P0 ramp bounds ===\n'); drawnow;
for t = 2:T
    for isense = 1:2
        sense = pick_sense(isense);
        fprintf('Solving temporal boundary: type=P0_ramp, sense=%s, target_t=%d/%d ...\n', sense, t, T); drawnow;
        soln = solve_temporal_boundary_problem('P0_ramp', t, sense, sys, ts, params);
        all_solutions = append_struct(all_solutions, soln);
    end
end

fprintf('\n=== Solving cumulative P0 exchange-energy bounds ===\n'); drawnow;
for t = 1:T
    for isense = 1:2
        sense = pick_sense(isense);
        fprintf('Solving temporal boundary: type=P0_cum_energy, sense=%s, target_t=%d/%d ...\n', sense, t, T); drawnow;
        soln = solve_temporal_boundary_problem('P0_cum_energy', t, sense, sys, ts, params);
        all_solutions = append_struct(all_solutions, soln);
    end
end

temporal_bounds = assemble_temporal_bounds_table(all_solutions, T);

%% ------------------------- Save and plot -------------------------------
save_temporal_results_csv(all_solutions, temporal_bounds, RESULT_DIR, params);
save(fullfile(RESULT_DIR, 'temporal_boundary_dispatch.mat'), 'all_solutions', 'temporal_bounds', 'sys', 'ts', 'params', '-v7.3');
fprintf('Saved MAT: %s\n', fullfile(RESULT_DIR, 'temporal_boundary_dispatch.mat')); drawnow;

plot_p0_power_bounds(temporal_bounds, params);
plot_p0_ramp_bounds(temporal_bounds, params);
plot_p0_cum_energy_bounds(temporal_bounds, params);
plot_ess_soc_examples(all_solutions, params);
plot_all_boundary_p0_trajectories(all_solutions, temporal_bounds, params);


%% ---------------- Matrix validation and VG/VB decoupling ---------------
all_solutions = annotate_vgvb_decomposition(all_solutions, ts, params);
netmat = build_network_matrices_lindisflow(sys, ts, params);
matrix_validation = validate_network_matrix_against_original(all_solutions, netmat, sys, ts, params);
save_network_matrix_validation_csv(matrix_validation, RESULT_DIR);
if ~matrix_validation.is_passed
    error('Network matrix validation failed; skip VG/VB decoupling.');
end

decoup = compute_robust_decoupled_limits(netmat, ts, params);
save_decoupled_network_margins_csv(decoup, netmat, RESULT_DIR);

fprintf('\n=== Solving coupled VG/VB projection bounds ===\n'); drawnow;
coupled_vgvb_solutions = solve_all_coupled_vgvb_projection_bounds(sys, ts, params);
coupled_vgvb_bounds = assemble_vgvb_bounds_table(coupled_vgvb_solutions, params.T);
save_vgvb_bounds_csv(coupled_vgvb_bounds, fullfile(RESULT_DIR, 'coupled_vgvb_projection_bounds.csv'));

fprintf('\n=== Solving decoupled VG/VB robust bounds ===\n'); drawnow;
decoupled_vgvb_solutions = solve_all_decoupled_vgvb_bounds(decoup, sys, ts, params);
decoupled_vgvb_bounds = assemble_vgvb_bounds_table(decoupled_vgvb_solutions, params.T);
save_vgvb_bounds_csv(decoupled_vgvb_bounds, fullfile(RESULT_DIR, 'decoupled_vgvb_bounds.csv'));
save_vgvb_solutions_csv(decoupled_vgvb_solutions, fullfile(RESULT_DIR, 'decoupled_vgvb_solutions.csv'), params);

comparison = assemble_vgvb_comparison(coupled_vgvb_bounds, decoupled_vgvb_bounds, params.T);
writetable(struct2table(comparison), fullfile(RESULT_DIR, 'vgvb_decoupling_comparison.csv'));
fprintf('Saved CSV: %s\n', fullfile(RESULT_DIR, 'vgvb_decoupling_comparison.csv')); drawnow;

plot_vgvb_power_bounds_compare(coupled_vgvb_bounds, decoupled_vgvb_bounds, params, 'VG');
plot_vg_ramp_bounds_compare(coupled_vgvb_bounds, decoupled_vgvb_bounds, params);
plot_vgvb_power_bounds_compare(coupled_vgvb_bounds, decoupled_vgvb_bounds, params, 'VB');
plot_vb_energy_bounds_compare(coupled_vgvb_bounds, decoupled_vgvb_bounds, params);
plot_network_decoupling_margin_summary(decoup, params);

save(fullfile(RESULT_DIR, 'decoupled_vgvb_dispatch.mat'), 'netmat', 'decoup', 'matrix_validation', ...
    'all_solutions', 'coupled_vgvb_solutions', 'decoupled_vgvb_solutions', ...
    'coupled_vgvb_bounds', 'decoupled_vgvb_bounds', 'comparison', 'sys', 'ts', 'params', '-v7.3');
fprintf('Saved MAT: %s\n', fullfile(RESULT_DIR, 'decoupled_vgvb_dispatch.mat')); drawnow;

fprintf('\nCompleted IEEE-13 temporal P0 boundary Step 1 at %s\n', datestr(now)); drawnow;
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
        error('Please install/add the missing dependencies before running this script. MATPOWER is not required.');
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

    % Simplified positive-sequence-equivalent line parameters. These are
    % deliberately approximate and are intended for Step 1 workflow validation.
    r = [0.010; 0.006; 0.006; 0.008; 0.006; 0.012; 0.004; 0.007; 0.005; 0.005; 0.006; 0.007];
    x = [0.030; 0.018; 0.018; 0.024; 0.018; 0.036; 0.012; 0.021; 0.015; 0.015; 0.018; 0.021];
    Smax = [4.0; 2.5; 1.8; 1.6; 1.3; 3.5; 1.0; 1.8; 0.8; 0.8; 1.8; 1.6];

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

    fprintf('Built IEEE-13 simplified balanced feeder with %d buses and %d radial branches.\n', nb, nl); drawnow;
    fprintf('Root/PCC=%s, root branch=%s, G1=%s, PV=%s, ESS=%s.\n', ...
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

function sense = pick_sense(isense)
    if isense == 1
        sense = 'min';
    else
        sense = 'max';
    end
end

function soln = solve_temporal_boundary_problem(objective_type, target_t, sense, sys, ts, params)
    if strcmpi(objective_type, 'P0_ramp') && target_t < 2
        error('P0_ramp objective requires target_t >= 2.');
    end
    [cons, vars] = build_temporal_lindisflow_constraints(sys, ts, params);
    switch objective_type
        case 'P0_power'
            obj = vars.P0(target_t);
        case 'P0_ramp'
            obj = vars.P0(target_t) - vars.P0(target_t - 1);
        case 'P0_cum_energy'
            obj = sum(vars.P0(1:target_t)) * params.DT;
        otherwise
            error('Unknown objective_type: %s', objective_type);
    end

    ops = sdpsettings('solver', 'gurobi', 'verbose', params.YALMIP_VERBOSE, 'gurobi.OutputFlag', params.GUROBI_VERBOSE);
    soln = init_temporal_solution(objective_type, target_t, sense, sys, params);
    try
        if strcmpi(sense, 'max')
            sol = optimize(cons, -obj, ops);
        elseif strcmpi(sense, 'min')
            sol = optimize(cons, obj, ops);
        else
            error('Unknown sense: %s', sense);
        end
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

    obj_val = value(obj);
    if isnan(obj_val)
        soln.solver_status = 99;
        soln.solver_info = 'YALMIP returned NaN objective value';
        fprintf('  Optimization returned NaN objective; marking as failed.\n'); drawnow;
        return;
    end
    soln = collect_temporal_solution(soln, vars, obj_val);
    if soln.solver_status ~= 0
        fprintf('  Optimization value collection failed: problem=%d, info=%s\n', soln.solver_status, soln.solver_info); drawnow;
        return;
    end
    fprintf('  Success: objective_value=%.6f, P0=[%s]\n', soln.objective_value, sprintf(' %.4f', soln.P0)); drawnow;
end

function [cons, vars] = build_temporal_lindisflow_constraints(sys, ts, params)
    nb = sys.nb;
    nl = sys.nl;
    T = params.T;

    Pbr = sdpvar(nl, T, 'full');
    Qbr = sdpvar(nl, T, 'full');
    v = sdpvar(nb, T, 'full');
    PG = sdpvar(1, T, 'full');
    QG = sdpvar(1, T, 'full');
    PPV = sdpvar(1, T, 'full');
    QPV = sdpvar(1, T, 'full');
    PESS = sdpvar(1, T, 'full');
    QESS = sdpvar(1, T, 'full');
    E = sdpvar(1, T + 1, 'full');
    P0 = -Pbr(sys.root_branch, :);
    Q0 = -Qbr(sys.root_branch, :);

    cons = [];
    cons = [cons, E(1) == params.E0];
    cons = [cons, params.E_min <= E, E <= params.E_max];
    if params.ENFORCE_TERMINAL_SOC
        cons = [cons, E(T + 1) == params.E0];
    end

    for t = 1:T
        cons = [cons, v(sys.root_bus, t) == sys.V0^2]; %#ok<AGROW>
        cons = [cons, params.VMIN <= v(:, t), v(:, t) <= params.VMAX]; %#ok<AGROW>

        cons = [cons, params.PG_min <= PG(t), PG(t) <= params.PG_max]; %#ok<AGROW>
        cons = [cons, params.QG_min <= QG(t), QG(t) <= params.QG_max]; %#ok<AGROW>
        if t == 1
            prev_pg = params.PG0;
        else
            prev_pg = PG(t - 1);
        end
        cons = [cons, -params.RG_down <= PG(t) - prev_pg, PG(t) - prev_pg <= params.RG_up]; %#ok<AGROW>

        cons = [cons, 0 <= PPV(t), PPV(t) <= ts.PV_AVAIL_PROFILE(t)]; %#ok<AGROW>
        cons = [cons, -params.QPV_max <= QPV(t), QPV(t) <= params.QPV_max]; %#ok<AGROW>
        if params.USE_PV_SOC
            cons = [cons, cone([PPV(t); QPV(t)], params.SPV)]; %#ok<AGROW>
        end

        cons = [cons, -params.PESS_ch_max <= PESS(t), PESS(t) <= params.PESS_dis_max]; %#ok<AGROW>
        cons = [cons, -params.QESS_max <= QESS(t), QESS(t) <= params.QESS_max]; %#ok<AGROW>
        if params.USE_ESS_SOC_CONE
            cons = [cons, cone([PESS(t); QESS(t)], params.SESS)]; %#ok<AGROW>
        end
        cons = [cons, E(t + 1) == E(t) - PESS(t) * params.DT]; %#ok<AGROW>

        for e = 1:nl
            i = sys.from(e);
            j = sys.to(e);
            child_edges = sys.children_branches{j};
            [p_gen_j, q_gen_j, p_pv_j, q_pv_j, p_ess_j, q_ess_j] = device_injection_at_bus_time(j, sys, PG(t), QG(t), PPV(t), QPV(t), PESS(t), QESS(t));
            cons = [cons, Pbr(e, t) == sum(Pbr(child_edges, t)) + ts.Pload(j, t) - p_gen_j - p_pv_j - p_ess_j]; %#ok<AGROW>
            cons = [cons, Qbr(e, t) == sum(Qbr(child_edges, t)) + ts.Qload(j, t) - q_gen_j - q_pv_j - q_ess_j]; %#ok<AGROW>
            cons = [cons, v(j, t) == v(i, t) - 2 * (sys.r(e) * Pbr(e, t) / params.baseMVA + sys.x(e) * Qbr(e, t) / params.baseMVA)]; %#ok<AGROW>
            cons = [cons, -sys.Smax(e) <= Pbr(e, t), Pbr(e, t) <= sys.Smax(e)]; %#ok<AGROW>
            cons = [cons, -sys.Smax(e) <= Qbr(e, t), Qbr(e, t) <= sys.Smax(e)]; %#ok<AGROW>
            if params.USE_BRANCH_SOC
                cons = [cons, cone([Pbr(e, t); Qbr(e, t)], sys.Smax(e))]; %#ok<AGROW>
            end
        end
    end

    vars = struct('Pbr', Pbr, 'Qbr', Qbr, 'v', v, 'PG', PG, 'QG', QG, ...
        'PPV', PPV, 'QPV', QPV, 'PESS', PESS, 'QESS', QESS, 'E', E, 'P0', P0, 'Q0', Q0);
end

function [p_gen, q_gen, p_pv, q_pv, p_ess, q_ess] = device_injection_at_bus_time(bus_idx, sys, PG_t, QG_t, PPV_t, QPV_t, PESS_t, QESS_t)
    p_gen = 0; q_gen = 0; p_pv = 0; q_pv = 0; p_ess = 0; q_ess = 0;
    if bus_idx == sys.gen_bus
        p_gen = PG_t;
        q_gen = QG_t;
    end
    if bus_idx == sys.pv_bus
        p_pv = PPV_t;
        q_pv = QPV_t;
    end
    if bus_idx == sys.ess_bus
        p_ess = PESS_t;
        q_ess = QESS_t;
    end
end

function soln = init_temporal_solution(objective_type, target_t, sense, sys, params)
    soln = struct();
    soln.objective_type = objective_type;
    soln.sense = sense;
    soln.target_t = target_t;
    soln.objective_value = NaN;
    soln.solver_status = NaN;
    soln.solver_info = '';
    soln.P0 = NaN(1, params.T);
    soln.Q0 = NaN(1, params.T);
    soln.PG = NaN(1, params.T);
    soln.QG = NaN(1, params.T);
    soln.PPV = NaN(1, params.T);
    soln.QPV = NaN(1, params.T);
    soln.PESS = NaN(1, params.T);
    soln.QESS = NaN(1, params.T);
    soln.E = NaN(1, params.T + 1);
    soln.Pbr = NaN(sys.nl, params.T);
    soln.Qbr = NaN(sys.nl, params.T);
    soln.v = NaN(sys.nb, params.T);
end

function soln = collect_temporal_solution(soln, vars, objective_value)
    vals = [objective_value; value(vars.P0(:)); value(vars.Q0(:)); value(vars.PG(:)); value(vars.QG(:)); ...
            value(vars.PPV(:)); value(vars.QPV(:)); value(vars.PESS(:)); value(vars.QESS(:)); value(vars.E(:)); ...
            value(vars.Pbr(:)); value(vars.Qbr(:)); value(vars.v(:))];
    if any(isnan(vals))
        soln.solver_status = 98;
        soln.solver_info = 'YALMIP returned NaN dispatch values';
        fprintf('  Optimization dispatch values contain NaN; marking as failed.\n'); drawnow;
        return;
    end
    soln.objective_value = objective_value;
    soln.P0 = value(vars.P0);
    soln.Q0 = value(vars.Q0);
    soln.PG = value(vars.PG);
    soln.QG = value(vars.QG);
    soln.PPV = value(vars.PPV);
    soln.QPV = value(vars.QPV);
    soln.PESS = value(vars.PESS);
    soln.QESS = value(vars.QESS);
    soln.E = value(vars.E);
    soln.Pbr = value(vars.Pbr);
    soln.Qbr = value(vars.Qbr);
    soln.v = value(vars.v);
end

function temporal_bounds = assemble_temporal_bounds_table(all_solutions, T)
    temporal_bounds = repmat(struct('t', NaN, 'P0_min', NaN, 'P0_max', NaN, ...
        'P0_min_status', NaN, 'P0_max_status', NaN, 'R0_min', NaN, 'R0_max', NaN, ...
        'R0_min_status', NaN, 'R0_max_status', NaN, 'P0_cum_energy_min', NaN, ...
        'P0_cum_energy_max', NaN, 'P0_cum_energy_min_status', NaN, 'P0_cum_energy_max_status', NaN), T, 1);
    for t = 1:T
        temporal_bounds(t).t = t;
    end
    for i = 1:numel(all_solutions)
        s = all_solutions(i);
        t = s.target_t;
        val = s.objective_value;
        status = s.solver_status;
        switch s.objective_type
            case 'P0_power'
                if strcmpi(s.sense, 'min')
                    temporal_bounds(t).P0_min = val;
                    temporal_bounds(t).P0_min_status = status;
                else
                    temporal_bounds(t).P0_max = val;
                    temporal_bounds(t).P0_max_status = status;
                end
            case 'P0_ramp'
                if strcmpi(s.sense, 'min')
                    temporal_bounds(t).R0_min = val;
                    temporal_bounds(t).R0_min_status = status;
                else
                    temporal_bounds(t).R0_max = val;
                    temporal_bounds(t).R0_max_status = status;
                end
            case 'P0_cum_energy'
                if strcmpi(s.sense, 'min')
                    temporal_bounds(t).P0_cum_energy_min = val;
                    temporal_bounds(t).P0_cum_energy_min_status = status;
                else
                    temporal_bounds(t).P0_cum_energy_max = val;
                    temporal_bounds(t).P0_cum_energy_max_status = status;
                end
        end
    end
end

function save_temporal_results_csv(all_solutions, temporal_bounds, result_dir, params)
    writetable(struct2table(temporal_bounds), fullfile(result_dir, 'temporal_p0_bounds.csv'));
    fprintf('Saved CSV: %s\n', fullfile(result_dir, 'temporal_p0_bounds.csv')); drawnow;

    rows = repmat(empty_solution_csv_row(params), numel(all_solutions), 1);
    for i = 1:numel(all_solutions)
        rows(i) = solution_to_csv_row(all_solutions(i), params);
    end
    writetable(struct2table(rows), fullfile(result_dir, 'temporal_boundary_solutions.csv'));
    fprintf('Saved CSV: %s\n', fullfile(result_dir, 'temporal_boundary_solutions.csv')); drawnow;
end

function row = empty_solution_csv_row(params)
    row = struct();
    row.objective_type = '';
    row.sense = '';
    row.target_t = NaN;
    row.objective_value = NaN;
    row.solver_status = NaN;
    row.solver_info = '';
    for t = 1:params.T
        row.(sprintf('P0_%d', t)) = NaN;
    end
    for t = 1:params.T
        row.(sprintf('PG_%d', t)) = NaN;
    end
    for t = 1:params.T
        row.(sprintf('PPV_%d', t)) = NaN;
    end
    for t = 1:params.T
        row.(sprintf('PESS_%d', t)) = NaN;
    end
    for t = 1:(params.T + 1)
        row.(sprintf('E_%d', t)) = NaN;
    end
end

function row = solution_to_csv_row(soln, params)
    row = empty_solution_csv_row(params);
    row.objective_type = soln.objective_type;
    row.sense = soln.sense;
    row.target_t = soln.target_t;
    row.objective_value = soln.objective_value;
    row.solver_status = soln.solver_status;
    row.solver_info = soln.solver_info;
    for t = 1:params.T
        row.(sprintf('P0_%d', t)) = soln.P0(t);
        row.(sprintf('PG_%d', t)) = soln.PG(t);
        row.(sprintf('PPV_%d', t)) = soln.PPV(t);
        row.(sprintf('PESS_%d', t)) = soln.PESS(t);
    end
    for t = 1:(params.T + 1)
        row.(sprintf('E_%d', t)) = soln.E(t);
    end
end

function plot_p0_power_bounds(temporal_bounds, params)
    tbl = struct2table(temporal_bounds);
    t = tbl.t;
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 900 520]);
    ax = axes(fig); %#ok<LAXES>
    hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
    hBand = plot_step_band(ax, t, tbl.P0_min, tbl.P0_max, [0.84 0.92 1.00], 'P0 feasible band');
    [xs_min, ys_min] = make_step_xy(t, tbl.P0_min, 'post');
    [xs_max, ys_max] = make_step_xy(t, tbl.P0_max, 'post');
    hMin = plot(ax, xs_min, ys_min, 'b-', 'LineWidth', 1.8, 'DisplayName', 'P0 min');
    hMax = plot(ax, xs_max, ys_max, 'r-', 'LineWidth', 1.8, 'DisplayName', 'P0 max');
    yline_compat(ax, 0, 'k--', 'P_0=0');
    xlabel(ax, 't');
    ylabel(ax, 'P_0 (MW)');
    title(ax, 'IEEE-13 simplified feeder | temporal P0 power bounds');
    xlim(ax, [min(t), max(t) + 1]);
    set(ax, 'XTick', min(t):(max(t) + 1));
    legend(ax, [hBand, hMin, hMax], {'P0 feasible band', 'P0 min', 'P0 max'}, 'Location', 'best');
    save_png_compat(fig, fullfile(params.RESULT_DIR, 'p0_power_bounds.png'));
    close(fig);
    fprintf('Saved figure: %s\n', fullfile(params.RESULT_DIR, 'p0_power_bounds.png')); drawnow;
end

function plot_p0_ramp_bounds(temporal_bounds, params)
    tbl = struct2table(temporal_bounds);
    mask = ~isnan(tbl.R0_min) | ~isnan(tbl.R0_max);
    t = tbl.t(mask);
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 900 520]);
    ax = axes(fig); %#ok<LAXES>
    hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
    handles = gobjects(0);
    labels = {};
    if ~isempty(t)
        hBand = plot_step_band(ax, t, tbl.R0_min(mask), tbl.R0_max(mask), [0.90 0.95 0.88], 'R0 feasible band');
        [xs_min, ys_min] = make_step_xy(t, tbl.R0_min(mask), 'post');
        [xs_max, ys_max] = make_step_xy(t, tbl.R0_max(mask), 'post');
        hMin = plot(ax, xs_min, ys_min, 'b-', 'LineWidth', 1.8, 'DisplayName', 'R0 min');
        hMax = plot(ax, xs_max, ys_max, 'r-', 'LineWidth', 1.8, 'DisplayName', 'R0 max');
        handles = [hBand, hMin, hMax]; %#ok<AGROW>
        labels = {'R0 feasible band', 'R0 min', 'R0 max'};
        xlim(ax, [min(t), max(t) + 1]);
        set(ax, 'XTick', min(t):(max(t) + 1));
    end
    yline_compat(ax, 0, 'k--', '0');
    xlabel(ax, 't');
    ylabel(ax, 'P_0(t)-P_0(t-1) (MW/step)');
    title(ax, 'IEEE-13 simplified feeder | temporal P0 ramp bounds');
    if ~isempty(handles)
        legend(ax, handles, labels, 'Location', 'best');
    end
    save_png_compat(fig, fullfile(params.RESULT_DIR, 'p0_ramp_bounds.png'));
    close(fig);
    fprintf('Saved figure: %s\n', fullfile(params.RESULT_DIR, 'p0_ramp_bounds.png')); drawnow;
end

function plot_p0_cum_energy_bounds(temporal_bounds, params)
    tbl = struct2table(temporal_bounds);
    t = tbl.t;
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 900 520]);
    ax = axes(fig); %#ok<LAXES>
    hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
    hBand = plot_step_band(ax, t, tbl.P0_cum_energy_min, tbl.P0_cum_energy_max, [0.95 0.90 1.00], 'cumulative exchange band');
    [xs_min, ys_min] = make_step_xy(t, tbl.P0_cum_energy_min, 'post');
    [xs_max, ys_max] = make_step_xy(t, tbl.P0_cum_energy_max, 'post');
    hMin = plot(ax, xs_min, ys_min, 'b-', 'LineWidth', 1.8, 'DisplayName', 'P0 cumulative min');
    hMax = plot(ax, xs_max, ys_max, 'r-', 'LineWidth', 1.8, 'DisplayName', 'P0 cumulative max');
    yline_compat(ax, 0, 'k--', '0');
    xlabel(ax, 't');
    ylabel(ax, 'P0 cumulative exchange energy (MWh)');
    title(ax, 'IEEE-13 simplified feeder | P0 cumulative exchange energy bounds');
    xlim(ax, [min(t), max(t) + 1]);
    set(ax, 'XTick', min(t):(max(t) + 1));
    legend(ax, [hBand, hMin, hMax], {'cumulative exchange band', 'P0 cumulative min', 'P0 cumulative max'}, 'Location', 'best');
    save_png_compat(fig, fullfile(params.RESULT_DIR, 'p0_cum_energy_bounds.png'));
    close(fig);
    fprintf('Saved figure: %s\n', fullfile(params.RESULT_DIR, 'p0_cum_energy_bounds.png')); drawnow;
end

function plot_ess_soc_examples(all_solutions, params)
    examples = { ...
        'P0_power', 'max', 1, 'P0 power max at t=1'; ...
        'P0_power', 'min', 1, 'P0 power min at t=1'; ...
        'P0_cum_energy', 'max', params.T, sprintf('P0 cumulative max at t=%d', params.T); ...
        'P0_cum_energy', 'min', params.T, sprintf('P0 cumulative min at t=%d', params.T)};
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 900 520]);
    ax = axes(fig); %#ok<LAXES>
    hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
    handles = gobjects(0);
    labels = {};
    for i = 1:size(examples, 1)
        soln = find_solution(all_solutions, examples{i, 1}, examples{i, 2}, examples{i, 3});
        if isempty(soln) || soln.solver_status ~= 0 || any(isnan(soln.E))
            fprintf('Skipping ESS SOC example: %s, sense=%s, target_t=%d (missing or failed).\n', examples{i,1}, examples{i,2}, examples{i,3}); drawnow;
            continue;
        end
        h = stairs(ax, 0:params.T, soln.E, 'LineWidth', 1.7, 'DisplayName', examples{i, 4});
        handles(end+1) = h; %#ok<AGROW>
        labels{end+1} = examples{i, 4}; %#ok<AGROW>
    end
    yline_compat(ax, params.E_min, 'k--', 'E_{min}');
    yline_compat(ax, params.E_max, 'k--', 'E_{max}');
    xl = xlim(ax);
    text(ax, xl(1), params.E_min, '  E_{min}', 'VerticalAlignment', 'bottom', 'Color', [0.1 0.1 0.1]);
    text(ax, xl(1), params.E_max, '  E_{max}', 'VerticalAlignment', 'top', 'Color', [0.1 0.1 0.1]);
    xlabel(ax, 'time index (0 initial, 1..T period ends)');
    ylabel(ax, 'ESS SOC / stored energy E (MWh)');
    title(ax, 'IEEE-13 simplified feeder | ESS SOC trajectories from temporal boundary examples');
    xlim(ax, [0, params.T]);
    set(ax, 'XTick', 0:params.T);
    if ~isempty(handles)
        legend(ax, handles, labels, 'Location', 'best');
    end
    save_png_compat(fig, fullfile(params.RESULT_DIR, 'ess_soc_boundary_examples.png'));
    close(fig);
    fprintf('Saved figure: %s\n', fullfile(params.RESULT_DIR, 'ess_soc_boundary_examples.png')); drawnow;
end

function plot_all_boundary_p0_trajectories(all_solutions, temporal_bounds, params)
    tbl = struct2table(temporal_bounds);
    t = tbl.t;
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 980 560]);
    ax = axes(fig); %#ok<LAXES>
    hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');

    hBand = plot_step_band(ax, t, tbl.P0_min, tbl.P0_max, [0.84 0.92 1.00], 'P0 feasible band');
    hAll = gobjects(0);
    for i = 1:numel(all_solutions)
        soln = all_solutions(i);
        if soln.solver_status ~= 0 || any(isnan(soln.P0))
            continue;
        end
        [xs, ys] = make_step_xy(t, soln.P0, 'post');
        h = plot(ax, xs, ys, '-', 'Color', [0.72 0.72 0.72], 'LineWidth', 0.75, 'HandleVisibility', 'off');
        if isempty(hAll)
            hAll = h;
        end
    end
    if ~isempty(hAll)
        set(hAll, 'DisplayName', 'all boundary trajectories', 'HandleVisibility', 'on');
    end

    highlight_specs = { ...
        'P0_power', 'max', 1, 'P0 power max at t=1', [0.85 0.10 0.10]; ...
        'P0_power', 'min', 1, 'P0 power min at t=1', [0.00 0.25 0.85]; ...
        'P0_power', 'max', params.T, sprintf('P0 power max at t=%d', params.T), [0.95 0.45 0.05]; ...
        'P0_power', 'min', params.T, sprintf('P0 power min at t=%d', params.T), [0.00 0.55 0.75]; ...
        'P0_cum_energy', 'max', params.T, sprintf('P0 cumulative max at t=%d', params.T), [0.45 0.10 0.75]; ...
        'P0_cum_energy', 'min', params.T, sprintf('P0 cumulative min at t=%d', params.T), [0.10 0.55 0.20]};
    hHighlights = gobjects(0);
    highlight_labels = {};
    for i = 1:size(highlight_specs, 1)
        soln = find_solution(all_solutions, highlight_specs{i, 1}, highlight_specs{i, 2}, highlight_specs{i, 3});
        if isempty(soln) || soln.solver_status ~= 0 || any(isnan(soln.P0))
            fprintf('Skipping highlighted P0 trajectory: %s, sense=%s, target_t=%d (missing or failed).\n', highlight_specs{i,1}, highlight_specs{i,2}, highlight_specs{i,3}); drawnow;
            continue;
        end
        [xs, ys] = make_step_xy(t, soln.P0, 'post');
        h = plot(ax, xs, ys, '-', 'Color', highlight_specs{i, 5}, 'LineWidth', 2.0, 'DisplayName', highlight_specs{i, 4});
        hHighlights(end+1) = h; %#ok<AGROW>
        highlight_labels{end+1} = highlight_specs{i, 4}; %#ok<AGROW>
    end

    yline_compat(ax, 0, 'k--', 'P_0=0');
    xlabel(ax, 't');
    ylabel(ax, 'P0 (MW)');
    title(ax, 'IEEE-13 simplified feeder | all temporal boundary P0 trajectories');
    xlim(ax, [min(t), max(t) + 1]);
    set(ax, 'XTick', min(t):(max(t) + 1));
    if isempty(hAll)
        legend(ax, [hBand, hHighlights], [{'P0 feasible band'}, highlight_labels], 'Location', 'best');
    else
        legend(ax, [hBand, hAll, hHighlights], [{'P0 feasible band', 'all boundary trajectories'}, highlight_labels], 'Location', 'best');
    end
    save_png_compat(fig, fullfile(params.RESULT_DIR, 'p0_all_boundary_trajectories.png'));
    close(fig);
    fprintf('Saved figure: %s\n', fullfile(params.RESULT_DIR, 'p0_all_boundary_trajectories.png')); drawnow;
end

function soln = find_solution(all_solutions, objective_type, sense, target_t)
    soln = [];
    for i = 1:numel(all_solutions)
        if strcmpi(all_solutions(i).objective_type, objective_type) && strcmpi(all_solutions(i).sense, sense) && all_solutions(i).target_t == target_t
            soln = all_solutions(i);
            return;
        end
    end
end

function [x_step, y_step] = make_step_xy(x, y, side)
    if nargin < 3
        side = 'post';
    end
    x = x(:).';
    y = y(:).';
    if numel(x) ~= numel(y)
        error('make_step_xy requires x and y to have the same length.');
    end
    if isempty(x)
        x_step = [];
        y_step = [];
        return;
    end
    if ~strcmpi(side, 'post')
        error('Only post step mode is implemented in this script.');
    end
    if numel(x) == 1
        dx = 1;
    else
        dx = median(diff(x));
    end
    x_edges = [x, x(end) + dx];
    x_step = zeros(1, 2 * numel(y));
    y_step = zeros(1, 2 * numel(y));
    for k = 1:numel(y)
        x_step(2*k - 1) = x_edges(k);
        x_step(2*k) = x_edges(k + 1);
        y_step(2*k - 1) = y(k);
        y_step(2*k) = y(k);
    end
end

function [x_poly, y_poly] = make_step_band_polygon(x, y_lower, y_upper)
    [x_low, y_low] = make_step_xy(x, y_lower, 'post');
    [x_up, y_up] = make_step_xy(x, y_upper, 'post');
    x_poly = [x_low, fliplr(x_up)];
    y_poly = [y_low, fliplr(y_up)];
end

function hBand = plot_step_band(ax, x, y_lower, y_upper, band_color, display_name)
    [x_poly, y_poly] = make_step_band_polygon(x, y_lower, y_upper);
    hBand = fill(ax, x_poly, y_poly, band_color, 'EdgeColor', 'none', 'FaceAlpha', 0.70, 'DisplayName', display_name);
end

function yline_compat(ax, y, style, label_text)
    try
        yline(ax, y, style, label_text, 'HandleVisibility', 'off');
    catch
        xl = xlim(ax);
        plot(ax, xl, [y y], style, 'HandleVisibility', 'off');
    end
end


function all_solutions = annotate_vgvb_decomposition(all_solutions, ts, params)
    max_residual = 0;
    for i = 1:numel(all_solutions)
        if all_solutions(i).solver_status ~= 0 || any(isnan(all_solutions(i).P0))
            all_solutions(i).PVG = NaN(1, params.T);
            all_solutions(i).QVG = NaN(1, params.T);
            all_solutions(i).PVB = NaN(1, params.T);
            all_solutions(i).QVB = NaN(1, params.T);
            all_solutions(i).P0_decomp_residual = NaN(1, params.T);
            continue;
        end
        totalP = sum(ts.Pload, 1);
        totalQ = sum(ts.Qload, 1);
        all_solutions(i).PVG = all_solutions(i).PG + all_solutions(i).PPV - totalP;
        all_solutions(i).QVG = all_solutions(i).QG + all_solutions(i).QPV - totalQ;
        all_solutions(i).PVB = all_solutions(i).PESS;
        all_solutions(i).QVB = all_solutions(i).QESS;
        all_solutions(i).P0_decomp_residual = all_solutions(i).P0 - all_solutions(i).PVG - all_solutions(i).PVB;
        max_residual = max(max_residual, max(abs(all_solutions(i).P0_decomp_residual)));
    end
    fprintf('P0 VG/VB decomposition max residual = %.6e\n', max_residual); drawnow;
    if max_residual > 1e-6
        warning('P0 decomposition residual is larger than 1e-6. Check sign conventions and losses.');
    end
end

function netmat = build_network_matrices_lindisflow(sys, ts, params)
    nb = sys.nb; nl = sys.nl; T = params.T;
    D = zeros(nl, nb);
    H = zeros(nb, nl);
    for e = 1:nl
        D(e, :) = downstream_buses_for_branch(e, sys);
    end
    for b = 1:nb
        H(b, :) = path_branches_to_bus(b, sys);
    end
    gen_down = D(:, sys.gen_bus);
    pv_down = D(:, sys.pv_bus);
    ess_down = D(:, sys.ess_bus);
    PcoefG = [-gen_down, -pv_down, zeros(nl,1), zeros(nl,1)];
    PcoefB = [-ess_down, zeros(nl,1)];
    QcoefG = [zeros(nl,1), zeros(nl,1), -gen_down, -pv_down];
    QcoefB = [zeros(nl,1), -ess_down];
    rrow = (sys.r(:)' / params.baseMVA);
    xrow = (sys.x(:)' / params.baseMVA);
    VcoefG = -2 * (H * (diag(rrow) * PcoefG + diag(xrow) * QcoefG));
    VcoefB = -2 * (H * (diag(rrow) * PcoefB + diag(xrow) * QcoefB));
    nrows = 4*nl + 2*nb;
    AG = zeros(nrows, 4); AB = zeros(nrows, 2); row_type = cell(nrows,1); row_name = cell(nrows,1); side = cell(nrows,1);
    row = 0;
    for e = 1:nl
        row=row+1; AG(row,:)=PcoefG(e,:); AB(row,:)=PcoefB(e,:); row_type{row}='branch_P_upper'; row_name{row}=sprintf('%s P upper', sys.branch_name{e}); side{row}='upper';
    end
    for e = 1:nl
        row=row+1; AG(row,:)=-PcoefG(e,:); AB(row,:)=-PcoefB(e,:); row_type{row}='branch_P_lower'; row_name{row}=sprintf('%s P lower', sys.branch_name{e}); side{row}='lower';
    end
    for e = 1:nl
        row=row+1; AG(row,:)=QcoefG(e,:); AB(row,:)=QcoefB(e,:); row_type{row}='branch_Q_upper'; row_name{row}=sprintf('%s Q upper', sys.branch_name{e}); side{row}='upper';
    end
    for e = 1:nl
        row=row+1; AG(row,:)=-QcoefG(e,:); AB(row,:)=-QcoefB(e,:); row_type{row}='branch_Q_lower'; row_name{row}=sprintf('%s Q lower', sys.branch_name{e}); side{row}='lower';
    end
    for b = 1:nb
        row=row+1; AG(row,:)=VcoefG(b,:); AB(row,:)=VcoefB(b,:); row_type{row}='voltage_upper'; row_name{row}=sprintf('%s V upper', sys.bus_name{b}); side{row}='upper';
    end
    for b = 1:nb
        row=row+1; AG(row,:)=-VcoefG(b,:); AB(row,:)=-VcoefB(b,:); row_type{row}='voltage_lower'; row_name{row}=sprintf('%s V lower', sys.bus_name{b}); side{row}='lower';
    end
    b_by_t = zeros(nrows, T);
    for t = 1:T
        Pload_down = D * ts.Pload(:,t);
        Qload_down = D * ts.Qload(:,t);
        Vconst = sys.V0^2 - 2 * (H * (sys.r .* Pload_down / params.baseMVA + sys.x .* Qload_down / params.baseMVA));
        row = 0;
        for e=1:nl, row=row+1; b_by_t(row,t)=sys.Smax(e)-Pload_down(e); end
        for e=1:nl, row=row+1; b_by_t(row,t)=sys.Smax(e)+Pload_down(e); end
        for e=1:nl, row=row+1; b_by_t(row,t)=sys.Smax(e)-Qload_down(e); end
        for e=1:nl, row=row+1; b_by_t(row,t)=sys.Smax(e)+Qload_down(e); end
        for b=1:nb, row=row+1; b_by_t(row,t)=params.VMAX-Vconst(b); end
        for b=1:nb, row=row+1; b_by_t(row,t)=Vconst(b)-params.VMIN; end
    end
    netmat = struct('AG',AG,'AB',AB,'b_by_t',b_by_t,'row_type',{row_type},'row_name',{row_name}, ...
        'constraint_side',{side},'downstream_matrix',D,'path_matrix',H,'varG_names',{{'PG','PPV','QG','QPV'}}, ...
        'varB_names',{{'PESS','QESS'}},'PcoefG',PcoefG,'PcoefB',PcoefB,'QcoefG',QcoefG,'QcoefB',QcoefB,'VcoefG',VcoefG,'VcoefB',VcoefB);
    fprintf('Network matrix rows: %d\n', nrows); drawnow;
    print_row_type_counts(row_type);
end

function down = downstream_buses_for_branch(e, sys)
    down = false(1, sys.nb);
    start_bus = sys.to(e);
    queue = start_bus;
    down(start_bus) = true;
    while ~isempty(queue)
        b = queue(1); queue(1) = [];
        child_edges = sys.children_branches{b};
        for k = 1:numel(child_edges)
            child_bus = sys.to(child_edges(k));
            if ~down(child_bus)
                down(child_bus) = true;
                queue(end+1) = child_bus; %#ok<AGROW>
            end
        end
    end
    down = double(down);
end

function path = path_branches_to_bus(bus_idx, sys)
    path = zeros(1, sys.nl);
    current = bus_idx;
    while current ~= sys.root_bus
        e = find(sys.to == current, 1);
        if isempty(e), error('No parent branch found for bus %d.', current); end
        path(e) = 1;
        current = sys.from(e);
    end
end

function print_row_type_counts(row_type)
    types = unique(row_type);
    for i = 1:numel(types)
        fprintf('  %s rows: %d\n', types{i}, sum(strcmp(row_type, types{i}))); drawnow;
    end
end

function validation = validate_network_matrix_against_original(all_solutions, netmat, sys, ts, params)
    max_diff = 0; max_viol_mat = 0; max_viol_orig = 0; nsol = 0; nrows = 0;
    for i = 1:numel(all_solutions)
        s = all_solutions(i);
        if s.solver_status ~= 0 || any(isnan(s.P0)), continue; end
        nsol = nsol + 1;
        for t = 1:params.T
            xG = [s.PG(t); s.PPV(t); s.QG(t); s.QPV(t)];
            xB = [s.PESS(t); s.QESS(t)];
            slack_mat = netmat.b_by_t(:,t) - netmat.AG*xG - netmat.AB*xB;
            slack_orig = original_network_slack_vector(s, t, sys, params);
            max_diff = max(max_diff, max(abs(slack_mat - slack_orig)));
            max_viol_mat = max(max_viol_mat, max(max(-slack_mat, 0)));
            max_viol_orig = max(max_viol_orig, max(max(-slack_orig, 0)));
            nrows = nrows + numel(slack_mat);
        end
    end
    validation = struct('max_abs_slack_diff', max_diff, 'max_network_violation_matrix', max_viol_mat, ...
        'max_network_violation_original', max_viol_orig, 'num_checked_solutions', nsol, ...
        'num_checked_rows', nrows, 'is_passed', max_diff < 1e-6 && max_viol_mat < 1e-6);
    if validation.is_passed
        fprintf('Network matrix validation: PASS\n'); drawnow;
    else
        fprintf('Network matrix validation: FAIL\n'); drawnow;
    end
    fprintf('max_abs_slack_diff = %.6e\n', max_diff); drawnow;
    fprintf('max_network_violation_matrix = %.6e\n', max_viol_mat); drawnow;
    fprintf('max_network_violation_original = %.6e\n', max_viol_orig); drawnow;
end

function slack = original_network_slack_vector(s, t, sys, params)
    slack = [sys.Smax - s.Pbr(:,t); sys.Smax + s.Pbr(:,t); sys.Smax - s.Qbr(:,t); sys.Smax + s.Qbr(:,t); params.VMAX - s.v(:,t); s.v(:,t) - params.VMIN];
end

function save_network_matrix_validation_csv(validation, result_dir)
    writetable(struct2table(validation), fullfile(result_dir, 'network_matrix_validation.csv'));
    fprintf('Saved CSV: %s\n', fullfile(result_dir, 'network_matrix_validation.csv')); drawnow;
end

function decoup = compute_robust_decoupled_limits(netmat, ts, params)
    T = params.T; nrow = size(netmat.AG,1);
    xBmin = [-params.PESS_ch_max; -params.QESS_max];
    xBmax = [ params.PESS_dis_max;  params.QESS_max];
    deltaB = worst_case_box(netmat.AB, xBmin, xBmax);
    deltaG_by_t = zeros(nrow, T);
    for t=1:T
        xGmin = [params.PG_min; 0; params.QG_min; -params.QPV_max];
        xGmax = [params.PG_max; ts.PV_AVAIL_PROFILE(t); params.QG_max; params.QPV_max];
        deltaG_by_t(:,t) = worst_case_box(netmat.AG, xGmin, xGmax);
    end
    nG_by_t = netmat.b_by_t - repmat(deltaB, 1, T);
    nB_by_t = netmat.b_by_t - deltaG_by_t;
    decoup = struct('NG', netmat.AG, 'nG_by_t', nG_by_t, 'NB', netmat.AB, 'nB_by_t', nB_by_t, ...
        'deltaB', deltaB, 'deltaG_by_t', deltaG_by_t, 'row_type', {netmat.row_type}, 'row_name', {netmat.row_name});
    fprintf('VG robust margin deltaB: min=%.6f max=%.6f mean=%.6f\n', min(deltaB), max(deltaB), mean(deltaB)); drawnow;
    fprintf('VB robust margin deltaG: min=%.6f max=%.6f mean=%.6f\n', min(deltaG_by_t(:)), max(deltaG_by_t(:)), mean(deltaG_by_t(:))); drawnow;
    fprintf('Negative nG rows over all t: %d\n', sum(nG_by_t(:) < -1e-9)); drawnow;
    fprintf('Negative nB rows over all t: %d\n', sum(nB_by_t(:) < -1e-9)); drawnow;
    if any(nG_by_t(:) < -1e-9) || any(nB_by_t(:) < -1e-9)
        warning('Some robust decoupled network RHS values are negative; decoupling may be very conservative or infeasible.');
    end
end

function val = worst_case_box(A, xmin, xmax)
    val = zeros(size(A,1),1);
    for r=1:size(A,1)
        for k=1:size(A,2)
            if A(r,k) >= 0
                val(r) = val(r) + A(r,k)*xmax(k);
            else
                val(r) = val(r) + A(r,k)*xmin(k);
            end
        end
    end
end

function save_decoupled_network_margins_csv(decoup, netmat, result_dir)
    rows = struct([]);
    for t=1:size(netmat.b_by_t,2)
        for r=1:size(netmat.b_by_t,1)
            row = struct('t',t,'row_id',r,'row_type',netmat.row_type{r},'row_name',netmat.row_name{r}, ...
                'b_original',netmat.b_by_t(r,t),'deltaB_for_VG',decoup.deltaB(r),'nG',decoup.nG_by_t(r,t), ...
                'deltaG_for_VB',decoup.deltaG_by_t(r,t),'nB',decoup.nB_by_t(r,t));
            rows = append_struct(rows, row);
        end
    end
    writetable(struct2table(rows), fullfile(result_dir, 'decoupled_network_margins.csv'));
    fprintf('Saved CSV: %s\n', fullfile(result_dir, 'decoupled_network_margins.csv')); drawnow;
end

function sols = solve_all_coupled_vgvb_projection_bounds(sys, ts, params)
    sols = struct([]);
    for t=1:params.T
        for is=1:2
            sense=pick_sense(is); fprintf('Solving coupled VG/VB projection: type=VG_power, sense=%s, target_t=%d/%d ...\n', sense, t, params.T); drawnow;
            sols=append_struct(sols, solve_coupled_vgvb_projection_boundary('VG_power',t,sense,sys,ts,params));
        end
    end
    for t=2:params.T
        for is=1:2
            sense=pick_sense(is); fprintf('Solving coupled VG/VB projection: type=VG_ramp, sense=%s, target_t=%d/%d ...\n', sense, t, params.T); drawnow;
            sols=append_struct(sols, solve_coupled_vgvb_projection_boundary('VG_ramp',t,sense,sys,ts,params));
        end
    end
    for t=1:params.T
        for is=1:2
            sense=pick_sense(is); fprintf('Solving coupled VG/VB projection: type=VB_power, sense=%s, target_t=%d/%d ...\n', sense, t, params.T); drawnow;
            sols=append_struct(sols, solve_coupled_vgvb_projection_boundary('VB_power',t,sense,sys,ts,params));
        end
    end
    for t=1:params.T
        for is=1:2
            sense=pick_sense(is); fprintf('Solving coupled VG/VB projection: type=VB_energy, sense=%s, target_t=%d/%d ...\n', sense, t, params.T); drawnow;
            sols=append_struct(sols, solve_coupled_vgvb_projection_boundary('VB_energy',t,sense,sys,ts,params));
        end
    end
end

function soln = solve_coupled_vgvb_projection_boundary(objective_type,target_t,sense,sys,ts,params)
    [cons, vars] = build_temporal_lindisflow_constraints(sys, ts, params);
    totalP = sum(ts.Pload,1);
    PVG = vars.PG + vars.PPV - totalP;
    PVB = vars.PESS;
    switch objective_type
        case 'VG_power', obj = PVG(target_t); rtype='VG';
        case 'VG_ramp', obj = PVG(target_t)-PVG(target_t-1); rtype='VG';
        case 'VB_power', obj = PVB(target_t); rtype='VB';
        case 'VB_energy', obj = vars.E(target_t+1); rtype='VB';
        otherwise, error('Unknown objective_type %s', objective_type);
    end
    soln = init_vgvb_solution(rtype, objective_type, target_t, sense, params);
    soln = optimize_and_collect_vgvb(soln, cons, obj, sense, vars, ts, params);
end

function sols = solve_all_decoupled_vgvb_bounds(decoup, sys, ts, params)
    sols = struct([]);
    for t=1:params.T
        for is=1:2
            sense=pick_sense(is); fprintf('Solving decoupled VG robust boundary: type=VG_power, sense=%s, target_t=%d/%d ...\n', sense, t, params.T); drawnow;
            sols=append_struct(sols, solve_vg_decoupled_boundary('VG_power',t,sense,decoup,sys,ts,params));
        end
    end
    for t=2:params.T
        for is=1:2
            sense=pick_sense(is); fprintf('Solving decoupled VG robust boundary: type=VG_ramp, sense=%s, target_t=%d/%d ...\n', sense, t, params.T); drawnow;
            sols=append_struct(sols, solve_vg_decoupled_boundary('VG_ramp',t,sense,decoup,sys,ts,params));
        end
    end
    for t=1:params.T
        for is=1:2
            sense=pick_sense(is); fprintf('Solving decoupled VB robust boundary: type=VB_power, sense=%s, target_t=%d/%d ...\n', sense, t, params.T); drawnow;
            sols=append_struct(sols, solve_vb_decoupled_boundary('VB_power',t,sense,decoup,sys,ts,params));
        end
    end
    for t=1:params.T
        for is=1:2
            sense=pick_sense(is); fprintf('Solving decoupled VB robust boundary: type=VB_energy, sense=%s, target_t=%d/%d ...\n', sense, t, params.T); drawnow;
            sols=append_struct(sols, solve_vb_decoupled_boundary('VB_energy',t,sense,decoup,sys,ts,params));
        end
    end
end

function soln = solve_vg_decoupled_boundary(objective_type,target_t,sense,decoup,sys,ts,params)
    T=params.T; PG=sdpvar(1,T,'full'); QG=sdpvar(1,T,'full'); PPV=sdpvar(1,T,'full'); QPV=sdpvar(1,T,'full');
    cons=[];
    for t=1:T
        cons=[cons, params.PG_min<=PG(t), PG(t)<=params.PG_max, params.QG_min<=QG(t), QG(t)<=params.QG_max]; %#ok<AGROW>
        cons=[cons, 0<=PPV(t), PPV(t)<=ts.PV_AVAIL_PROFILE(t), -params.QPV_max<=QPV(t), QPV(t)<=params.QPV_max]; %#ok<AGROW>
        if t==1, prev=params.PG0; else, prev=PG(t-1); end
        cons=[cons, -params.RG_down<=PG(t)-prev, PG(t)-prev<=params.RG_up]; %#ok<AGROW>
        xG=[PG(t);PPV(t);QG(t);QPV(t)];
        cons=[cons, decoup.NG*xG <= decoup.nG_by_t(:,t)]; %#ok<AGROW>
    end
    PVG=PG+PPV-sum(ts.Pload,1);
    switch objective_type
        case 'VG_power', obj=PVG(target_t);
        case 'VG_ramp', obj=PVG(target_t)-PVG(target_t-1);
        otherwise, error('Unknown VG objective %s', objective_type);
    end
    vars=struct('PG',PG,'QG',QG,'PPV',PPV,'QPV',QPV,'PESS',NaN(1,T),'QESS',NaN(1,T),'E',NaN(1,T+1));
    soln=init_vgvb_solution('VG',objective_type,target_t,sense,params);
    soln=optimize_and_collect_vgvb(soln,cons,obj,sense,vars,ts,params);
end

function soln = solve_vb_decoupled_boundary(objective_type,target_t,sense,decoup,sys,ts,params) %#ok<INUSD>
    T=params.T; PESS=sdpvar(1,T,'full'); QESS=sdpvar(1,T,'full'); E=sdpvar(1,T+1,'full');
    cons=[E(1)==params.E0, params.E_min<=E, E<=params.E_max];
    if params.ENFORCE_TERMINAL_SOC, cons=[cons, E(T+1)==params.E0]; end
    for t=1:T
        cons=[cons, -params.PESS_ch_max<=PESS(t), PESS(t)<=params.PESS_dis_max, -params.QESS_max<=QESS(t), QESS(t)<=params.QESS_max]; %#ok<AGROW>
        cons=[cons, E(t+1)==E(t)-PESS(t)*params.DT]; %#ok<AGROW>
        xB=[PESS(t);QESS(t)]; cons=[cons, decoup.NB*xB <= decoup.nB_by_t(:,t)]; %#ok<AGROW>
    end
    switch objective_type
        case 'VB_power', obj=PESS(target_t);
        case 'VB_energy', obj=E(target_t+1);
        otherwise, error('Unknown VB objective %s', objective_type);
    end
    vars=struct('PG',NaN(1,T),'QG',NaN(1,T),'PPV',NaN(1,T),'QPV',NaN(1,T),'PESS',PESS,'QESS',QESS,'E',E);
    soln=init_vgvb_solution('VB',objective_type,target_t,sense,params);
    soln=optimize_and_collect_vgvb(soln,cons,obj,sense,vars,ts,params);
end

function soln = init_vgvb_solution(resource_type, objective_type, target_t, sense, params)
    soln=struct('resource_type',resource_type,'objective_type',objective_type,'sense',sense,'target_t',target_t,'objective_value',NaN,'solver_status',NaN,'solver_info','', ...
        'PVG',NaN(1,params.T),'QVG',NaN(1,params.T),'PVB',NaN(1,params.T),'QVB',NaN(1,params.T),'EVB',NaN(1,params.T+1), ...
        'PG',NaN(1,params.T),'QG',NaN(1,params.T),'PPV',NaN(1,params.T),'QPV',NaN(1,params.T),'PESS',NaN(1,params.T),'QESS',NaN(1,params.T));
end

function soln = optimize_and_collect_vgvb(soln, cons, obj, sense, vars, ts, params)
    ops=sdpsettings('solver','gurobi','verbose',params.YALMIP_VERBOSE,'gurobi.OutputFlag',params.GUROBI_VERBOSE);
    try
        if strcmpi(sense,'max'), r=optimize(cons,-obj,ops); else, r=optimize(cons,obj,ops); end
        soln.solver_status=r.problem; soln.solver_info=r.info;
    catch ME
        soln.solver_status=100; soln.solver_info=sprintf('optimize exception: %s',ME.message); fprintf('  Optimization exception: %s\n',ME.message); drawnow; return;
    end
    if soln.solver_status~=0, fprintf('  Optimization failed: problem=%d, info=%s\n', soln.solver_status, soln.solver_info); drawnow; return; end
    val=value(obj); if isnan(val), soln.solver_status=99; soln.solver_info='YALMIP returned NaN objective value'; return; end
    soln.objective_value=val;
    if isa(vars.PG,'sdpvar'), soln.PG=value(vars.PG); soln.QG=value(vars.QG); soln.PPV=value(vars.PPV); soln.QPV=value(vars.QPV); end
    if isa(vars.PESS,'sdpvar'), soln.PESS=value(vars.PESS); soln.QESS=value(vars.QESS); end
    if isa(vars.E,'sdpvar'), soln.EVB=value(vars.E); end
    totalP=sum(ts.Pload,1); totalQ=sum(ts.Qload,1);
    if strcmpi(soln.resource_type,'VG') || isa(vars.PG,'sdpvar')
        soln.PVG=soln.PG+soln.PPV-totalP; soln.QVG=soln.QG+soln.QPV-totalQ;
    end
    if strcmpi(soln.resource_type,'VB') || isa(vars.PESS,'sdpvar')
        soln.PVB=soln.PESS; soln.QVB=soln.QESS;
    end
    fprintf('  Success: objective_value=%.6f\n', soln.objective_value); drawnow;
end

function bounds = assemble_vgvb_bounds_table(solutions, T)
    bounds=repmat(empty_vgvb_bound_row(),T,1); for t=1:T, bounds(t).t=t; end
    for i=1:numel(solutions)
        s=solutions(i); t=s.target_t; val=s.objective_value; st=s.solver_status; key='';
        switch s.objective_type
            case 'VG_power', key='PVG';
            case 'VG_ramp', key='RVG';
            case 'VB_power', key='PVB';
            case 'VB_energy', key='EVB';
        end
        if isempty(key), continue; end
        suffix='max'; if strcmpi(s.sense,'min'), suffix='min'; end
        bounds(t).([key '_' suffix])=val; bounds(t).([key '_' suffix '_status'])=st;
    end
end

function row=empty_vgvb_bound_row()
    row=struct('t',NaN,'PVG_min',NaN,'PVG_max',NaN,'RVG_min',NaN,'RVG_max',NaN,'PVB_min',NaN,'PVB_max',NaN,'EVB_min',NaN,'EVB_max',NaN, ...
        'PVG_min_status',NaN,'PVG_max_status',NaN,'RVG_min_status',NaN,'RVG_max_status',NaN,'PVB_min_status',NaN,'PVB_max_status',NaN,'EVB_min_status',NaN,'EVB_max_status',NaN);
end

function save_vgvb_bounds_csv(bounds, fname)
    writetable(struct2table(bounds), fname); fprintf('Saved CSV: %s\n', fname); drawnow;
end

function save_vgvb_solutions_csv(solutions, fname, params)
    rows=repmat(empty_vgvb_solution_csv_row(params),numel(solutions),1);
    for i=1:numel(solutions), rows(i)=vgvb_solution_to_csv_row(solutions(i),params); end
    writetable(struct2table(rows), fname); fprintf('Saved CSV: %s\n', fname); drawnow;
end

function row=empty_vgvb_solution_csv_row(params)
    row=struct('resource_type','','objective_type','','sense','','target_t',NaN,'objective_value',NaN,'solver_status',NaN,'solver_info','');
    for t=1:params.T, row.(sprintf('PVG_%d',t))=NaN; end
    for t=1:params.T, row.(sprintf('PVB_%d',t))=NaN; end
    for t=1:params.T+1, row.(sprintf('EVB_%d',t))=NaN; end
    for t=1:params.T, row.(sprintf('PG_%d',t))=NaN; end
    for t=1:params.T, row.(sprintf('PPV_%d',t))=NaN; end
    for t=1:params.T, row.(sprintf('PESS_%d',t))=NaN; end
end

function row=vgvb_solution_to_csv_row(s,params)
    row=empty_vgvb_solution_csv_row(params); row.resource_type=s.resource_type; row.objective_type=s.objective_type; row.sense=s.sense; row.target_t=s.target_t; row.objective_value=s.objective_value; row.solver_status=s.solver_status; row.solver_info=s.solver_info;
    for t=1:params.T, row.(sprintf('PVG_%d',t))=s.PVG(t); row.(sprintf('PVB_%d',t))=s.PVB(t); row.(sprintf('PG_%d',t))=s.PG(t); row.(sprintf('PPV_%d',t))=s.PPV(t); row.(sprintf('PESS_%d',t))=s.PESS(t); end
    for t=1:params.T+1, row.(sprintf('EVB_%d',t))=s.EVB(t); end
end

function comparison=assemble_vgvb_comparison(coupled,decoupled,T)
    comparison=repmat(struct(),T,1);
    fields={'PVG_min','PVG_max','RVG_min','RVG_max','PVB_min','PVB_max','EVB_min','EVB_max'};
    for t=1:T
        comparison(t).t=t;
        for i=1:numel(fields)
            f=fields{i}; comparison(t).([f '_coupled'])=coupled(t).(f); comparison(t).([f '_decoupled'])=decoupled(t).(f);
        end
        comparison(t).PVG_max_shrink=coupled(t).PVG_max-decoupled(t).PVG_max;
        comparison(t).PVG_min_shrink=decoupled(t).PVG_min-coupled(t).PVG_min;
        comparison(t).RVG_max_shrink=coupled(t).RVG_max-decoupled(t).RVG_max;
        comparison(t).RVG_min_shrink=decoupled(t).RVG_min-coupled(t).RVG_min;
        comparison(t).PVB_max_shrink=coupled(t).PVB_max-decoupled(t).PVB_max;
        comparison(t).PVB_min_shrink=decoupled(t).PVB_min-coupled(t).PVB_min;
        comparison(t).EVB_max_shrink=coupled(t).EVB_max-decoupled(t).EVB_max;
        comparison(t).EVB_min_shrink=decoupled(t).EVB_min-coupled(t).EVB_min;
    end
end

function plot_vgvb_power_bounds_compare(coupled, decoupled, params, resource_type)
    tblC=struct2table(coupled); tblD=struct2table(decoupled); t=tblC.t;
    fig=figure('Visible','off','Color','w','Position',[100 100 900 520]); ax=axes(fig); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    if strcmpi(resource_type,'VG'), minf='PVG_min'; maxf='PVG_max'; fname='vg_power_bounds_compare.png'; ttl='VG power bounds compare'; ylab='P^{VG} (MW)';
    else, minf='PVB_min'; maxf='PVB_max'; fname='vb_power_bounds_compare.png'; ttl='VB power bounds compare'; ylab='P^{VB} (MW)'; end
    [x,y]=make_step_xy(t,tblC.(minf),'post'); h1=plot(ax,x,y,'b-','LineWidth',1.6,'DisplayName','coupled min');
    [x,y]=make_step_xy(t,tblC.(maxf),'post'); h2=plot(ax,x,y,'r-','LineWidth',1.6,'DisplayName','coupled max');
    [x,y]=make_step_xy(t,tblD.(minf),'post'); h3=plot(ax,x,y,'b--','LineWidth',1.8,'DisplayName','decoupled min');
    [x,y]=make_step_xy(t,tblD.(maxf),'post'); h4=plot(ax,x,y,'r--','LineWidth',1.8,'DisplayName','decoupled max');
    yline_compat(ax,0,'k--','0'); xlabel(ax,'t'); ylabel(ax,ylab); title(ax,ttl); xlim(ax,[min(t),max(t)+1]); set(ax,'XTick',min(t):(max(t)+1)); legend(ax,[h1,h2,h3,h4],'Location','best');
    save_png_compat(fig,fullfile(params.RESULT_DIR,fname)); close(fig); fprintf('Saved figure: %s\n', fullfile(params.RESULT_DIR,fname)); drawnow;
end

function plot_vg_ramp_bounds_compare(coupled, decoupled, params)
    tblC=struct2table(coupled); tblD=struct2table(decoupled); mask=~isnan(tblC.RVG_min)|~isnan(tblD.RVG_min); t=tblC.t(mask);
    fig=figure('Visible','off','Color','w','Position',[100 100 900 520]); ax=axes(fig); hold(ax,'on'); grid(ax,'on'); box(ax,'on'); hs=gobjects(0);
    if ~isempty(t)
        [x,y]=make_step_xy(t,tblC.RVG_min(mask),'post'); hs(end+1)=plot(ax,x,y,'b-','LineWidth',1.6,'DisplayName','coupled min');
        [x,y]=make_step_xy(t,tblC.RVG_max(mask),'post'); hs(end+1)=plot(ax,x,y,'r-','LineWidth',1.6,'DisplayName','coupled max');
        [x,y]=make_step_xy(t,tblD.RVG_min(mask),'post'); hs(end+1)=plot(ax,x,y,'b--','LineWidth',1.8,'DisplayName','decoupled min');
        [x,y]=make_step_xy(t,tblD.RVG_max(mask),'post'); hs(end+1)=plot(ax,x,y,'r--','LineWidth',1.8,'DisplayName','decoupled max');
        xlim(ax,[min(t),max(t)+1]); set(ax,'XTick',min(t):(max(t)+1));
    end
    yline_compat(ax,0,'k--','0'); xlabel(ax,'t'); ylabel(ax,'P^{VG}(t)-P^{VG}(t-1) (MW/step)'); title(ax,'VG ramp bounds compare'); if ~isempty(hs), legend(ax,hs,'Location','best'); end
    save_png_compat(fig,fullfile(params.RESULT_DIR,'vg_ramp_bounds_compare.png')); close(fig); fprintf('Saved figure: %s\n', fullfile(params.RESULT_DIR,'vg_ramp_bounds_compare.png')); drawnow;
end

function plot_vb_energy_bounds_compare(coupled, decoupled, params)
    tblC=struct2table(coupled); tblD=struct2table(decoupled); t=tblC.t;
    fig=figure('Visible','off','Color','w','Position',[100 100 900 520]); ax=axes(fig); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    [x,y]=make_step_xy(t,tblC.EVB_min,'post'); h1=plot(ax,x,y,'b-','LineWidth',1.6,'DisplayName','coupled min');
    [x,y]=make_step_xy(t,tblC.EVB_max,'post'); h2=plot(ax,x,y,'r-','LineWidth',1.6,'DisplayName','coupled max');
    [x,y]=make_step_xy(t,tblD.EVB_min,'post'); h3=plot(ax,x,y,'b--','LineWidth',1.8,'DisplayName','decoupled min');
    [x,y]=make_step_xy(t,tblD.EVB_max,'post'); h4=plot(ax,x,y,'r--','LineWidth',1.8,'DisplayName','decoupled max');
    xlabel(ax,'t'); ylabel(ax,'E^{VB}_{t+1} (MWh)'); title(ax,'VB energy bounds compare'); xlim(ax,[min(t),max(t)+1]); set(ax,'XTick',min(t):(max(t)+1)); legend(ax,[h1,h2,h3,h4],'Location','best');
    save_png_compat(fig,fullfile(params.RESULT_DIR,'vb_energy_bounds_compare.png')); close(fig); fprintf('Saved figure: %s\n', fullfile(params.RESULT_DIR,'vb_energy_bounds_compare.png')); drawnow;
end

function plot_network_decoupling_margin_summary(decoup, params)
    types=unique(decoup.row_type); valsB=zeros(numel(types),1); valsG=zeros(numel(types),1);
    for i=1:numel(types)
        mask=strcmp(decoup.row_type,types{i}); valsB(i)=mean(decoup.deltaB(mask)); valsG(i)=mean(mean(decoup.deltaG_by_t(mask,:),2));
    end
    fig=figure('Visible','off','Color','w','Position',[100 100 980 520]); ax=axes(fig); bar(ax,[valsB valsG]); grid(ax,'on'); box(ax,'on'); set(ax,'XTick',1:numel(types),'XTickLabel',types,'XTickLabelRotation',30); ylabel(ax,'mean robust margin'); title(ax,'Network decoupling margin summary'); legend(ax,{'deltaB for VG','deltaG for VB'},'Location','best');
    save_png_compat(fig,fullfile(params.RESULT_DIR,'network_decoupling_margin_summary.png')); close(fig); fprintf('Saved figure: %s\n', fullfile(params.RESULT_DIR,'network_decoupling_margin_summary.png')); drawnow;
end

function save_png_compat(fig, fname)
    try
        exportgraphics(fig, fname, 'Resolution', 200);
    catch
        set(fig, 'PaperPositionMode', 'auto');
        print(fig, fname, '-dpng', '-r200');
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
