%% ieee13_temporal_p0_boundary_lindisflow_yalmip_gurobi_v1.m
% IEEE-13 simplified balanced-feeder temporal P0 trajectory-boundary Step 1.
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
%   - This version does NOT implement VG/VB decomposition, bound shrinking,
%     Farkas/KKT infeasible-point search, or neural-network training.
%   - This is only the first temporal boundary step and is not the complete
%     second stage of the paper workflow.
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
RESULT_DIR = 'results_ieee13_temporal_p0_boundary_lindisflow_v1';

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

fprintf('IEEE-13 temporal P0 boundary LinDistFlow Step 1 started at %s\n', datestr(now)); drawnow;
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
    fig = figure('Visible', 'off', 'Color', 'w');
    hold on; grid on; box on;
    fill([t; flipud(t)], [tbl.P0_min; flipud(tbl.P0_max)], [0.84 0.92 1.00], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.7, 'DisplayName', 'P_0 feasible band');
    plot(t, tbl.P0_min, 'b-o', 'LineWidth', 1.5, 'DisplayName', 'P_0 min');
    plot(t, tbl.P0_max, 'r-o', 'LineWidth', 1.5, 'DisplayName', 'P_0 max');
    yline_compat(0, 'k--', 'P_0=0');
    xlabel('t');
    ylabel('P_0 (MW)');
    title('IEEE-13 simplified feeder | temporal P0 power bounds | LinDistFlow');
    legend('Location', 'best');
    save_png_compat(fig, fullfile(params.RESULT_DIR, 'p0_power_bounds.png'));
    close(fig);
    fprintf('Saved figure: %s\n', fullfile(params.RESULT_DIR, 'p0_power_bounds.png')); drawnow;
end

function plot_p0_ramp_bounds(temporal_bounds, params)
    tbl = struct2table(temporal_bounds);
    mask = ~isnan(tbl.R0_min) | ~isnan(tbl.R0_max);
    t = tbl.t(mask);
    fig = figure('Visible', 'off', 'Color', 'w');
    hold on; grid on; box on;
    if ~isempty(t)
        fill([t; flipud(t)], [tbl.R0_min(mask); flipud(tbl.R0_max(mask))], [0.90 0.95 0.88], ...
            'EdgeColor', 'none', 'FaceAlpha', 0.7, 'DisplayName', 'R_0 feasible band');
        plot(t, tbl.R0_min(mask), 'b-o', 'LineWidth', 1.5, 'DisplayName', 'R_0 min');
        plot(t, tbl.R0_max(mask), 'r-o', 'LineWidth', 1.5, 'DisplayName', 'R_0 max');
    end
    yline_compat(0, 'k--', '0');
    xlabel('t');
    ylabel('P_0(t)-P_0(t-1) (MW/step)');
    title('IEEE-13 simplified feeder | temporal P0 ramp bounds | LinDistFlow');
    legend('Location', 'best');
    save_png_compat(fig, fullfile(params.RESULT_DIR, 'p0_ramp_bounds.png'));
    close(fig);
    fprintf('Saved figure: %s\n', fullfile(params.RESULT_DIR, 'p0_ramp_bounds.png')); drawnow;
end

function plot_p0_cum_energy_bounds(temporal_bounds, params)
    tbl = struct2table(temporal_bounds);
    t = tbl.t;
    fig = figure('Visible', 'off', 'Color', 'w');
    hold on; grid on; box on;
    fill([t; flipud(t)], [tbl.P0_cum_energy_min; flipud(tbl.P0_cum_energy_max)], [0.95 0.90 1.00], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.7, 'DisplayName', 'cumulative exchange band');
    plot(t, tbl.P0_cum_energy_min, 'b-o', 'LineWidth', 1.5, 'DisplayName', 'P0 cumulative min');
    plot(t, tbl.P0_cum_energy_max, 'r-o', 'LineWidth', 1.5, 'DisplayName', 'P0 cumulative max');
    yline_compat(0, 'k--', '0');
    xlabel('t');
    ylabel('P0 cumulative exchange energy (MWh)');
    title('IEEE-13 simplified feeder | P0 cumulative exchange energy bounds | LinDistFlow');
    legend('Location', 'best');
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
    fig = figure('Visible', 'off', 'Color', 'w');
    hold on; grid on; box on;
    plotted = 0;
    for i = 1:size(examples, 1)
        soln = find_solution(all_solutions, examples{i, 1}, examples{i, 2}, examples{i, 3});
        if isempty(soln) || soln.solver_status ~= 0 || any(isnan(soln.E))
            fprintf('Skipping ESS SOC example: %s, sense=%s, target_t=%d (missing or failed).\n', examples{i,1}, examples{i,2}, examples{i,3}); drawnow;
            continue;
        end
        stairs(0:params.T, soln.E, 'LineWidth', 1.5, 'DisplayName', examples{i, 4});
        plotted = plotted + 1;
    end
    yline_compat(params.E_min, 'k--', 'E_{min}');
    yline_compat(params.E_max, 'k--', 'E_{max}');
    xlabel('time index');
    ylabel('ESS SOC / stored energy E (MWh)');
    title('IEEE-13 simplified feeder | ESS SOC trajectories from temporal boundary examples');
    if plotted > 0
        legend('Location', 'best');
    end
    save_png_compat(fig, fullfile(params.RESULT_DIR, 'ess_soc_boundary_examples.png'));
    close(fig);
    fprintf('Saved figure: %s\n', fullfile(params.RESULT_DIR, 'ess_soc_boundary_examples.png')); drawnow;
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

function yline_compat(y, style, label_text)
    try
        yline(y, style, label_text);
    catch
        xl = xlim;
        plot(xl, [y y], style, 'DisplayName', label_text);
    end
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
