%% ieee9_pq_flex_scan_yalmip_gurobi_v1.m
% Small-scale IEEE-9 (P-Q) flexibility-domain scan using MATLAB + YALMIP + Gurobi.
%
% This script is a first-stage reproduction/verification example for scanning
% the per-period PCC active/reactive exchange section. It deliberately avoids
% nonlinear AC-OPF and instead builds finite-difference AC power-flow
% sensitivities around period-dependent base cases, then solves linear/SOCP
% YALMIP models for many scanning directions.
%
% Sign convention used in every CSV/figure produced by this script:
%   P0 = -(net active injection at PCC/slack bus 1 from MATPOWER PF result)
%   Q0 = -(net reactive injection at PCC/slack bus 1 from MATPOWER PF result)
% Therefore P0 > 0 and Q0 > 0 mean the IEEE-9 internal system exports active
% and reactive power to the upstream/TSO grid. MATPOWER slack injection itself
% has the opposite sign: a positive slack injection means import from upstream.
%
% Important modeling note:
%   The PCC/slack generator at bus 1 is treated as the external grid equivalent,
%   not as an internal controllable unit. The non-PCC IEEE-9 generators are used
%   as controllable internal units. Generator/ESS injection changes are applied
%   as equivalent net-injection changes in the finite-difference and AC
%   validation power-flow cases.
%
% Interpretation note:
%   The obtained (P0,Q0) regions are per-period sections. In multi_period mode,
%   the objective scans only one period's alpha*P0+beta*Q0, but constraints span
%   all periods through generator ramping and ESS SOC recursion. These sections
%   are still not the full high-dimensional temporal flexibility region; future
%   extensions should combine them with active-power trajectory boundaries,
%   VG/VB decomposition, and bound shrinking.

clear; clc; close all;

%% --------------------------- User parameters ---------------------------
RESULT_DIR = 'results_ieee9_pq_flex_scan_v1';
SCAN_MODE = 'fixed_angle';              % 'fixed_angle' or 'adaptive_icpf'
RUN_TIME_COUPLING_MODES = {'single_period','multi_period'};
DEFAULT_TIME_COUPLING_MODE = 'multi_period'; %#ok<NASGU> % documented default
N_THETA = 36;
ADAPTIVE_INITIAL_N_THETA = 12;
ADAPTIVE_MAX_ROUNDS = 3;
ADAPTIVE_MAX_POINTS = 72;
ADAPTIVE_DIST_TOL_MVA = 15.0;
ADAPTIVE_Q_TOL_MVAR = 10.0;
YALMIP_VERBOSE = 0;
GUROBI_VERBOSE = 0;

T = 4;
LOAD_SCALE = [0.90, 1.00, 1.10, 0.95];
DT = 1.0;                               % h
EPS_SENS = 1e-3;                        % MW/MVAr finite-difference step
VMIN = 0.95;
VMAX = 1.05;
BRANCH_RATE_DEFAULT = 250;              % MVA/MW proxy when MATPOWER rateA is 0
BRANCH_RATE_SCALE = 1.00;
RAMP_FRACTION_UP = 0.30;                % default Pg capacity fraction per period
RAMP_FRACTION_DOWN = 0.30;

% Energy storage system, connected as an equivalent controllable injection.
ESS_BUS = 5;
ESS_E_MAX = 80;                         % MWh
ESS_E_MIN = 10;
ESS_E0 = 40;
ESS_P_CH_MAX = 30;                      % MW
ESS_P_DIS_MAX = 30;
ESS_Q_MAX = 25;                         % MVAr
ESS_S_MAX = 35;                         % MVA
ETA_CH = 0.95;
ETA_DIS = 0.95;
USE_ESS_SOC_CONE = false;               % false -> robust linear box approximation
USE_ESS_BINARY = false;                 % reserved; false keeps continuous relaxation

% Future-extension placeholders for paper-style second-stage refinements.
ENABLE_RAMP_BOUNDARY_EXTENSION = false; %#ok<NASGU>
ENABLE_ENERGY_BOUNDARY_EXTENSION = false; %#ok<NASGU>
ENABLE_VG_VB_DECOMPOSITION = false; %#ok<NASGU>
ENABLE_BOUND_SHRINKING = false; %#ok<NASGU>

%% ------------------------- Initialization ------------------------------
if exist(RESULT_DIR, 'dir')
    warning('Result directory already exists. New run outputs may overwrite files inside %s.', RESULT_DIR);
else
    mkdir(RESULT_DIR);
end
log_file = fullfile(RESULT_DIR, 'run_log.txt');
if exist(log_file, 'file')
    delete(log_file);
end
diary(log_file);
diary on;
cleanup_diary = onCleanup(@() diary('off'));

fprintf('IEEE-9 P-Q flexibility scan started at %s\n', datestr(now)); drawnow;
fprintf('Result directory: %s\n', RESULT_DIR); drawnow;
fprintf('Scan mode: %s, N_THETA=%d, T=%d\n', SCAN_MODE, N_THETA, T); drawnow;

check_dependencies();
define_constants;
mpc0 = load_ieee9_case();
[time_cases, base_pf, base, sys] = build_time_series_cases(mpc0, T, LOAD_SCALE, ESS_BUS, BRANCH_RATE_DEFAULT, BRANCH_RATE_SCALE);
params = struct();
params.RESULT_DIR = RESULT_DIR;
params.SCAN_MODE = SCAN_MODE;
params.N_THETA = N_THETA;
params.ADAPTIVE_INITIAL_N_THETA = ADAPTIVE_INITIAL_N_THETA;
params.ADAPTIVE_MAX_ROUNDS = ADAPTIVE_MAX_ROUNDS;
params.ADAPTIVE_MAX_POINTS = ADAPTIVE_MAX_POINTS;
params.ADAPTIVE_DIST_TOL_MVA = ADAPTIVE_DIST_TOL_MVA;
params.ADAPTIVE_Q_TOL_MVAR = ADAPTIVE_Q_TOL_MVAR;
params.YALMIP_VERBOSE = YALMIP_VERBOSE;
params.GUROBI_VERBOSE = GUROBI_VERBOSE;
params.T = T;
params.DT = DT;
params.EPS_SENS = EPS_SENS;
params.VMIN = VMIN;
params.VMAX = VMAX;
params.RAMP_FRACTION_UP = RAMP_FRACTION_UP;
params.RAMP_FRACTION_DOWN = RAMP_FRACTION_DOWN;
params.ESS_BUS = ESS_BUS;
params.ESS_E_MAX = ESS_E_MAX;
params.ESS_E_MIN = ESS_E_MIN;
params.ESS_E0 = ESS_E0;
params.ESS_P_CH_MAX = ESS_P_CH_MAX;
params.ESS_P_DIS_MAX = ESS_P_DIS_MAX;
params.ESS_Q_MAX = ESS_Q_MAX;
params.ESS_S_MAX = ESS_S_MAX;
params.ETA_CH = ETA_CH;
params.ETA_DIS = ETA_DIS;
params.USE_ESS_SOC_CONE = USE_ESS_SOC_CONE;
params.USE_ESS_BINARY = USE_ESS_BINARY;

fprintf('Computing finite-difference sensitivities for all periods ...\n'); drawnow;
sens = compute_sensitivities_for_all_t(time_cases, base_pf, base, sys, params);

all_rows = struct([]);
all_solutions = struct();
summary_rows = struct([]);

%% ------------------------- Main scanning loop --------------------------
for imode = 1:numel(RUN_TIME_COUPLING_MODES)
    mode = RUN_TIME_COUPLING_MODES{imode};
    fprintf('\n=== Running TIME_COUPLING_MODE=%s ===\n', mode); drawnow;
    all_solutions.(mode) = cell(T, 1);

    for target_t = 1:T
        fprintf('\n--- Scanning mode=%s, period t=%d/%d ---\n', mode, target_t, T); drawnow;
        theta_list = get_scan_angles(params, []);
        period_solutions = struct([]);
        period_rows = struct([]);

        if strcmpi(SCAN_MODE, 'adaptive_icpf')
            [period_rows, period_solutions] = run_adaptive_scan(mode, target_t, theta_list, time_cases, base, sys, sens, params);
        else
            for k = 1:numel(theta_list)
                theta = theta_list(k);
                fprintf('Solving mode=%s, t=%d/%d, theta=%.1f deg, idx=%d/%d ...\n', ...
                    mode, target_t, T, rad2deg(theta), k, numel(theta_list)); drawnow;
                soln = solve_pq_boundary_point_yalmip(mode, target_t, theta, time_cases, base, sys, sens, params);
                [row, soln] = build_result_row_and_validate(mode, target_t, theta, soln, time_cases, base, sys, params);
                period_rows = append_struct(period_rows, row);
                period_solutions = append_struct(period_solutions, soln);
            end
        end

        all_rows = append_struct(all_rows, period_rows);
        all_solutions.(mode){target_t} = period_solutions;
        summary_row = summarize_period(mode, target_t, period_rows);
        summary_rows = append_struct(summary_rows, summary_row);

        plot_pq_flex_region(period_rows, mode, target_t, params);
    end
end

%% ------------------------- Save and compare ----------------------------
save_results_csv(all_rows, summary_rows, RESULT_DIR);
save(fullfile(RESULT_DIR, 'dispatch_solutions.mat'), 'all_solutions', 'summary_rows', 'all_rows', 'params', 'sys', 'base', 'sens', '-v7.3');
plot_compare_all_periods(all_rows, RUN_TIME_COUPLING_MODES, T, params);

fprintf('\nCompleted IEEE-9 P-Q flexibility scan at %s\n', datestr(now)); drawnow;
fprintf('CSV, MAT, PNG, and log outputs saved under %s\n', RESULT_DIR); drawnow;

%% ============================= Local functions =========================
function check_dependencies()
    missing = {};
    if exist('case9', 'file') ~= 2
        missing{end+1} = 'MATPOWER case9 (install MATPOWER and add it to the MATLAB path, or keep a compatible case9.m on the path)'; %#ok<AGROW>
    end
    if exist('runpf', 'file') ~= 2 || exist('mpoption', 'file') ~= 2
        missing{end+1} = 'MATPOWER runpf/mpoption (install MATPOWER and addpath(genpath(matpower_root)))'; %#ok<AGROW>
    end
    if exist('sdpvar', 'file') ~= 2 || exist('optimize', 'file') ~= 2 || exist('sdpsettings', 'file') ~= 2
        missing{end+1} = 'YALMIP (install YALMIP and add it to the MATLAB path)'; %#ok<AGROW>
    end
    if exist('gurobi', 'file') ~= 3 && exist('gurobi', 'file') ~= 2 && exist('gurobi_mex', 'file') ~= 3 && exist('gurobi_mex', 'file') ~= 2
        missing{end+1} = 'Gurobi MATLAB interface (install Gurobi, activate a license, and add its MATLAB folder to the path)'; %#ok<AGROW>
    end
    if ~isempty(missing)
        fprintf('Dependency check failed. Missing components:\n');
        for i = 1:numel(missing)
            fprintf('  - %s\n', missing{i});
        end
        error('Please install/add the missing dependencies before running this script.');
    end
    fprintf('Dependency check passed: MATPOWER, YALMIP, and Gurobi are visible on the MATLAB path.\n'); drawnow;
end

function mpc = load_ieee9_case()
    mpc = case9();
    fprintf('Loaded IEEE-9 case: %d buses, %d branches, %d generators.\n', size(mpc.bus,1), size(mpc.branch,1), size(mpc.gen,1)); drawnow;
end

function [time_cases, base_pf, base, sys] = build_time_series_cases(mpc0, T, load_scale, ess_bus, branch_rate_default, branch_rate_scale)
    define_constants;
    if numel(load_scale) ~= T
        error('LOAD_SCALE length (%d) must equal T (%d).', numel(load_scale), T);
    end
    sys = struct();
    sys.pcc_bus = 1;
    sys.ess_bus = ess_bus;
    sys.nb = size(mpc0.bus, 1);
    sys.ng_all = size(mpc0.gen, 1);
    sys.nl = size(mpc0.branch, 1);
    sys.gen_bus_all = mpc0.gen(:, GEN_BUS);
    sys.ctrl_gen_idx = find(mpc0.gen(:, GEN_BUS) ~= sys.pcc_bus);
    sys.ng = numel(sys.ctrl_gen_idx);
    sys.ctrl_gen_bus = mpc0.gen(sys.ctrl_gen_idx, GEN_BUS);
    if sys.ng == 0
        error('No non-PCC internal generators found. This script treats bus-1 slack as external exchange.');
    end
    sys.branch_rate = mpc0.branch(:, RATE_A);
    sys.branch_rate(sys.branch_rate <= 0) = branch_rate_default;
    sys.branch_rate = branch_rate_scale * sys.branch_rate;

    mpopt = mpoption('verbose', 0, 'out.all', 0);
    time_cases = cell(T, 1);
    base_pf = cell(T, 1);
    base = repmat(struct(), T, 1);
    for t = 1:T
        mpc = mpc0;
        mpc.bus(:, PD) = mpc0.bus(:, PD) * load_scale(t);
        mpc.bus(:, QD) = mpc0.bus(:, QD) * load_scale(t);
        pf = runpf(mpc, mpopt);
        if ~pf.success
            error('Base AC PF failed for period t=%d. Check case data and MATPOWER installation.', t);
        end
        time_cases{t} = mpc;
        base_pf{t} = pf;
        base(t).Pg = pf.gen(sys.ctrl_gen_idx, PG);
        base(t).Qg = pf.gen(sys.ctrl_gen_idx, QG);
        base(t).P0 = calc_pcc_exchange(pf, sys.pcc_bus, 'P');
        base(t).Q0 = calc_pcc_exchange(pf, sys.pcc_bus, 'Q');
        base(t).V = pf.bus(:, VM);
        base(t).Pf = pf.branch(:, PF);
        base(t).Pt = pf.branch(:, PT);
        base(t).Qf = pf.branch(:, QF);
        base(t).Qt = pf.branch(:, QT);
        fprintf('Base PF t=%d: success=1, P0=%.4f MW, Q0=%.4f MVAr.\n', t, base(t).P0, base(t).Q0); drawnow;
    end
end

function sens = compute_sensitivities_for_all_t(time_cases, base_pf, base, sys, params)
    T = numel(time_cases);
    sens = repmat(struct(), T, 1);
    for t = 1:T
        fprintf('  Sensitivity t=%d/%d ...\n', t, T); drawnow;
        sens(t) = compute_sensitivity_one_t(time_cases{t}, base_pf{t}, base(t), sys, params);
    end
end

function st = compute_sensitivity_one_t(mpc_base, pf_base, base_t, sys, params)
    define_constants;
    nctrl = 2 * sys.ng + 2;
    ny = 2 + sys.nb + 2 * sys.nl;
    y0 = pack_outputs(pf_base, sys);
    S = zeros(ny, nctrl);
    mpopt = mpoption('verbose', 0, 'out.all', 0);

    for j = 1:nctrl
        mpc = mpc_base;
        if j <= sys.ng
            bus = sys.ctrl_gen_bus(j);
            mpc.bus(bus, PD) = mpc.bus(bus, PD) - params.EPS_SENS;
        elseif j <= 2 * sys.ng
            bus = sys.ctrl_gen_bus(j - sys.ng);
            mpc.bus(bus, QD) = mpc.bus(bus, QD) - params.EPS_SENS;
        elseif j == 2 * sys.ng + 1
            mpc.bus(sys.ess_bus, PD) = mpc.bus(sys.ess_bus, PD) - params.EPS_SENS;
        else
            mpc.bus(sys.ess_bus, QD) = mpc.bus(sys.ess_bus, QD) - params.EPS_SENS;
        end
        pf = runpf(mpc, mpopt);
        if ~pf.success
            warning('Sensitivity PF failed for control j=%d; using zero column.', j);
            S(:, j) = 0;
        else
            S(:, j) = (pack_outputs(pf, sys) - y0) / params.EPS_SENS;
        end
    end
    st = struct();
    st.y0 = y0;
    st.S = S;
    st.SP0 = S(1, :);
    st.SQ0 = S(2, :);
    st.SV = S(3:(2+sys.nb), :);
    st.SPf = S((3+sys.nb):(2+sys.nb+sys.nl), :);
    st.SPt = S((3+sys.nb+sys.nl):(2+sys.nb+2*sys.nl), :);
    st.base = base_t;
end

function y = pack_outputs(pf, sys)
    define_constants;
    y = [calc_pcc_exchange(pf, sys.pcc_bus, 'P'); ...
         calc_pcc_exchange(pf, sys.pcc_bus, 'Q'); ...
         pf.bus(:, VM); ...
         pf.branch(:, PF); ...
         pf.branch(:, PT)];
end

function val = calc_pcc_exchange(pf, pcc_bus, pqflag)
    define_constants;
    gen_at_pcc = find(pf.gen(:, GEN_BUS) == pcc_bus);
    bus_row = find(pf.bus(:, BUS_I) == pcc_bus, 1);
    if isempty(bus_row)
        error('PCC bus %d not found in PF result.', pcc_bus);
    end
    if strcmpi(pqflag, 'P')
        net_slack_injection = sum(pf.gen(gen_at_pcc, PG)) - pf.bus(bus_row, PD);
    else
        net_slack_injection = sum(pf.gen(gen_at_pcc, QG)) - pf.bus(bus_row, QD);
    end
    val = -net_slack_injection;
end

function theta_list = get_scan_angles(params, seed_angles)
    if strcmpi(params.SCAN_MODE, 'adaptive_icpf')
        if nargin >= 2 && ~isempty(seed_angles)
            theta_list = unique(mod(seed_angles(:), 2*pi)).';
        else
            theta_list = linspace(0, 2*pi, params.ADAPTIVE_INITIAL_N_THETA + 1);
            theta_list(end) = [];
        end
    else
        theta_list = linspace(0, 2*pi, params.N_THETA + 1);
        theta_list(end) = [];
    end
end

function [period_rows, period_solutions] = run_adaptive_scan(mode, target_t, theta_list, time_cases, base, sys, sens, params)
    solved = containers.Map('KeyType', 'char', 'ValueType', 'any');
    period_rows = struct([]);
    period_solutions = struct([]);
    current = unique(mod(theta_list(:), 2*pi)).';
    for round_id = 1:params.ADAPTIVE_MAX_ROUNDS
        current = sort(current);
        for k = 1:numel(current)
            key = sprintf('%.12f', current(k));
            if ~isKey(solved, key)
                fprintf('Solving adaptive mode=%s, t=%d/%d, round=%d, theta=%.1f deg, idx=%d/%d ...\n', ...
                    mode, target_t, params.T, round_id, rad2deg(current(k)), k, numel(current)); drawnow;
                soln = solve_pq_boundary_point_yalmip(mode, target_t, current(k), time_cases, base, sys, sens, params);
                [row, soln] = build_result_row_and_validate(mode, target_t, current(k), soln, time_cases, base, sys, params);
                solved(key) = {row, soln};
            end
        end
        rows_tmp = collect_rows_from_map(solved);
        new_angles = propose_adaptive_angles(rows_tmp, params);
        if isempty(new_angles) || numel(current) >= params.ADAPTIVE_MAX_POINTS
            break;
        end
        current = unique(mod([current(:); new_angles(:)], 2*pi)).';
    end
    keys = solved.keys;
    theta_vals = zeros(numel(keys), 1);
    for i = 1:numel(keys), theta_vals(i) = str2double(keys{i}); end
    [~, order] = sort(theta_vals);
    for ii = 1:numel(order)
        pair = solved(keys{order(ii)});
        period_rows = append_struct(period_rows, pair{1});
        period_solutions = append_struct(period_solutions, pair{2});
    end
end

function rows = collect_rows_from_map(mp)
    rows = struct([]);
    keys = mp.keys;
    for i = 1:numel(keys)
        pair = mp(keys{i});
        rows = append_struct(rows, pair{1});
    end
end

function new_angles = propose_adaptive_angles(rows, params)
    new_angles = [];
    if isempty(rows), return; end
    tbl = struct2table(rows);
    tbl = tbl(~isnan(tbl.P0_pred) & ~isnan(tbl.Q0_pred), :);
    if height(tbl) < 2, return; end
    tbl = sortrows(tbl, 'theta_deg');
    theta = deg2rad(tbl.theta_deg);
    P = tbl.P0_pred;
    Q = tbl.Q0_pred;
    n = numel(theta);
    for i = 1:n
        j = i + 1;
        if j > n, j = 1; end
        d = hypot(P(i) - P(j), Q(i) - Q(j));
        dq = abs(Q(i) - Q(j));
        if d > params.ADAPTIVE_DIST_TOL_MVA || dq > params.ADAPTIVE_Q_TOL_MVAR
            th_mid = angle_mid(theta(i), theta(j));
            new_angles(end+1) = th_mid; %#ok<AGROW>
        end
    end
    new_angles = unique(new_angles);
end

function th = angle_mid(a, b)
    if b < a
        b = b + 2*pi;
    end
    th = mod((a + b) / 2, 2*pi);
end

function soln = solve_pq_boundary_point_yalmip(mode, target_t, theta, time_cases, base, sys, sens, params)
    define_constants;
    alpha = cos(theta);
    beta = sin(theta);
    Tmodel = params.T;
    if strcmpi(mode, 'single_period')
        periods = target_t;
    else
        periods = 1:Tmodel;
    end
    np = numel(periods);
    ng = sys.ng;
    nl = sys.nl;
    nb = sys.nb;

    Pg = sdpvar(ng, np, 'full');
    Qg = sdpvar(ng, np, 'full');
    Pch = sdpvar(1, np, 'full');
    Pdis = sdpvar(1, np, 'full');
    Qess = sdpvar(1, np, 'full');
    if params.USE_ESS_BINARY
        zch = binvar(1, np, 'full');
        zdis = binvar(1, np, 'full');
    else
        zch = [];
        zdis = [];
    end
    if strcmpi(mode, 'multi_period')
        E = sdpvar(1, Tmodel + 1, 'full');
    else
        E = sdpvar(1, np, 'full');
    end

    P0_expr = cell(Tmodel, 1);
    Q0_expr = cell(Tmodel, 1);
    V_expr = cell(Tmodel, 1);
    Pf_expr = cell(Tmodel, 1);
    Pt_expr = cell(Tmodel, 1);
    cons = [];

    Pg_min = time_cases{1}.gen(sys.ctrl_gen_idx, PMIN);
    Pg_max = time_cases{1}.gen(sys.ctrl_gen_idx, PMAX);
    Qg_min = time_cases{1}.gen(sys.ctrl_gen_idx, QMIN);
    Qg_max = time_cases{1}.gen(sys.ctrl_gen_idx, QMAX);
    ramp_up = params.RAMP_FRACTION_UP * max(Pg_max - Pg_min, 1);
    ramp_down = params.RAMP_FRACTION_DOWN * max(Pg_max - Pg_min, 1);

    for kk = 1:np
        t = periods(kk);
        Pess = Pdis(kk) - Pch(kk);
        u = [Pg(:, kk) - base(t).Pg; Qg(:, kk) - base(t).Qg; Pess; Qess(kk)];
        P0_expr{t} = base(t).P0 + sens(t).SP0 * u;
        Q0_expr{t} = base(t).Q0 + sens(t).SQ0 * u;
        V_expr{t} = base(t).V + sens(t).SV * u;
        Pf_expr{t} = base(t).Pf + sens(t).SPf * u;
        Pt_expr{t} = base(t).Pt + sens(t).SPt * u;

        % Generator constraints.
        cons = [cons, Pg_min <= Pg(:, kk), Pg(:, kk) <= Pg_max]; %#ok<AGROW>
        cons = [cons, Qg_min <= Qg(:, kk), Qg(:, kk) <= Qg_max]; %#ok<AGROW>

        % Storage power/SOC constraints.
        cons = [cons, 0 <= Pch(kk), Pch(kk) <= params.ESS_P_CH_MAX]; %#ok<AGROW>
        cons = [cons, 0 <= Pdis(kk), Pdis(kk) <= params.ESS_P_DIS_MAX]; %#ok<AGROW>
        cons = [cons, -params.ESS_Q_MAX <= Qess(kk), Qess(kk) <= params.ESS_Q_MAX]; %#ok<AGROW>
        cons = [cons, -params.ESS_S_MAX <= Pess, Pess <= params.ESS_S_MAX]; %#ok<AGROW>
        if params.USE_ESS_SOC_CONE
            cons = [cons, cone([Pess; Qess(kk)], params.ESS_S_MAX)]; %#ok<AGROW>
        end
        if params.USE_ESS_BINARY
            cons = [cons, Pch(kk) <= params.ESS_P_CH_MAX * zch(kk), Pdis(kk) <= params.ESS_P_DIS_MAX * zdis(kk), zch(kk) + zdis(kk) <= 1]; %#ok<AGROW>
        end

        % Linearized network constraints. Branch limits use signed active-flow
        % limits on both branch ends as a first-version linear approximation.
        cons = [cons, params.VMIN <= V_expr{t}, V_expr{t} <= params.VMAX]; %#ok<AGROW>
        cons = [cons, -sys.branch_rate <= Pf_expr{t}, Pf_expr{t} <= sys.branch_rate]; %#ok<AGROW>
        cons = [cons, -sys.branch_rate <= Pt_expr{t}, Pt_expr{t} <= sys.branch_rate]; %#ok<AGROW>
    end

    % Time-coupling constraints.
    if strcmpi(mode, 'multi_period')
        cons = [cons, E(1) == params.ESS_E0];
        for t = 1:Tmodel
            kk = find(periods == t, 1);
            cons = [cons, params.ESS_E_MIN <= E(t), E(t) <= params.ESS_E_MAX]; %#ok<AGROW>
            cons = [cons, E(t+1) == E(t) + params.ETA_CH * Pch(kk) * params.DT - (1 / params.ETA_DIS) * Pdis(kk) * params.DT]; %#ok<AGROW>
            if t == 1
                prev_pg = base(1).Pg;
            else
                prev_pg = Pg(:, kk - 1);
            end
            cons = [cons, -ramp_down <= Pg(:, kk) - prev_pg, Pg(:, kk) - prev_pg <= ramp_up]; %#ok<AGROW>
        end
        cons = [cons, params.ESS_E_MIN <= E(Tmodel + 1), E(Tmodel + 1) <= params.ESS_E_MAX];
    else
        cons = [cons, params.ESS_E_MIN <= E(1), E(1) <= params.ESS_E_MAX];
        % The single-period case is independent and keeps no cross-period SOC recursion.
    end

    if strcmpi(mode, 'single_period')
        obj_expr = alpha * P0_expr{target_t} + beta * Q0_expr{target_t};
    else
        obj_expr = alpha * P0_expr{target_t} + beta * Q0_expr{target_t};
    end
    ops = sdpsettings('solver', 'gurobi', 'verbose', params.YALMIP_VERBOSE, 'gurobi.OutputFlag', params.GUROBI_VERBOSE);

    sol = optimize(cons, -obj_expr, ops);
    soln = struct();
    soln.mode = mode;
    soln.t = target_t;
    soln.theta_deg = rad2deg(theta);
    soln.alpha = alpha;
    soln.beta = beta;
    soln.problem = sol.problem;
    soln.solver_info = sol.info;
    soln.objective = NaN;
    soln.P0_pred = NaN;
    soln.Q0_pred = NaN;
    soln.Pg = NaN(ng, Tmodel);
    soln.Qg = NaN(ng, Tmodel);
    soln.Pch = NaN(1, Tmodel);
    soln.Pdis = NaN(1, Tmodel);
    soln.Pess = NaN(1, Tmodel);
    soln.Qess = NaN(1, Tmodel);
    soln.E = NaN(1, Tmodel + 1);
    soln.V_pred = NaN(nb, Tmodel);
    soln.Pf_pred = NaN(nl, Tmodel);
    soln.Pt_pred = NaN(nl, Tmodel);

    if sol.problem ~= 0
        fprintf('  Optimization failed: problem=%d, info=%s\n', sol.problem, sol.info); drawnow;
        return;
    end

    val_obj = value(obj_expr);
    if isnan(val_obj)
        soln.problem = 99;
        soln.solver_info = 'YALMIP returned NaN objective/value';
        fprintf('  Optimization returned NaN values; marking as failed.\n'); drawnow;
        return;
    end
    soln.objective = val_obj;
    soln.P0_pred = value(P0_expr{target_t});
    soln.Q0_pred = value(Q0_expr{target_t});
    val_Pg = value(Pg);
    val_Qg = value(Qg);
    val_Pch = value(Pch);
    val_Pdis = value(Pdis);
    val_Qess = value(Qess);
    if any(isnan([val_Pg(:); val_Qg(:); val_Pch(:); val_Pdis(:); val_Qess(:)]))
        soln.problem = 98;
        soln.solver_info = 'YALMIP returned NaN dispatch values';
        fprintf('  Optimization dispatch values contain NaN; marking as failed.\n'); drawnow;
        return;
    end
    for kk = 1:np
        t = periods(kk);
        soln.Pg(:, t) = val_Pg(:, kk);
        soln.Qg(:, t) = val_Qg(:, kk);
        soln.Pch(t) = val_Pch(kk);
        soln.Pdis(t) = val_Pdis(kk);
        soln.Pess(t) = val_Pdis(kk) - val_Pch(kk);
        soln.Qess(t) = val_Qess(kk);
        soln.V_pred(:, t) = value(V_expr{t});
        soln.Pf_pred(:, t) = value(Pf_expr{t});
        soln.Pt_pred(:, t) = value(Pt_expr{t});
    end
    if strcmpi(mode, 'multi_period')
        soln.E = value(E);
    else
        soln.E(target_t) = value(E(1));
    end
    fprintf('  Success: P0_pred=%.4f, Q0_pred=%.4f, objective=%.4f\n', soln.P0_pred, soln.Q0_pred, soln.objective); drawnow;
end

function [row, soln] = build_result_row_and_validate(mode, target_t, theta, soln, time_cases, base, sys, params)
    row = empty_result_row();
    row.mode = mode;
    row.t = target_t;
    row.theta_deg = rad2deg(theta);
    row.alpha = cos(theta);
    row.beta = sin(theta);
    row.P0_pred = soln.P0_pred;
    row.Q0_pred = soln.Q0_pred;
    row.objective = soln.objective;
    row.solver_status = soln.problem;
    row.solver_info = soln.solver_info;
    if soln.problem == 0
        ac = validate_solution_ac_pf(soln, target_t, time_cases{target_t}, base(target_t), sys, params);
        soln.ac_validation = ac;
        row.P0_ac_validation = ac.P0_ac;
        row.Q0_ac_validation = ac.Q0_ac;
        row.ac_pf_success = ac.success;
        row.ac_validation_violation = ac.violation;
        row.max_voltage_violation = ac.max_voltage_violation;
        row.max_branch_violation = ac.max_branch_violation;
    else
        row.P0_ac_validation = NaN;
        row.Q0_ac_validation = NaN;
        row.ac_pf_success = false;
        row.ac_validation_violation = true;
        row.max_voltage_violation = NaN;
        row.max_branch_violation = NaN;
    end
end

function ac = validate_solution_ac_pf(soln, target_t, mpc_base, base_t, sys, params)
    define_constants;
    mpc = mpc_base;
    for i = 1:sys.ng
        bus = sys.ctrl_gen_bus(i);
        dP = soln.Pg(i, target_t) - base_t.Pg(i);
        dQ = soln.Qg(i, target_t) - base_t.Qg(i);
        mpc.bus(bus, PD) = mpc.bus(bus, PD) - dP;
        mpc.bus(bus, QD) = mpc.bus(bus, QD) - dQ;
    end
    mpc.bus(sys.ess_bus, PD) = mpc.bus(sys.ess_bus, PD) - soln.Pess(target_t);
    mpc.bus(sys.ess_bus, QD) = mpc.bus(sys.ess_bus, QD) - soln.Qess(target_t);
    mpopt = mpoption('verbose', 0, 'out.all', 0);
    pf = runpf(mpc, mpopt);
    ac = struct();
    ac.success = logical(pf.success);
    if pf.success
        ac.P0_ac = calc_pcc_exchange(pf, sys.pcc_bus, 'P');
        ac.Q0_ac = calc_pcc_exchange(pf, sys.pcc_bus, 'Q');
        ac.P0_error = ac.P0_ac - soln.P0_pred;
        ac.Q0_error = ac.Q0_ac - soln.Q0_pred;
        v_low = max(params.VMIN - pf.bus(:, VM), 0);
        v_high = max(pf.bus(:, VM) - params.VMAX, 0);
        ac.max_voltage_violation = max([v_low; v_high]);
        branch_abs = max(abs([pf.branch(:, PF), pf.branch(:, PT)]), [], 2);
        ac.max_branch_violation = max(max(branch_abs - sys.branch_rate, 0));
        ac.violation = ac.max_voltage_violation > 1e-5 || ac.max_branch_violation > 1e-5;
        fprintf('  AC validation: success=1, P0=%.4f, Q0=%.4f, dP0=%.4e, dQ0=%.4e, viol=%d\n', ...
            ac.P0_ac, ac.Q0_ac, ac.P0_error, ac.Q0_error, ac.violation); drawnow;
    else
        ac.P0_ac = NaN;
        ac.Q0_ac = NaN;
        ac.P0_error = NaN;
        ac.Q0_error = NaN;
        ac.max_voltage_violation = Inf;
        ac.max_branch_violation = Inf;
        ac.violation = true;
        fprintf('  AC validation PF failed.\n'); drawnow;
    end
end

function row = empty_result_row()
    row = struct('mode', '', 't', NaN, 'theta_deg', NaN, 'alpha', NaN, 'beta', NaN, ...
        'P0_pred', NaN, 'Q0_pred', NaN, 'P0_ac_validation', NaN, 'Q0_ac_validation', NaN, ...
        'objective', NaN, 'solver_status', NaN, 'solver_info', '', 'ac_pf_success', false, ...
        'ac_validation_violation', true, 'max_voltage_violation', NaN, 'max_branch_violation', NaN);
end

function summary = summarize_period(mode, target_t, rows)
    summary = struct('mode', mode, 't', target_t, 'num_successful_points', 0, 'polygon_area', NaN, ...
        'num_failed_optimizations', 0, 'num_ac_validation_violations', 0);
    if isempty(rows)
        return;
    end
    tbl = struct2table(rows);
    ok = tbl.solver_status == 0 & ~isnan(tbl.P0_pred) & ~isnan(tbl.Q0_pred);
    summary.num_successful_points = sum(ok);
    summary.num_failed_optimizations = sum(tbl.solver_status ~= 0);
    summary.num_ac_validation_violations = sum(tbl.ac_validation_violation);
    if summary.num_successful_points >= 3
        [P, Q] = sort_boundary_points(tbl.P0_pred(ok), tbl.Q0_pred(ok));
        summary.polygon_area = polygon_area_pq(P, Q);
    end
    fprintf('Summary mode=%s, t=%d: success=%d, failed=%d, AC violations=%d, area=%.4f\n', ...
        mode, target_t, summary.num_successful_points, summary.num_failed_optimizations, ...
        summary.num_ac_validation_violations, summary.polygon_area); drawnow;
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

function plot_pq_flex_region(rows, mode, target_t, params)
    if isempty(rows), return; end
    tbl = struct2table(rows);
    ok = tbl.solver_status == 0 & ~isnan(tbl.P0_pred) & ~isnan(tbl.Q0_pred);
    if sum(ok) < 1, return; end
    P = tbl.P0_pred(ok);
    Q = tbl.Q0_pred(ok);
    area = NaN;
    if sum(ok) >= 3
        [Ps, Qs] = sort_boundary_points(P, Q);
        area = polygon_area_pq(Ps, Qs);
    else
        Ps = P; Qs = Q;
    end
    fig = figure('Visible', 'off', 'Color', 'w');
    hold on; grid on; box on;
    if numel(Ps) >= 3
        fill([Ps; Ps(1)], [Qs; Qs(1)], [0.80 0.90 1.00], 'EdgeColor', [0 0.25 0.85], 'LineWidth', 1.5, 'FaceAlpha', 0.55, 'DisplayName', 'polygon');
    end
    plot(P, Q, 'o', 'Color', [0.1 0.1 0.1], 'MarkerFaceColor', [0.1 0.4 0.9], 'DisplayName', 'boundary points');
    xlabel('P_0 export to upstream grid (MW)');
    ylabel('Q_0 export to upstream grid (MVAr)');
    title(sprintf('IEEE-9 P-Q flex region, t=%02d, %s, %s, N=%d, area=%.3f', ...
        target_t, strrep(mode, '_', '\_'), strrep(params.SCAN_MODE, '_', '\_'), sum(ok), area));
    legend('Location', 'best');
    fname = fullfile(params.RESULT_DIR, sprintf('pq_flex_t%02d_%s.png', target_t, mode));
    save_png_compat(fig, fname);
    close(fig);
    fprintf('Saved figure: %s\n', fname); drawnow;
end

function plot_compare_all_periods(all_rows, modes, T, params)
    if isempty(all_rows), return; end
    tbl_all = struct2table(all_rows);
    colors = lines(numel(modes));
    for t = 1:T
        fig = figure('Visible', 'off', 'Color', 'w');
        hold on; grid on; box on;
        for im = 1:numel(modes)
            mode = modes{im};
            mask = strcmp(tbl_all.mode, mode) & tbl_all.t == t & tbl_all.solver_status == 0 & ~isnan(tbl_all.P0_pred) & ~isnan(tbl_all.Q0_pred);
            if sum(mask) < 1, continue; end
            P = tbl_all.P0_pred(mask);
            Q = tbl_all.Q0_pred(mask);
            if sum(mask) >= 3
                [Ps, Qs] = sort_boundary_points(P, Q);
                area = polygon_area_pq(Ps, Qs);
                plot([Ps; Ps(1)], [Qs; Qs(1)], '-', 'LineWidth', 1.7, 'Color', colors(im, :), ...
                    'DisplayName', sprintf('%s area=%.3f', mode, area));
            end
            plot(P, Q, 'o', 'Color', colors(im, :), 'MarkerFaceColor', colors(im, :), 'HandleVisibility', 'off');
        end
        xlabel('P_0 export to upstream grid (MW)');
        ylabel('Q_0 export to upstream grid (MVAr)');
        title(sprintf('IEEE-9 P-Q flex comparison, t=%02d, scan=%s', t, strrep(params.SCAN_MODE, '_', '\_')));
        legend('Location', 'best');
        fname = fullfile(params.RESULT_DIR, sprintf('pq_flex_compare_t%02d.png', t));
        save_png_compat(fig, fname);
        close(fig);
        fprintf('Saved comparison figure: %s\n', fname); drawnow;
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
