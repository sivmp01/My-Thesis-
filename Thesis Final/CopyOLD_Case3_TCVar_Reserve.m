%Reserves can be promised from 1000 till 1500kW else it is 0 no less than
%100kWh 
%% Does it work ? No as of
% 23/09 
%% ---- Input DATA ---
%clear; clc;
%Tin_Day;
%yalmip('clear')

% --- Load data ---
%Price_data_filteration;
load('JanFeb_Prices.mat');
Price = TT_full.price_values / 1000;% €/kWh

%LoadWeather;
load('weather_hourly.mat');
T_out = Temperature_clean;% °C

%Tin_Day;
load('RCinput.mat');
load('thermal.mat');

% --- Parameters ---
n         = length(T_out);          % timesteps
dt        = 600; delta_t = dt/3600; % time [h]
T_initial = 20; T_min = 18; T_max = 23;
COP       = 2.5;
Q_int_3   = 500; 
Q_max_3 = 30000;% each floor 2kWh max (1000m2)from Q_max file 

S_max3  = 150;  % kWh 
eta_c3 = 0.95; eta_d3 = 0.95;
P_charge_max3   = 5;
Discharge_limit3= S_max3 / num_f;
P_floor_max3    = Discharge_limit3 / delta_t;
S0_3  = 0;

% --- Market participation---
B  = 8;      % buildings
Ns  = 3;      % scenarios
pi_s  = [0.3 0.5 0.2]; % ones(Ns,1) / Ns ------------------------> changed before executing so havent run the code yet !!
% 50% of chance that the actual forecast could happen and the cold case
% could be 30% and warm case could be 20% 
pen_short = 300;    % €/MWh
cap_price = 30;      % €/MW/h
SlackPenalty = 1e6;
%% Presolve sanity check
Pel1_building_cap = (num_f * Q_max_3) / (COP*1000) + P_charge_max3;  % kW
Pel_pool_cap = B * Pel1_building_cap;                           % kW
fprintf('Pel1 cap per building ~ %.1f kW; Pooled capacity  ~ %.1f kW\n', ...
        Pel1_building_cap, Pel_pool_cap);
%% Scenario decision variables
T_out_s = cell(Ns,1);
Price_s = cell(Ns,1);
Q_int_s = cell(Ns,1);

for w = 1:Ns
    T_out_s{w} = T_out + 2.0*(w - ceil(Ns/2)); % change in T_out for different Scenarios  % ±0.5 °C
    Price_s{w} = Price .* (1 + 0.2*(w - ceil(Ns/2))); % change in Price for different scenarios  % ±5%
    Q_int_s{w} = Q_int_3 * ones(n,1);
end

%%  Stage 1 Variables & Constraints
Q1 = sdpvar(n+1, num_f);  
T1 = sdpvar(n+1, num_f);
S1 = sdpvar(n+1, 1); 
Pch1 = sdpvar(n+1, 1);
Pdis1_f = sdpvar(n+1, num_f); 
slk1 = sdpvar(n+1, num_f);
Res_RA = sdpvar(n,1); % reserve commitment
z_res = binvar(n,1);  % binary: 1 if reserve committed at time t
Res_min =1000;%1MW
Res_max=3000;
conStg1=[];
conStg1 = [conStg1;T1(1,:) == T_initial, S1(1) == S0_3];

for t = 1:n
    total_dis1 = sum(Pdis1_f(t,:)) * delta_t;
    for f = 1:num_f
        a = dt/(C_zone(f)*R_zone(f));
        b = dt/C_zone(f);
        Qstor = Pdis1_f(t,f)*1000;
        conStg1 = [conStg1;
            T1(t+1,f) == T1(t,f) + a*(T_out(t)-T1(t,f)) + b*(Q1(t,f) + Q_int_3 + Qstor)];
            conStg1 = [conStg1;0 <= Pdis1_f(t,f) <= P_floor_max3];
            %conStg1 = [conStg1;Pdis1_f(t,f)*delta_t <= Discharge_limit3];
    end
    conStg1 = [conStg1;S1(t+1) == S1(t) + eta_c3*Pch1(t)*delta_t - total_dis1/eta_d3];
    conStg1 = [conStg1 ;total_dis1 <= S1(t)];
end

% --- Electric power (per building, decisions 1..n) ---
Pel1 = (sum(Q1(1:n,:),2)/(COP*1000)) + Pch1(1:n);   % kW

% --- Comfort & slack (apply once, not inside loops) ---
conStg1 = [conStg1; ...
  T1(2:end,:) + slk1(2:end,:) >= T_min, ...
  T1(2:end,:) <= T_max, ...
  0 <= slk1 <= (T_max - T_min), ...
  0 <= S1 <= S_max3, 0 <= Q1 <= Q_max_3, 0 <= Pch1 <= P_charge_max3];

% --- Storage charging bound (vectorized once over all t) ---
conStg1 = [conStg1; Pch1(1:n)*delta_t <= S_max3 - S1(1:n)];

% --- Reserve gating & caps (apply ONCE, outside loops) ---
conStg1 = [conStg1; ...
  0 <= Res_RA <= Res_max, ...         % base bounds
  Res_RA >= Res_min * z_res, ...      % min 1 MW if committed
  Res_RA <= Res_max * z_res, ...      % 0 if not committed
  Res_RA <= B * Pel1];                % physics cap: can’t reduce > you consume

%conStg1 = [conStg1;T1 + slk1 >= T_min, T1 <= T_max, slk1 >= 0];
% conStg1 = [conStg1; T1(2:end,:) + slk1(2:end,:) >= T_min, T1(2:end,:) <= T_max, 0 <= slk1 <= (T_max - T_min)];
% conStg1 = [conStg1; 0 <= S1 <= S_max3, 0 <= Q1 <= Q_max_3, 0 <= Pch1 <= P_charge_max3];
% 
% % for t=1:n
% %     conStg1 = [conStg1; Pch1(t)*delta_t <= S_max3 - S1(t)];
% %     conStg1 = [conStg1; T1(2:end,:) + slk1(2:end,:) >= T_min, T1(2:end,:) <= T_max, 0 <= slk1 <= (T_max - T_min)];
% %     %conStg1=[conStg1;Res_RA >=Res_min *z_res];
% %     %conStg1=[conStg1;Res_RA <=Res_max *z_res];
% % end
% % %Electric power 
% conStg1 = [conStg1; Res_RA <= B * Pel1];  % fleet can’t reduce > it consumes
% 
% Pel1 = (sum(Q1(1:n,:),2)/(COP*1000)) + Pch1(1:n);
%% ----Stage 2 Variables & Constraints----
Q2 = cell(Ns,1); 
T2 = cell(Ns,1);
S2 = cell(Ns,1);
Pch2 = cell(Ns,1); Pdis2_f = cell(Ns,1);
slk2 = cell(Ns,1);
%r = cell(Ns,1);
sshort = cell(Ns,1); nu = cell(Ns,1);
conRec = cell(Ns,1);

for w = 1:Ns
    Q2{w} = sdpvar(n+1, num_f);      T2{w} = sdpvar(n+1, num_f);
    S2{w} = sdpvar(n+1, 1);          Pch2{w} = sdpvar(n+1, 1);
    Pdis2_f{w} = sdpvar(n+1, num_f); slk2{w} = sdpvar(n+1, num_f);
    %r{w} = sdpvar(n,1);              
    sshort{w} = sdpvar(n,1);% short fall between Reserve in stage 1 and stage 2 
    nu{w} = sdpvar(1);% Tail losses in TCVaR

    conRec{w} = [conRec{w};T2{w}(1,:) == T_initial];
    conRec{w}=[conRec{w};S2{w}(1) == S1(end)];

    for t = 1:n
        total_dis2 = sum(Pdis2_f{w}(t,:)) * delta_t;
        for f = 1:num_f
            a = dt/(C_zone(f)*R_zone(f));
            b = dt/C_zone(f);
            Qstor2 = Pdis2_f{w}(t,f)*1000;
            conRec{w} = [conRec{w};T2{w}(t+1,f) == T2{w}(t,f) + a*(T_out_s{w}(t)-T2{w}(t,f)) ...
                             + b*(Q2{w}(t,f) + Q_int_s{w}(t) + Qstor2)];
            conRec{w} = [conRec{w};   0 <= Pdis2_f{w}(t,f) <= P_floor_max3];
            %conRec{w} = [conRec{w};Pdis2_f{w}(t,f)*delta_t <=
            %Discharge_limit3]; %repeated the before line 
        end
        conRec{w} = [conRec{w};S2{w}(t+1) == S2{w}(t) + eta_c3*Pch2{w}(t)*delta_t - total_dis2/eta_d3];
        conRec{w} = [conRec{w};total_dis2 <= S2{w}(t)];
    end

    %conRec{w} = [conRec{w};T2{w} + slk2{w} >= T_min, T2{w} <= T_max, slk2{w} >= 0];
    conRec{w} = [conRec{w};T2{w}(2:end,:) + slk2{w}(2:end,:) >= T_min, T2{w}(2:end,:) <= T_max,0 <= slk2{w} <= (T_max - T_min)];
    conRec{w} = [conRec{w}; 0 <= S2{w} <= S_max3, 0 <= Q2{w} <= Q_max_3, 0 <= Pch2{w} <= P_charge_max3];

    for t=1:n
        conRec{w} = [conRec{w}; Pch2{w}(t)*delta_t <= S_max3 - S2{w}(t)];
    end

    %conRec{w} = [conRec{w}; S2{w}(end) == S0_2];%for day wise but n is not
    %144 which is 1 day ; hence it is multiday horizon!! we do it for 60
    %days hence change it in initialisation also if youre switching to day
    %wise  change  it to S2{w}(1) == S0_2;


    Pel2 = (sum(Q2{w}(1:n,:),2)/(COP*1000)) + Pch2{w}(1:n);

    conRec{w} = [conRec{w};
        Pel2 <= Pel1 + 1e-6;
        %r{w} >= 0; 
        sshort{w} >= 0; %shortfall limits, cannot be negative 
        %r{w} <= B*(Pel1 - Pel2); % difference of the two stage Power is Commited as reserve in Stage 2 after aggregating
        %r{w} +
        sshort{w} >= Res_RA-B*(Pel1 - Pel2)];
end

%% Objective with TCVaR
alpha = 0.8; % considering there will be 20% of worst cases 
z = sdpvar(1); % VaR threshold
TCVaR_obj = z;
conRA = [];%constraints of Risk aversion 

for w = 1:Ns
    E_hvac = (Q2{w}(1:n,:)*dt) / (3.6e6 * COP);
    E_ch   = Pch2{w}(1:n)*delta_t;
    E_tot  = B * (sum(E_hvac,2) + E_ch);

    cost_energy = sum(E_tot .* Price_s{w}(1:n));
    cost_penalty = pen_short * sum(sshort{w} * delta_t / 1000);
    scenario_cost = cost_energy + cost_penalty;

    conRA = [conRA;
        nu{w} >= scenario_cost - z;
        nu{w} >= 0];

    TCVaR_obj = TCVaR_obj + (1/(1 - alpha)) * pi_s(w) * sum(nu{w});
end

% Add slack and subtract revenue
TCVaR_obj = TCVaR_obj - cap_price * sum(Res_RA) / 1000;
%TCVaR_obj = TCVaR_obj + SlackPenalty * (sum(slk1(:)) + sum(cellfun(@(s) sum(s(:)), slk2)));
sum_slk2 = 0; for w=1:Ns, sum_slk2 = sum_slk2 + sum(slk2{w}(:)); end
TCVaR_obj = TCVaR_obj + SlackPenalty * (sum(slk1(:)) + sum_slk2);

%% SOLVE -- 2 stage Stochastic Optimisation model with Risk avaersion TCVaR
conAll = [conStg1; conRA];
for w = 1:Ns
    conAll = [conAll; conRec{w}];
end

%ops = sdpsettings('solver','gurobi','verbose',1);
% Safety assertions for shapes
% assert(iscolumn(Price) && numel(Price)==n, 'Price must be n×1 column.');
% %assert(iscolumn(hod)   && numel(hod)==n  , 'hod must be n×1 column.');
% assert(iscolumn(y)     && numel(y)==n    , 'y must be n×1 column.');
% assert(iscolumn(Res_RA)&& numel(Res_RA)==n, 'Res_RA must be n×1 column.');

ops_3t = sdpsettings('solver','gurobi','verbose',1, ...
  'gurobi.MIPGap',       0.005, ...    % 0.1% -> 0.5% often 2x faster, same result quality here
  'gurobi.TimeLimit',    2000,  ...
  'gurobi.MIPFocus',     1,     ...    % focus on bound
  'gurobi.Presolve',     2,     ...
  'gurobi.Heuristics',   0.1,   ...    % reduce a bit to save time
  'gurobi.Cuts',         2,     ...
  'gurobi.Threads',      max(1, feature('numcores')-1), ... % leave 1 core free
  'gurobi.Aggregate',    1, ...
  'gurobi.PreSparsify',  1, ...
  'gurobi.Symmetry',     2);

% ops = sdpsettings('solver','gurobi','verbose',1, ...
%   'gurobi.MIPGap',       0.001, ...
%   'gurobi.TimeLimit',    2000,  ...
%   'gurobi.MIPFocus',     1,     ...
%   'gurobi.Presolve',     2,     ...
%   'gurobi.Heuristics',   0.2,   ...
%   'gurobi.NumericFocus', 1,     ...
%   'gurobi.Cuts',         2);
% ops.gurobi.LogFile = 'case3_gurobi.log';
 ops_3t.gurobi.NodefileStart = 0.5;
diagnostics = optimize(conAll, TCVaR_obj, ops_3t);
solver_time = diagnostics.solvertime
if diagnostics.problem ~= 0
    disp('[ERROR] Optimization failed:')
    disp(diagnostics.info)
else
    fprintf('\n T-CVaR Two-Stage Optimization SUCCESS! You did it SIVA ! Case 3 with CVaR is done. \n');
    R_opt = value(Res_RA);
    fprintf('Max committed reserve: %.2f kW\n', max(R_opt));
    fprintf('Avg committed reserve: %.2f kW\n', mean(R_opt));
end

%save("Case3_TCVAR_Reserve.mat");

% if diagnostics.problem == 0
%    % Stage 1
%     T1v   = value(T1);            % (n+1) x num_f
%     S1v   = value(S1);            % (n+1) x 1
%     Pch1v = value(Pch1);          % (n+1) x 1
%     Q1v   = value(Q1);            % (n+1) x num_f
%     Pel1v = (sum(Q1v(1:n,:),2)/(COP*1000)) + Pch1v(1:n);  % kW
% 
%     % Stage 2 (mid scenario)
%     T2v_mid   = value(T2{mid});
%     S2v_mid   = value(S2{mid});
%     Pch2v_mid = value(Pch2{mid});
%     Q2v_mid   = value(Q2{mid});
%     Pel2v_mid = (sum(Q2v_mid(1:n,:),2)/(COP*1000)) + Pch2v_mid(1:n);
% 
%     % Other
%     R_opt   = value(Res_RA);
%     Price_v = Price(:);
%     t_h     = (0:n) * delta_t;
%     t_h_n   = (1:n) * delta_t;
%     [~, num_f_detect] = size(T2v_mid);
%     f_idx   = unique([1, max(1, round(num_f_detect/2)), num_f_detect]); % ground/middle/top
% 
%     save('Case3_basic.mat', ...
%          't_h','t_h_n','R_opt','Price_v','S2v_mid','Pel2v_mid', ...
%          'T2v_mid','f_idx','T_min','T_max');
% end

if diagnostics.problem == 0
    mid = ceil(Ns/2);

    % ---- Stage 1 (per-building) ----
    T1v   = value(T1);                        % (n+1) x num_f
    S1v   = value(S1);                        % (n+1) x 1
    Pch1v = value(Pch1);                      % (n+1) x 1
    Q1v   = value(Q1);                        % (n+1) x num_f
    Pel1v = (sum(Q1v(1:n,:),2)/(COP*1000)) + Pch1v(1:n);   % kW

    % ---- Stage 2 (per scenario, per-building) ----
    Pel2_all = zeros(n, Ns);                  % kW
    S2_all   = zeros(n+1, Ns);                % kWh
    sshort_all = zeros(n, Ns);                % kW shortfall-equivalent (aggregated)
    for w = 1:Ns
        S2w   = value(S2{w});
        Pchw  = value(Pch2{w});
        Q2w   = value(Q2{w});
        Pel2w = (sum(Q2w(1:n,:),2)/(COP*1000)) + Pchw(1:n); % kW

        S2_all(:,w)   = S2w;
        Pel2_all(:,w) = Pel2w;
        sshort_all(:,w) = value(sshort{w});
    end

    % ---- Prices for each scenario (n x Ns) ----
    Price_all = zeros(n, Ns);
    for w = 1:Ns
        Price_all(:,w) = Price_s{w}(1:n);
    end

    % ---- Other items ----
    R_opt   = value(Res_RA);                  % kW (aggregate commitment)
    t_h     = (0:n) * delta_t;                % hours, states
    t_h_n   = (1:n) * delta_t;                % hours, decisions

    % (optional) ground/middle/top indices if you ever need temps later
    [~, num_f_detect] = size(T1v);
    f_idx = unique([1, max(1, round(num_f_detect/2)), num_f_detect]);

    save('Case3_dump.mat', ...
        'n','Ns','pi_s','delta_t','B','cap_price','pen_short', ...
        't_h','t_h_n','R_opt', ...
        'S1v','Pel1v', ...
        'S2_all','Pel2_all','Price_all','sshort_all', ...
        'T_min','T_max','f_idx');
    fprintf('Saved all results for fast plotting -> Case3_dump.mat\n');
end
