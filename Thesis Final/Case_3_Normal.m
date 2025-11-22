
%% ------Simulation parameters-----
% loaded from Thermal model Tin_Day.m
yalmip('clear')
%% Given Data 
%% ---price---
Price_data_filteration;
load('JanFeb_Prices.mat');
Price = TT_full.price_values/1000;% €/kWh

%% --weather data--
LoadWeather;
load('weather_hourly.mat');
T_out = Temperature_clean; % °C

%% --RC model inputs--
Tin_Day;
load('RCinput.mat'); %inputs
load('thermal.mat');% output of Tin_Day.m

n        = length(T_out);% time steps
dt       = 600;%  Time steps in s (10 minutes)
delta_t  = dt/3600;%  time steps in h
T_initial = 20;  T_min = 18;  T_max = 23;% °C
COP       = 2.5;
Q_int_3n     = 500; % W / floor can be changed 
Q_max_3n     = 15000; % W / floor -- 7500 is assummed value for 15 floors 

if length(Price) == n-1
    Price(n) = Price(end); 
end


% Storage parameters
S_max3n = 200;% kWh
eta_c2n= 0.95;
eta_d2n= 0.95;
Discharge_limit3n = S_max3n/num_f;% kWh per floor per step (as in your code)
P_Disfloor_max3n = Discharge_limit3n/delta_t;% kW
S0_3n= 0; % kWh
P_charge_max3n= 7.5;% kW

%  Reserve requirements 
Bn= 8;% number of identical buildings aggregated
Nsn = 3;% number of scenarios
pi_sn = ones(Nsn,1)/Nsn; % scenario probabilities equal but can change it to random
cap_pricen= 30;% €/MW/h capacity revenue (set >0 if you model it)
pen_shortn= 300;% €/MWh shortfall penalty (tune as needed)
rng(1);% reproducible scenario perbutations // to avoid change of scenarios every run

% Presolve sanity check
Pel1_building_capn = (num_f * Q_max_3) / (COP*1000) + P_charge_max3n;  % kW
Pel_pool_capn = Bn * Pel1_building_capn;                           % kW
fprintf('Pel1 cap per building ~ %.1f kW; Pooled capacity  ~ %.1f kW\n', ...
        Pel1_building_capn, Pel_pool_capn);
%% Second Stage Data - for each scenarios there will be a matrix or vector 
T_out_sn = cell(Nsn,1);
Price_sn = cell(Nsn,1);
Q_int_sn = cell(Nsn,1);

for wn = 1:Nsn
    % Example: small perturbations; replace with your scenario generation
    T_out_sn{wn} = T_out + 0.5*(wn-ceil(Nsn/2));% ±0.5°C, etc.
    Price_sn{wn} = Price .* (1 + 0.05*(wn-ceil(Nsn/2)));% ±5%
    Q_int_sn{wn} = Q_int_3n*ones(n,1);% keep fixed for all floors 
end

%% STAGE 1 Here and Now - for only one building 
%Decision variables 
Q1n      = sdpvar(n+1, num_f);
T1n      = sdpvar(n+1, num_f);
S1n      = sdpvar(n+1, 1);
Pch1n    = sdpvar(n+1, 1);
Pdis1_fn = sdpvar(n+1, num_f);
slk1n    = sdpvar(n+1, num_f);
Resn       = sdpvar(n,1);                 % kW committed reserve (aggregate)
slack_small_resn = sdpvar(n,1);
Res_minn =1000; % 0.1 MW is minimum value required to enter flexibility market 
%Constraints stage 1
conStg1n = [];
conStg1n = [conStg1n; T1n(1,:) == T_initial; S1n(1) == S0_3n];

for t = 1:n
    tot_dis1n = sum(Pdis1_fn(t,:))*delta_t;
    for f = 1:num_f
        a = dt/(C_zone(f)*R_zone(f));
        b = dt/C_zone(f);
        Qstor = Pdis1_fn(t,f)*1000;% W
        conStg1n = [conStg1n;T1n(t+1,f) == T1n(t,f) + a*(T_out(t)-T1n(t,f)) + b*(Q1n(t,f) + Q_int_3n + Qstor)];
        conStg1n = [conStg1n; 0 <= Pdis1_fn(t,f) <= P_Disfloor_max3n];
        conStg1n = [conStg1n; Pdis1_fn(t,f)*delta_t <= Discharge_limit3n];
    end
    conStg1n = [conStg1n;S1n(t+1) == S1n(t) + eta_c2n*Pch1n(t)*delta_t - tot_dis1n/eta_d2n];
    conStg1n = [conStg1n;tot_dis1n <= S1n(t)];

    %min reserve implement
    conStg1n = [conStg1n;Resn(t) >= 0;slack_small_resn(t) >= Res_minn - Resn(t);slack_small_resn(t) >= 0];
end

conStg1n = [conStg1n; T1n+slk1n >= T_min; T1n <= T_max; slk1n >= 0];
conStg1n = [conStg1n; 0 <= S1n <= S_max3n; 0 <= Q1n <= Q_max_3n; 0 <= Pch1n <= P_charge_max3n];

for t=1:n
    conStg1n = [conStg1n; Pch1n(t)*delta_t <= S_max3n - S1n(t)]; 
end

% Baseline electric power for one building
Pel1n = (sum(Q1n(1:n,:),2)/(COP*1000)) + Pch1n(1:n);   % kW_e per step

%% Stage 2 Decison Variables  Recourse Decision Variables
% all are 3 dimensions 
Q2n = cell(Nsn,1); 
T2n = cell(Nsn,1); 
S2n = cell(Nsn,1);
Pch2n = cell(Nsn,1); 
Pdis2_fn = cell(Nsn,1); 
slk2n = cell(Nsn,1);
rn = cell(Nsn,1); 
sshortn = cell(Nsn,1);
conRecn = cell(Nsn,1);

for wn = 1:Nsn
    Q2n{wn}      = sdpvar(n+1, num_f);
    T2n{wn}      = sdpvar(n+1, num_f);
    S2n{wn}      = sdpvar(n+1, 1);
    Pch2n{wn}    = sdpvar(n+1, 1);
    Pdis2_fn{wn} = sdpvar(n+1, num_f);
    slk2n{wn}    = sdpvar(n+1, num_f);
    rn{wn}       = sdpvar(n,1);     % kW delivered reserve (aggregate, across B)
    sshortn{wn}  = sdpvar(n,1);     % kW shortfall

    conRecn{wn} = [];
    conRecn{wn} = [conRecn{wn}; T2n{wn}(1,:) == T_initial; S2n{wn}(1) == S0_3n];

    for t = 1:n
        total_dis2n = sum(Pdis2_fn{wn}(t,:))*delta_t;
        for f = 1:num_f
            a = dt/(C_zone(f)*R_zone(f));
            b = dt/C_zone(f);
            Qstor2 = Pdis2_fn{wn}(t,f)*1000;% W
            conRecn{wn} = [conRecn{wn};T2n{wn}(t+1,f) == T2n{wn}(t,f) + a*(T_out_sn{wn}(t)-T2n{wn}(t,f)) ...
                + b*(Q2n{wn}(t,f) + Q_int_sn{wn}(t) + Qstor2)];
            conRecn{wn} = [conRecn{wn}; 0 <= Pdis2_fn{wn}(t,f) <= P_Disfloor_max3n];
            conRecn{wn} = [conRecn{wn}; Pdis2_fn{wn}(t,f)*delta_t <= Discharge_limit3n];
        end
        conRecn{wn} = [conRecn{wn};
            S2n{wn}(t+1) == S2n{wn}(t) + eta_c2n*Pch2n{wn}(t)*delta_t - total_dis2n/eta_d2n;
            total_dis2n <= S2n{wn}(t)];
    end

    conRecn{wn} = [conRecn{wn}; T2n{wn}+slk2n{wn} >= T_min; T2n{wn} <= T_max; slk2n{wn} >= 0];
    conRecn{wn} = [conRecn{wn}; 0 <= S2n{wn} <= S_max3n; 0 <= Q2n{wn} <= Q_max_3n; 0 <= Pch2n{wn} <= P_charge_max3n];
    for t=1:n
        conRecn{wn} = [conRecn{wn}; Pch2n{wn}(t)*delta_t <= S_max3n - S2n{wn}(t)]; 
    end
    conRecn{wn} = [conRecn{wn}; S2n{wn}(end) == S0_3n]; % initial condition of Battery 

    % Electric power during recourse stage 
    Pel2n = (sum(Q2n{wn}(1:n,:),2)/(COP*1000)) + Pch2n{wn}(1:n);   % kW_e

    % To avoid reduction in second stage 
    conRecn{wn} = [conRecn{wn}; Pel2n <= Pel1n + 1e-9];

    % Delivered reserve and shortfall (AGGREGATED over B identical buildings)
    conRecn{wn} = [conRecn{wn}; rn{wn} >= 0; sshortn{wn} >= 0];
    conRecn{wn} = [conRecn{wn}; rn{wn} <= Bn*(Pel1n - Pel2n)];     % linear, nonnegative due to Pel2<=Pel1
    conRecn{wn} = [conRecn{wn}; rn{wn} + sshortn{wn} >= Resn];       % meet commitment or pay shortfall
end

%%  Objective: minimize expected cost - capacity revenue 
% Capacity revenue: €/MW/h * MW
CapRevn = cap_pricen* delta_t * sum(Resn)/1000;
%Pre allocate 
ExpEnergyCostn = 0;
ExpPenaltyn    = 0;
SlackPenaltyn  = 1e6*sum(slk1n(:));

for wn = 1:Nsn
    % Energy cost for B buildings under scenario w
    E_hvac_wn = (Q2n{wn}(1:n,:)*dt)/(3.6e6*COP);    % kWh single building
    E_ch_wn   =  Pch2n{wn}(1:n)*delta_t;            % kWh single building
    E_tot_wn  = Bn*( sum(E_hvac_wn,2) + E_ch_wn );   % [n×1] kWh Aggregate of B no of buildings 
    ExpEnergyCostn = ExpEnergyCostn + pi_sn(wn)*sum( E_tot_wn .* Price_sn{wn}(1:n) );

    % Shortfall penalty: pen_short [€/MWh] * (kW*hours)/1000
    ExpPenaltyn = ExpPenaltyn + pi_sn(wn) * pen_shortn * sum( sshortn{wn} * delta_t / 1000 );

    SlackPenaltyn = SlackPenaltyn + 1e6*sum(slk2n{wn}(:));
end

objRecstagen = ExpEnergyCostn + ExpPenaltyn - CapRevn + SlackPenaltyn;

%% YALMIP 2 STAGE STOCHASTIC OPTIMIZATION PROBLEM SOLVE 
conAlln = [conStg1n];
for wn = 1:Nsn
    conAlln = [conAlln; conRecn{wn}]; 
end

ops2n = sdpsettings('solver','gurobi','verbose',1);
diag2n = optimize(conAlln, objRecstagen, ops2n);
if diag2n.problem ~= 0
    disp('Two-stage optimization issue:'), disp(diag2n.info);
end
diag2n.solvertime
%% ---------- Report ----------
R_valuen = value(Resn);
% fprintf('\n[TWO-STAGE] Aggregation B = %d, scenarios = %d\n', B, Ns);
fprintf('Committed reserve R(t): peak = %.1f kW, average = %.1f kW\n', max(R_valuen), mean(R_valuen));

% Delivered reserve statistics (expected, approximate)
r_expn = zeros(n,1);
for wn = 1:Nsn
    r_expn = r_expn + pi_sn(wn)*value(rn{wn});
end
fprintf('Delivered reserve (expected): peak = %.1f kW, average = %.1f kW\n', max(r_expn), mean(r_expn));

% Optional plots
% time = (1:n)*delta_t;
% figure; plot(time, R_opt, 'LineWidth',1.2); hold on; plot(time, r_exp, 'LineWidth',1.2);
% grid on; xlabel('Time (h)'); ylabel('kW'); legend('Committed R','Expected delivered r');
% title(sprintf('Two-Stage Reserve (B=%d, Ns=%d)', B, Ns));
