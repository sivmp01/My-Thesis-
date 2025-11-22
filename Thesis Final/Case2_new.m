% %% Simulation parameters-----
% % loaded from Thermal model Tin_Day.m
% %yalmip('clear')
n = length(T_out);% time steps
dt= 3600;%  Time steps in s (10 minutes)
delta_t  = dt/3600;%  time steps in h
T_initial = 20;  T_min = 18;  T_max = 23;% °C
COP= 2.5;
Q_int_2= 500; % W / floor can be changed 
Q_max_2 = 25000; % W / floor -- 7500 is assummed value for 15 floors 1000m2 

if length(Price) == n-1
    Price(n) = Price(end); 
end

%% Storage parameters 2 here represents case 2
S_max2 = 200;%  SOC kWh
eta_c2 = 0.95;
eta_d2  = 0.95;
Discharge_limit2 = S_max2/num_f; % kWh per floor per step
P_Disfloor_max2  = Discharge_limit2/delta_t; % kW
S0_2= 0;  % kWh
P_charge_max2 = 7.5; % kW

%% Decision variables
Q_2  = sdpvar(n+1, num_f);% W
T_2= sdpvar(n+1, num_f);% °C
slack_2= sdpvar(n+1, num_f);% °C
S_2 = sdpvar(n+1, 1);% kWh
P_ch_2= sdpvar(n+1, 1);% kW
P_dis_f2 = sdpvar(n+1, num_f);% kW per floor % f2 represents case 2 

%% Constraints 
con2 = []; % pre allocate constraint matrix 

% initial conditions
con2 = [con2; T_2(1,:) == T_initial];%initial temperature
con2 = [con2; S_2(1)   == S0_2];% initial storage level

% dynamics & power limits
for t = 1:n
    tot_dis2 = sum(P_dis_f2(t,:)) * delta_t; % kWh
    
    for f = 1:num_f
        a = dt/(C_zone(f)*R_zone(f));
        b = dt/C_zone(f);
        Q_storage2 = P_dis_f2(t,f) * 1000;% kW to W  Floorwise
        
        con2 = [con2;
            T_2(t+1,f) == T_2(t,f) + a*(T_out(t)-T_2(t,f)) ...
                        + b*(Q_2(t,f) + Q_int_2 + Q_storage2)]; % thermal constraint
        
        con2 = [con2;  0 <= P_dis_f2(t,f) <= P_Disfloor_max2]; % discharge constraint
        con2 = [con2; P_dis_f2(t,f)*delta_t <= Discharge_limit2];
    end
    
    con2 = [con2;
        S_2(t+1) == S_2(t) + eta_c2*P_ch_2(t)*delta_t - tot_dis2/eta_d2;
        tot_dis2 <= S_2(t)]; %storage SoC constraint
end

% comfort & capacity limits
con2 = [con2; T_2 + slack_2 >= T_min;   T_2 <= T_max; slack_2 >= 0]; % indoor temperature ;
con2 = [con2; 0 <= S_2 <= S_max2];% SoC
con2 = [con2; 0 <= Q_2 <= Q_max_2]; % Heating
con2 = [con2; 0 <= P_ch_2 <= P_charge_max2]; % chanrging
for t = 1:n
    con2 = [con2; P_ch_2(t)*delta_t <= S_max2 - S_2(t)];% overcharge constraint 
end

%% ---------------- objective ----------------
penalty = 1e6;
E_HVAC_e = (Q_2(1:n,:)*dt)./(3.6e6*COP);% kWh electrical energy 
E_ch_e   = P_ch_2(1:n)*delta_t;% kWh heatpump stores so, electrical
cost_2     = Price'*sum(E_HVAC_e,2)  + Price'*E_ch_e;% €
objective_2 = cost_2 + penalty*sum(slack_2(:));

%% Solving the Optimisation problem 
ops2= sdpsettings('solver','gurobi','verbose',1);
diagnostics_2 = optimize(con2, objective_2, ops2);
%check for problems 
if diagnostics_2.problem ~= 0
    disp('Optimization issue:'), disp(diagnostics_2.info), check(con2);
end

%% -- Check point for where the optimisation failed
if diagnostics_2.problem ~= 0
    disp('Optimization failed:');
    disp(diagnostics_2.info)
end

%% Results Extraction
Q2_values   = value(Q_2);
T2_values   = value(T_2);
S2_values   = value(S_2);
Pch2_values = value(P_ch_2);
Pdis2_values= value(P_dis_f2);
obj2_value  = value(objective_2); 
slack2_val  = value(slack_2);
%comfort violations 
if all(slack2_val < 1e-6)
    disp('Model is feasible! No comfort violations. You did it SIVA !!! CASE 2 is done');
else
    disp('Comfort violation occurred at:');
    viol2=find(slack2_val > 1e-3);
end
% Quick whole-horizon electricity & cost from raw values
rowsQ  = size(Q2_values,1);
rowsT  = size(T2_values,1);
rowsS  = size(S2_values,1);
rowsPd = size(Pdis2_values,1);
rowsPc = size(Pch2_values,1);
rowsP  = numel(Price);
hasDT2 = (exist('DateTime','var') && ~isempty(DateTime));
rowsDT = hasDT2 * numel(DateTime) + ~hasDT2 * inf;

% Common safe horizon (never exceed any array)
N2 = min([n, rowsQ, rowsT, rowsS, rowsPd, rowsPc, rowsP, rowsDT]);

% Slice everything to N2
Q2      = Q2_values(1:N2,:);        % W per floor
T2      = T2_values(1:N2,:);        % °C per floor
S2      = S2_values(1:N2,1);        % kWh
Pch2    = Pch2_values(1:N2,1);      % kW electric charging
Pdis2   = Pdis2_values(1:N2,:);     % kW - thermal discharge per floor
Price2  = Price(1:N2);   % €/kWh

% Time vector (prefer DateTime if long enough; else numeric hours)
dtC2     = dt;  % [s/step], e.g., 600
delta_t2 = delta_t;  % [h/step], e.g., 1/6
if hasDT2 && numel(DateTime) >= N2
    time2 = DateTime(1:N2);
else
    time2 = (0:N2-1)' * (dtC2/3600);    % hours from start (numeric)
end
fprintf('Objective value = €%.2f ', obj2_value);
% Per-step electricity (kWh): HVAC + charging
COPC2       = COP;
Q2_sum_W    = sum(Q2,2);% W total across floors
E2_hvac_kWh = Q2_sum_W * dtC2/(3.6e6*COPC2);     % kWh per step (electric)
E2_ch_kWh   = Pch2 * delta_t2;% kWh per step (electric)
E2_tot_kWh  = E2_hvac_kWh + E2_ch_kWh;% kWh per step

% Costs
cost2_hour  = E2_tot_kWh .* Price2;% € per step
cost2_cum   = cumsum(cost2_hour); % €
cost2_total = sum(cost2_hour); % €
E2_total_kWh= sum(E2_tot_kWh); % kWh

% Average powers over actual simulated time
hours_total2 = N2 * delta_t2; % h
avg_Pele_kW2 = E2_total_kWh / hours_total2;% kW avg electric
avg_Pth_kW2  = avg_Pele_kW2 * COPC2; % kW avg thermal
peak_Qth_kW2 = max(Q2_sum_W)/1000; % kW_th peak thermal

fprintf('Case 2: Total electricity = %.1f kWh, Total cost = €%.2f\n', ...
        E2_total_kWh, cost2_total);
fprintf('Case 2: Avg Pel = %.2f kW_e, Avg Pth = %.2f kW_th, Peak Pth = %.2f kW_th\n', ...
        avg_Pele_kW2, avg_Pth_kW2, peak_Qth_kW2);


%% SAVE OUTPUTS
num_fC2 = num_f;
save('Case2_out.mat', ...
    'Q2','T2','S2','Pch2','Pdis2','Price2','time2','dtC2','delta_t2','COPC2','num_fC2', ...
    'E2_hvac_kWh','E2_ch_kWh','E2_tot_kWh','E2_total_kWh', ...
    'cost2_hour','cost2_cum','cost2_total', ...
    'avg_Pele_kW2','avg_Pth_kW2','peak_Qth_kW2');

%% Plot file 
PlotCASE2;




