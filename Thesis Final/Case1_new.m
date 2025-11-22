 
%% ------Simulation parameters-----
% loaded from Thermal model Tin_Day.m
yalmip('clear')
% ---price---
% in euro/MWh 
%% Price_data_filteration;
load('JanFeb_Prices.mat');
Price = TT_full.price_values/1000; % convert MWh to kWh
%% --weather data--
% LoadWeather;
load('weather_hourly.mat');
weather_T=zeros(length(DateTime),1);
T_out = Temperature_clean;
load('RCinput.mat');
%Tin_Day;
load('thermal.mat');
n = numel(T_out);% instead of length(weather_T) for more accurate counting
%n = length(weather_T);  % 1440 time steps 
dt = 3600;     % in seconds %%check with input_file 
tot_t = 24*60*dt;   % Simulate for 60 day
T_initial = 20;
T_min= 18;
T_max=23;
%load('Q_max.mat'); % can be done through this from Q_needed file  or manually add in the Q_max
%Q_max - maximum Heating Power 
Q_max1=25000;

if length(Price) == n-1
    Price(n) = Price(end);  % Extend to match length
end
%% decision variables

%heating input at each hour in W
Q_1 = sdpvar(n+1,num_f,'full'); 

%indoor temperature at each hour in C
T_1 = sdpvar(n+1,num_f,'full');
%% constraints
constraints=[];
constraints=[constraints,T_1(1,:)==T_initial];

%floor wise 
for t = 1:n-1
    for f = 1:num_f
        a = dt / (C_zone(f) * R_zone(f));   % precomputed constant to avoid non linearity
        b = dt / C_zone(f); % precomputed constant to avoid non linearity
        constraints = [constraints;T_1(t+1,f) == T_1(t,f) + a * (T_out(t) - T_1(t,f)) + b * (Q_1(t,f) + Q_int)];
    end
end
  
slack=sdpvar(n+1,num_f); % to override the constraint 
constraints = [constraints; T_1 + slack >= T_min];
constraints = [constraints; T_1 <= T_max + 1e-6];
constraints = [constraints; slack >= 0];

constraints = [constraints; 0 <= Q_1];
constraints = [constraints; Q_1 <= Q_max1];
penalty =1e6; % penalize if the constraint is overriden

E_kWh = ((Q_1(1:n,:) * dt)/(3.6e6*COP));% Convert Q from W to kWh: (Q in W * 1 hour) / 1000 = kWh electricity 
cost = sum(sum( E_kWh.* repmat(Price, 1, num_f)));
objective = cost + penalty * sum(sum(slack(1:n,:)));
%objective =0; % to check if the model is unbounded 

%% solve the problem
ops=sdpsettings('solver','gurobi','verbose',1,'debug',1);
diagnostics = optimize(constraints, objective, ops);
slack_val = value(slack);  % slack is n×1 vector
% Check where the optimization failed
if diagnostics.problem ~= 0
    disp('Optimization failed:');
    disp(diagnostics.info)
elseif  diagnostics.problem ==0
    disp("Your Model is feasible. You did it SIVA !!! CASE 1 is done")
end

if all(slack_val < 1e-6)
    disp('No comfort violations. ');
else
    disp('Comfort violation occurred at:');
    viol=find(slack_val > 1e-3);
end

%%  Case 1 numeric results
% aligned everything before save 

nC1     = n;% rewritten for results 
dtC1    = dt;% seconds per step 
COPC1   = COP;
num_fC1 = num_f;

% Pull numeric values from sdpvars
Q1_values = value(Q_1);
T1_values = value(T_1);
obj_value = value(objective); 

% Available lengths in each series
rowsQ  = size(Q1_values,1);
rowsT  = size(T1_values,1);
rowsP  = numel(Price);
hasDT  = (exist('DateTime','var') && ~isempty(DateTime));
rowsDT = hasDT * numel(DateTime) + ~hasDT * inf;  % inf if no DateTime

% Common safe horizon
N1 = min([nC1, rowsQ, rowsT, rowsP, rowsDT]);

% Slice to horizon
Q1     = Q1_values(1:N1,:);% [W] per floor
T1     = T1_values(1:N1,:);% [°C] per floor
Price1 = Price(1:N1);% [€/kWh]

% Time vector - used numeric hours for better navigation in data;
if hasDT && numel(DateTime) >= N1
    time1 = DateTime(1:N1);
else
    time1 = (0:N1-1)' * (dtC1/3600);  % hours from start
end

%% Per step electricity (kWh) used by HVAC only
Q1_sum        = sum(Q1,2);% [W] total across floors
E1_hvac_kWh   = Q1_sum * dtC1/(3.6e6*COPC1);% [kWh] per step (electric)
cost1_hour    = E1_hvac_kWh .* Price1; % [€] per step
cost1_cum     = cumsum(cost1_hour);% [€] cumulative
cost1_total   = sum(cost1_hour); % [€] total horizon
E_kWh1_total  = sum(E1_hvac_kWh);% [kWh] total horizon
mean_price    = mean(Price1);% [€/kWh] over horizon

% Averages over actual simulation time to get precise plots 
hours_total   = N1 * (dtC1/3600);
avg_Pele_kW   = E_kWh1_total / hours_total;% [kW_e]
avg_Pth_kW    = avg_Pele_kW * COPC1; % [kW_th]

fprintf('Total electricity = %.1f kWh, Total cost = €%.2f, Avg price = %.4f €/kWh\n', ...
        E_kWh1_total, cost1_total, mean_price);
fprintf('Avg electric power = %.2f kW_e, Avg thermal power = %.2f kW_th\n', ...
        avg_Pele_kW, avg_Pth_kW);

save('Case1_out.mat','Q1','T1','Price1','time1','dtC1','COPC1','num_fC1', ...
     'E1_hvac_kWh','cost1_hour','cost1_cum','cost1_total', ...
     'E_kWh1_total','avg_Pele_kW','avg_Pth_kW');

%% --- Plotting ----
PlotCASE1;
