%doesn't work 

%if you are 


%% ------Simulation parameters-----
% loaded from Thermal model Tin_Day.m
yalmip('clear')
%% ---price---
% in euro/MWh 
Price_data_filteration;
load('JanFeb_Prices.mat');
Price = TT_full.price_values/1000; % convert MWh to kWh
%% --weather data--
LoadWeather;
load('weather_hourly.mat');
weather_T=zeros(length(DateTime),1);
T_out = Temperature_clean;
load('RCinput.mat');
%Tin_Day;
load('thermal.mat');
%n = length(weather_T);  % 1440 time steps 
n = numel(T_out);   % instead of length(weather_T)
dt=3600;% 60 minutes in seconds
delta_t = dt/3600; %time step in hours 
tot_t = 24*60*dt; % Simulate for 60 day in seconds 
T_initial = 20;T_min=18;T_max=23;
Q_max_1=ceil(Q_max*1e-3)*1e3;
%Q_max = 25000; % HVAC heating limit in Watts found from Q_neededmax file
if length(Price) == n-1
    Price(n) = Price(end);  % Extend to match length
end
%% ----decision variables----
%heating input at each hour in W
Q_1 = sdpvar(n+1,num_f,'full'); 
%Q_int=0;
%indoor temperature at each hour in C
T_1 = sdpvar(n+1,num_f,'full');
% % symmetric slacks (both sides) to guarantee feasibility
% s_lo = sdpvar(n, num_f, 'full');  % °C, relax lower bound
% s_hi = sdpvar(n, num_f, 'full');  % °C, relax upper bound
%constraints
constraints_1=[];
constraints_1=[constraints_1,T_1(1,:)==T_initial];

%floor wise 
for t = 1:n-1
    for f = 1:num_f
        a = dt / (C_zone(f) * R_zone(f));   % precomputed constant to avoid non linearity
        b = dt / C_zone(f);                 % precomputed constant to avoid non linearity
        constraints_1 = [constraints_1;T_1(t+1,f) == T_1(t,f) + a * (T_out(t) - T_1(t,f)) + b * (Q_1(t,f) + Q_int)];
        %constraints = [constraints;T_1(t+1,f) == T_1(t,f) + a * (T_out(t) - T_1(t,f)) + b * (Q_1(t,f) + Q_int)];
    end
end
% for t = 1:n-1
%     for f = 1:num_f
%         b = dt / C_zone(f);  % [seconds / (J/°C)] 
%         % Heat flow rate [W]
%         dTdt_1 = (T_out(t) - T_1(t,f))/R_zone(f) + Q_1(t,f) + Q_int;
%         constraints_1 = [constraints_1; T_1(t+1,f) == T_1(t,f) + b * dTdt_1];
%     end
% end
  
slack=sdpvar(n+1,num_f); % to override the constraint 
constraints_1 = [constraints_1; T_1 + slack >= T_min];
constraints_1 = [constraints_1; T_1 <= T_max + 1e-6];
constraints_1 = [constraints_1; slack >= 0];

constraints_1 = [constraints_1; 0 <= Q_1];
constraints_1 = [constraints_1; Q_1 <= Q_max_1];
penalty =1e6; % penalize if the constraint is overriden

E_kWh = ((Q_1(1:n,:)*dt)/(3.6e6*COP));% Convert Q from W to kWh: (Q in W * 1 hour) / 1000 = kWh electricity 
cost = sum(sum(E_kWh.* repmat(Price, 1, num_f)));
%cost = sum(sum(E_kWh .* repmat(Price(1:size(E_kWh,1)), 1, num_f)));
objective = cost + penalty * sum(sum(slack(1:n,:)));
%objective =0; % to check if the model is unbounded 
%% solve the problem
ops=sdpsettings('solver','gurobi','verbose',2);
diagnostics_1 = optimize(constraints_1, objective, ops);
slack_val_1 = value(slack);  % slack is n×1 vector
%% -- Check point for where the optimisation failed
if diagnostics_1.problem ~= 0
    disp('Optimization failed:');
    disp(diagnostics_1.info)
elseif  diagnostics_1.problem ==0
    disp("Your Model is feasible. You did it SIVA !!! CASE 1 is done")
end

if all(slack_val_1 < 1e-6)
    disp('No comfort violations. ');
else
    disp('Comfort violation occurred at:');
    viol=find(slack_val_1 > 1e-3);
end

%%  Case 1 numeric results
% aligned everything before save 

nC1     = n;% rewritten for results 
dtC1    = delta_t;% seconds per step 
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
Q1     = Q1_values(1:N1,:);           % [W] per floor
T1     = T1_values(1:N1,:);           % [°C] per floor
Price1 = Price(1:N1);                 % [€/kWh]

% Time vector - used numeric hours for better navigation in data;
if hasDT && numel(DateTime) >= N1
    time1 = DateTime(1:N1);
else
    time1 = (0:N1-1)' * (dtC1/3600);  % hours from start
end

% Per step electricity (kWh) used by HVAC only
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

% save('Case1_out.mat', ...
%     'Q1','T1','Price1','time1','dtC1','COPC1','num_fC1', ...
%     'E1_hvac_kWh','cost1_hour','cost1_cum','cost1_total', ...
%     'E_kWh1_total','avg_Pele_kW','avg_Pth_kW');


%%
fprintf('Max slack used = %.4f°C\n', max(slack_val_1(:)));
fprintf('Total slack penalty = €%.2f\n', penalty * sum(sum(slack_val_1(1:n,:))));
fprintf('Optimal Q_1: min=%.1f, max=%.1f, mean=%.1f W\n', ...
        min(Q1_values(:)), max(Q1_values(:)), mean(Q1_values(:)));

%edit("Case2_new.m");
%% --- Plotting ----
%PlotCASE1;
