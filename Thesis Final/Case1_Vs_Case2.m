%% Compare_Case1_Case2.m

% --- Load outputs from all cases except case 3 ---
c1 = load('Case1_out.mat');
c2 = load('Case2_out.mat');
b = load('Baseline_out.mat');
Baseline = b.Baseline;


% --- Trim to a common horizon just in case ---
N = min(length(c1.Price1), length(c2.Price2));
Q1 = c1.Q1(1:N,:);     T1 = c1.T1(1:N,:);     Price1 = c1.Price1(1:N);
Q2 = c2.Q2(1:N,:);     T2 = c2.T2(1:N,:);     Price2 = c2.Price2(1:N);

% Time vectors
if isdatetime(c1.time1)
    time = c1.time1(1:N);
elseif isdatetime(c2.time2)
    time = c2.time2(1:N);
else
    time = (0:N-1)';  % hour index
end

% Electricity per step 
E1_hvac_kWh = c1.E1_hvac_kWh(1:N);  % HVAC only
E2_hvac_kWh = c2.E2_hvac_kWh(1:N);
E2_ch_kWh   = c2.E2_ch_kWh(1:N);
E2_tot_kWh  = c2.E2_tot_kWh(1:N);

% Cumulative costs
cost1_cum = c1.cost1_cum(1:N);
cost2_cum = c2.cost2_cum(1:N);

%% ---------- PLOTS (Whole Horizon) ----------
% 1) Indoor temperature
T1_avg = mean(T1,2);
T2_avg = mean(T2,2);

figure('Name','Indoor Temperature - Whole Horizon','Color','w');
plot(time, T1_avg, 'LineWidth', 1); hold on; grid on;
plot(time, T2_avg, 'LineWidth', 1);
yline(18,'--'); yline(23,'--');
xlabel('Time'); ylabel('Indoor Temp [°C]');
title('Indoor Temperature (Avg over floors) — Case 1 vs Case 2');
legend('Case 1','Case 2','T_{min}','T_{max}','Location','best');

% 2) HVAC heating power profile
Q1_total = sum(Q1,2);  % W
Q2_total = sum(Q2,2);  % W

figure('Name','HVAC Heating Power - Whole Horizon','Color','w');
plot(time, Q1_total/1000, 'LineWidth', 1); hold on; grid on;
plot(time, Q2_total/1000, 'LineWidth', 1);
xlabel('Time'); ylabel('HVAC Thermal Power [kW_{th}]');
title('HVAC Thermal Power — Case 1 vs Case 2');
legend('Case 1','Case 2','Location','best');

% 3) Total electric consumption (per step)
figure('Name','Electric Consumption - Whole Horizon','Color','w');
plot(time, E1_hvac_kWh, 'LineWidth', 1); hold on; grid on;
plot(time, E2_tot_kWh, 'LineWidth', 1);
xlabel('Time'); ylabel('Electricity [kWh per step]');
title('Electricity Use — Case 1 vs Case 2');
legend('Case 1 (HVAC only)','Case 2 (HVAC + Charging)','Location','best');

% 4) Storage SoC (Case 2 only, whole horizon)
if hasSoC
    figure('Name','Storage SoC - Whole Horizon','Color','w');
    plot(time, S2, 'LineWidth', 1); grid on;
    xlabel('Time'); ylabel('SoC [kWh]');
    title('Thermal Storage SoC — Case 2');
end

% 5) Cumulative cost (whole horizon) [this part you already had]

% 5) Cumulative cost (whole horizon)
figure('Name','Cumulative Cost - Whole Horizon','Color','w');
plot(time, cost1_cum, 'LineWidth', 1); hold on; grid on;
plot(time, cost2_cum, 'LineWidth', 1);
xlabel('Time'); ylabel('Cumulative Cost [€]');
title('Cumulative Cost — Case 1 vs Case 2 (Whole Horizon)');
legend(sprintf('Case 1 (Total=€%.2f)', c1.cost1_total), ...
       sprintf('Case 2 (Total=€%.2f)', c2.cost2_total), 'Location','best');

%% SUMMARY TABLE - Case 1 and Case 2 
% Peak HVAC - convert W to kW
peak_Q1_kW = max(Q1_total)/1000;
peak_Q2_kW = max(Q2_total)/1000;

% Total energy (kWh electric)
E1_total_kWh = sum(E1_hvac_kWh);
E2_total_kWh = sum(E2_tot_kWh);

% Summary = table( ...
%     E1_total_kWh,          c1.cost1_total,     peak_Q1_kW, ...
%     E2_total_kWh,          c2.cost2_total,     peak_Q2_kW, ...
%     'VariableNames', {'E_kWh_Case1','Cost_EUR_Case1','Peak_kWth_Case1', ...
%                       'E_kWh_Case2','Cost_EUR_Case2','Peak_kWth_Case2'});
disp('Summary--->');
% disp(Summary);
Summary = table( ...
    Baseline.energy_kWh, Baseline.cost_total, Baseline.peak_kWth, ...
    E1_total_kWh,        c1.cost1_total,     peak_Q1_kW, ...
    E2_total_kWh,        c2.cost2_total,     peak_Q2_kW, ...
    'VariableNames', {'E_kWh_BL','Cost_EUR_BL','Peak_kWth_BL', ...
                      'E_kWh_Case1','Cost_EUR_Case1','Peak_kWth_Case1', ...
                      'E_kWh_Case2','Cost_EUR_Case2','Peak_kWth_Case2'});
disp(Summary);

%% --- Coldest & Warmest Day by T_out (Baseline + Case1 + Case2) ---

stepsPerDay = 24;   % adjust if not hourly
numDays = floor(N/stepsPerDay);

% --- daily averages of T_out
Tout_daily = zeros(numDays,1);
for d = 1:numDays
    idx = (d-1)*stepsPerDay+1 : d*stepsPerDay;
    Tout_daily(d) = mean(T_out(idx));
end

[~,coldDay] = min(Tout_daily);
[~,warmDay] = max(Tout_daily);

% --- Precompute daily totals for each case ---
dailyE_BL = zeros(numDays,1); dailyC_BL = zeros(numDays,1);
dailyE1   = zeros(numDays,1); dailyC1   = zeros(numDays,1);
dailyE2   = zeros(numDays,1); dailyC2   = zeros(numDays,1);

for d = 1:numDays
    idx = (d-1)*stepsPerDay+1 : d*stepsPerDay;
    dailyE_BL(d) = sum(Eb_step_kWh(idx));
    dailyC_BL(d) = sum(costb_step(idx));
    dailyE1(d)   = sum(E1_hvac_kWh(idx));
    dailyC1(d)   = sum(E1_hvac_kWh(idx).*Price1(idx));
    dailyE2(d)   = sum(E2_tot_kWh(idx));
    dailyC2(d)   = sum(E2_tot_kWh(idx).*Price2(idx));
end

% --- Helper print function ---
printDay = @(label,d) fprintf([ ...
'\n=== %s (Day %d, Hours %d–%d) ===\n' ...
'T_out range: %.1f – %.1f °C (avg %.1f °C)\n' ...
'Energy: Baseline %.1f kWh | Case1 %.1f kWh | Case2 %.1f kWh\n' ...
'Cost:   Baseline %.2f €   | Case1 %.2f €   | Case2 %.2f €\n'], ...
label, d, (d-1)*stepsPerDay, d*stepsPerDay-1, ...
min(T_out((d-1)*stepsPerDay+1:d*stepsPerDay)), ...
max(T_out((d-1)*stepsPerDay+1:d*stepsPerDay)), ...
Tout_daily(d), ...
dailyE_BL(d), dailyE1(d), dailyE2(d), ...
dailyC_BL(d), dailyC1(d), dailyC2(d));

% --- Print coldest & warmest ---
printDay('Coldest Day', coldDay);
printDay('Warmest Day', warmDay);
