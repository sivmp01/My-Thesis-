% Quick diagnostics: are we preheating before high prices?

% Load saved outputs
c1 = load('Case1_out.mat');   % T1 [N×F], E1_hvac_kWh, Price1
c2 = load('Case2_out.mat');   % T2 [N×F], E2_* , Price2, S2 (optional)

Tmin = 18; Tmax = 23;    % set to your comfort limits if different

% Avg indoor temp across floors
T1_avg = mean(c1.T1, 2);
T2_avg = mean(c2.T2, 2);

% How often are we hugging the bounds?
near = @(x,thr,eps) mean(x > thr - eps);
epsC = 0.2;   % within 0.2°C of the bound counts as "near"
frac_T1_nearTop = near(T1_avg, Tmax, epsC);
frac_T2_nearTop = near(T2_avg, Tmax, epsC);
frac_T1_nearBot = mean(T1_avg < Tmin + epsC);
frac_T2_nearBot = mean(T2_avg < Tmin + epsC);

fprintf('Case1: near T_max %.1f%%, near T_min %.1f%%\n', 100*frac_T1_nearTop, 100*frac_T1_nearBot);
fprintf('Case2: near T_max %.1f%%, near T_min %.1f%%\n', 100*frac_T2_nearTop, 100*frac_T2_nearBot);

% Lead-lag check: does higher T precede high prices?
% Correlate T(t) with Price(t+lead)
leadH = 3;  % hours lead (adjust)
m1 = round(3600 / c1.dtC1);
m2 = round(3600 / c2.dtC2);
k1 = leadH * m1;  k2 = leadH * m2;

corr1 = corr(T1_avg(1:end-k1), c1.Price1(1+k1:end));
corr2 = corr(T2_avg(1:end-k2), c2.Price2(1+k2:end));
fprintf('Lead-%dh corr(T, future Price): Case1=%.3f, Case2=%.3f\n', leadH, corr1, corr2);

% Visual for one day (pick day 10)
day = 10;
steps1 = round(86400 / c1.dtC1);
steps2 = round(86400 / c2.dtC2);
i1 = (day-1)*steps1 + 1; j1 = min(i1+steps1-1, numel(T1_avg));
i2 = (day-1)*steps2 + 1; j2 = min(i2+steps2-1, numel(T2_avg));

t1 = (0:(j1-i1))' * (c1.dtC1/3600);
t2 = (0:(j2-i2))' * (c2.dtC2/3600);

figure('Color','w','Name','One-day temp vs price');
yyaxis left;  plot(t1, T1_avg(i1:j1), 'LineWidth', 1); hold on; grid on;
yline(Tmin,'--'); yline(Tmax,'--');
ylabel('Indoor Temp [°C]'); xlabel('Hour of day');
yyaxis right; plot(t1, c1.Price1(i1:j1), 'LineWidth', 1);
ylabel('Price [€/kWh]');
title('Case 1: Preheating (Temp) vs Price (same day)');

figure('Color','w','Name','One-day temp vs price (Case 2)');
yyaxis left;  plot(t2, T2_avg(i2:j2), 'LineWidth', 1); hold on; grid on;
yline(Tmin,'--'); yline(Tmax,'--');
ylabel('Indoor Temp [°C]'); xlabel('Hour of day');
yyaxis right; plot(t2, c2.Price2(i2:j2), 'LineWidth', 1);
ylabel('Price [€/kWh]');
title('Case 2: Preheating (Temp) vs Price (same day)');
