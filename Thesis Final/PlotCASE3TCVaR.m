% %% Plotting CASE 3

load("RCinput.mat");
load('thermal.mat');
load('JanFeb_Prices.mat');
load('Case3_TCVaR_3600.mat');
%%  STAGE 1 
t_all = TT_full.Properties.RowTimes;   % 1440×1 datetime from the timetable
% 1440-sample series (hourly):
t_days_hn   = t_all; % for t_h_n length = 1440 (R_opt, hr_rev, R_deliv)
t_days_Pel1 = t_all;% for Pel1v (1440)
t_days_Pel2 = t_all; % for Pel2_all (1440)
t_days_1440 = t_all;  % generic 1440 vector if needed

% 1441-sample series (hourly edges):
t_days_h = [t_all; t_all(end)+hours(1)]; % for t_h length = 1441 (S1v, S2_all)

% 1) Reserve commitment (Stage 1 only)
figure('Color','w');
stairs(t_days_hn, R_opt, 'LineWidth', 1.4); grid on;
title('CASE 3 — Stage 1: Committed Reserve (Res_{RA})');
xlabel('Time (h)'); ylabel('Reserve (kW)');
ax = gca;
ax.YLim=[-10,1500];
ax.XLim = [t_days_h(1) t_days_h(end)]; % was: [t(1) t(end)]
ax.XTick = t_days_h(1:24:end);
%ax.XTickLabelRotation = 45;
ax.GridAlpha = 0.3; ax.LineWidth = 0.8;

%% 2) Storage SoC (Stage 1 only)
figure('Color','w');
plot(t_days_h, S1v, 'LineWidth', 0.9); grid on;
title('CASE 3 — Stage 1: Storage SoC');
xlabel('Time (date)'); ylabel('SoC (kWh)');
ax = gca;
ax.YLim=[-10,200];
ax.XLim = [t_days_h(1) t_days_h(end)]; % was: [t(1) t(end)]
ax.XTick = t_days_h(1:24:end);
%ax.XTickLabelRotation = 45;
ax.GridAlpha = 0.3; ax.LineWidth = 0.8;

%% 3) Electric power per building (Stage 1 only)
%t_days = DateTime_clean(1:length(Pel1v));
figure('Color','w');
plot(t_days_Pel1, Pel1v, 'LineWidth', 0.9); grid on;
title('CASE 3 — Stage 1: Electric Power (per building)');
xlabel('Time '); ylabel('Power (kW)');
ax = gca;
ax.XLim = [t_days_Pel1(1) t_days_Pel1(end)]; % was: [t(1) t(end)]
ax.XTick = t_days_Pel1(1:24:end);
%ax.XTickLabelRotation = 45;
ax.GridAlpha = 0.3; ax.LineWidth = 0.8;

%% 4) Hourly capacity revenue (Stage 1 only)
hr_rev = cap_price * (R_opt/1000);   % €/MW/h * MW = €/h
figure('Color','w');
plot(t_days_hn, hr_rev, 'LineWidth', 0.9); grid on;
title('CASE 3 — Stage 1: Hourly Capacity Revenue');
xlabel('Time'); ylabel('EUR/h');
ax = gca;
ax.XLim = [t_days_Pel1(1) t_days_Pel1(end)]; % was: [t(1) t(end)]
ax.XTick = t_days_Pel1(1:24:end);
%ax.XTickLabelRotation = 45;
ax.GridAlpha = 0.3; ax.LineWidth = 0.8;
% (optional) show total in console:
fprintf('Total capacity revenue (EUR): %.2f\n', cap_price*sum(R_opt)/1000);

%% STAGE 2 

% Helper: grid size for subplots
rows = ceil(sqrt(Ns)); cols = ceil(Ns/rows);

% Delivered reserve per scenario (kW): B * (Pel1 - Pel2), clipped at 0
R_deliv = B * max(0, Pel1v - Pel2_all);   % n x Ns

% Rebuild scenario costs (energy + shortfall penalty)
E_cost = zeros(Ns,1); PEN = zeros(Ns,1);
for w = 1:Ns
    E_cost(w) = sum( B * Pel2_all(:,w) .* delta_t .* Price_all(:,w) );
    PEN(w)    = pen_short * sum( sshort_all(:,w) .* delta_t / 1000 );
end
SC = E_cost + PEN;  % €/scenario
[SCs, idx] = sort(SC,'ascend');  % VaR for reference
cumP = cumsum(pi_s(idx));
alpha = 0.8;
VaR = SCs(find(cumP >= alpha, 1, 'first'));

%% 5) Stage 2 — Delivered reserve (subplot per scenario)
figure('Color','w');
tiledlayout(rows, cols, 'TileSpacing','compact', 'Padding','compact');
for w = 1:Ns
    nexttile; plot(t_days_hn, R_deliv(:,w), 'LineWidth', 0.8); grid on;
    title(sprintf('Scenario %d', w));
    xlabel('Time (h)'); ylabel('Reserve (kW)');
    ax = gca;
    ax.YLim=[-10,1600];
    xlim([t_days_Pel2(1) t_days_Pel2(end)]);
    ax.XTick = t_days_Pel2(1:24:end);
    xtickformat('dd-MMM');% use this with datetime
    ax.XTickLabelRotation = 45;
    ax.GridAlpha = 0.3; ax.LineWidth = 0.8;
end
sgtitle('CASE 3 — Stage 2: Delivered Reserve (per scenario)');

%% 6) Stage 2 — Storage SoC (subplot per scenario)
figure('Color','w');
tiledlayout(rows, cols, 'TileSpacing','compact', 'Padding','compact');
for w = 1:Ns
    nexttile; plot(t_days_h, S2_all(:,w), 'LineWidth', 0.8); grid on;
    title(sprintf('Scenario %d', w));
    xlabel('Time (date)'); ylabel('SoC (kWh)');
    ax = gca;
    ax.YLim=[-10,200];
    ax.XLim = [t_days_h(1) t_days_h(end)]; % was: [t(1) t(end)]
    ax.XTick = t_days_h(1:24:end);
    %ax.XTickLabelRotation = 45;
    ax.GridAlpha = 0.3; ax.LineWidth = 0.8;
end
sgtitle('CASE 3 — Stage 2: Storage SoC (per scenario)');

%% 7) Stage 2 — Electric power per building (subplot per scenario)
figure('Color','w');
tiledlayout(rows, cols, 'TileSpacing','compact', 'Padding','compact');
for w = 1:Ns
    nexttile; 
    plot(t_days_Pel2, Pel2_all(:,w), 'LineWidth', 0.8); grid on;
    title(sprintf('Scenario %d', w));
    xlabel('Time'); ylabel('Power (kW)');
    ax = gca;
    xlim([t_days_Pel2(1) t_days_Pel2(end)]);
    ax.XTick = t_days_Pel2(1:24:end);
    xtickformat('dd-MMM');        % use this with datetime
    ax.XTickLabelRotation = 45;
    ax.GridAlpha = 0.3; ax.LineWidth = 0.8;
end
sgtitle('CASE 3 — Stage 2: Electric Power (per scenario)');
