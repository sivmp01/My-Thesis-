load("RCinput.mat");
load('thermal.mat');
load('JanFeb_Prices.mat');
load('plotCase1.mat');
load('Case2_out.mat');

t= DateTime_clean(1:n);

n     = size(T2_values,1) - 1;   % time steps used for decisions
num_f = size(T2_values,2);
t_h   = (1:n) * delta_t;         % hours
%% Run Baseline function 
[Tb, Tb_avg, Qb, Qb_avg, Eb_tot_kWh, Eb_zone_kWh, E_J_e] = ...
    Building_RC_Thermostat( T_out(1:n), dt,n*dt,num_f, C_zone, R_total, R_zone,T_initial, T_set, Q_int, COP);

% Floors to plot
ground_f = 1;                 
middle_f = ceil(num_f/2);     
top_f    = num_f;             

figure('Name','CASE 2 - Ground, Middle, and Top Floor Temps','Color','w');
% Ground floor
subplot(3,1,1);
plot(t, T2_values(1:n, ground_f), 'LineWidth', 1.2); hold on;
ylabel('Indoor Temperature °C'); grid on;
%legend('Ground floor','Middle floor','Top floor', 'Location','best');
title(sprintf('Ground Floor (F%d)', ground_f));

ax = gca;
ax.XLim = [t(1) t(end)];
ax.XTick = t(1:24:end);              % one tick per day
datetick('x','dd-mmm','keepticks');
ax.XGrid = 'on'; ax.YGrid = 'on';
ax.GridAlpha = 0.3; ax.LineWidth = 0.6;

% Middle floor
subplot(3,1,2);
plot(t, T2_values(1:n, middle_f), 'LineWidth', 1.2); hold on;
ylabel('Indoor Temperature °C');grid on;
title(sprintf('Middle Floor (F%d)', middle_f));

ax = gca;
ax.XLim = [t(1) t(end)];
ax.XTick = t(1:24:end);              % one tick per day
datetick('x','dd-mmm','keepticks');
ax.XGrid = 'on'; ax.YGrid = 'on';
ax.GridAlpha = 0.3; ax.LineWidth = 0.6;


% Top floor
subplot(3,1,3);
plot(t, T2_values(1:n, top_f), 'LineWidth', 1.2); hold on;
ylabel('Indoor Temperature°C'); 
xlabel('Time (h)'); grid on;
title(sprintf('Top Floor (F%d)', top_f));

ax = gca;
ax.XLim = [t(1) t(end)];
ax.XTick = t(1:24:end);              % one tick per day
datetick('x','dd-mmm','keepticks');
ax.XGrid = 'on'; ax.YGrid = 'on';
ax.GridAlpha = 0.3; ax.LineWidth = 0.6;
sgtitle('Indoor Temperatures – CASE 2 (Ground, Middle, Top Floors)');

%% 2) HVAC heating profile (electrical kW)
P_HVAC_kW  = sum(Q2_values(1:n,:),2)/(1000*COP);% total HVAC but for energy *dt/ 3.6e3
P_dis_kW   = sum(Pdis2_values(1:n,:),2);% total discharge
P_ch_kW    = Pch2_values(1:n); % charging power (if available)
cost_P_HVAC_kW= P_HVAC_kW.* Price(1:n);
cost_storage_kW= P_ch_kW.*Price(1:n);

%fprintf('HVAC cost = €%.2f | Storage cost = €%.2f ', cost_P_HVAC_kW, cost_storage_kW );
figure('Name','CASE 2 - HVAC, Charging & Discharging','Color','w');

% --- Subplot 1: HVAC only ---
subplot(2,1,1);
plot(t, P_HVAC_kW, '-b', 'LineWidth', 1.2); hold on;
xlabel(''); ylabel('kW'); grid on;
title('HVAC Electrical Power – CASE 2');
legend('HVAC (kW)','Location','best');
ax = gca;
%te=t(end);
ax.YLim = [-5, max(P_HVAC_kW)*1.15]; % add some headroom
ax.XLim = [t(1) t(end)];
ax.XTick = t(1:24:end);  % one tick per day (if datetime)
datetick('x','dd-mmm','keepticks');
ax.GridAlpha = 0.3; ax.LineWidth = 0.8;

% --- Subplot 2: Storage charging & discharging ---
subplot(2,1,2);
plot(t, P_ch_kW, '-b', 'LineWidth', 1.2); hold on;      % charging (green)
plot(t, P_dis_kW, '--r', 'LineWidth', 1.2);             % discharging (red)
xlabel('Time (h)'); ylabel('kW');
legend('Charge (kW)','Discharge (kW)','Location','best');
title('Thermal Storage Power – CASE 2');
grid on;
ax = gca; %te=t(end);
ax.YLim = [-5, max(max(P_ch_kW), max(P_dis_kW))*1.15];
ax.XLim = [t(1) t(end)];
ax.XTick = t(1:24:end);
datetick('x','dd-mmm','keepticks');
ax.GridAlpha = 0.3; ax.LineWidth = 0.8;

% --- Global title and font ---
sgtitle('CASE 2 – HVAC and Storage Power Profiles');
set(findall(gcf,'-property','FontSize'),'FontSize',14);

%% SoC
figure('Name','CASE 2 - Storage SoC','Color','w'); clf;
plot(t, S2_values(1:n), 'LineWidth', 0.9);
xlabel('Time (h)'); ylabel('kWh');
title('Storage SoC – CASE 2');
grid on;

ax = gca;
ax.YLim = [-10, max(S2_values)*1.1]; % small headroom (use *1.1, not *10)
ax.XLim = [t(1) t(end)];
ax.XTick = t(1:24:end); % one tick per day
datetick('x','dd-mmm','keepticks');
ax.GridAlpha = 0.3;
ax.LineWidth = 0.8;

