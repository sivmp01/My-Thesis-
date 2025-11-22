
%% -- READ ME !
% -- This file contains inputs required to run this model, you don't have to specifically run this file. Just run Tin_Day.m file which runs this file first! %%
%clc;
%clear all;
%% Simulation parameters
dt = 3600;  % Time step in seconds 
tot_t = 24*60*dt;   % Simulate for 60 day % seconds

%% Initial indoor temp (°C)
T_initial = 20;

%% Building parameters
%% ---------- Building geometry----------%%%%%%%%%%
%% wall to window ratio

W2W_ratio = 0.3; % from BBR complaince

air_changes_per_hour = 0.5; %ACH to calculate ventilation this is an assumption

% number of floors
num_f = 15

% Floor Area 
f_A = 1000 % m²

% Floor Height 
f_h = 3% m 
%Total height of the building is f_h * num_f

% if youre changing the num_f, f_A,f_h then Q_max in all cases should be
% changed or the prblem can become infeasible 

% Volume of total building - this is only if the whole building is modelled as one thermal zone
V=f_A*f_h*num_f; %m3 

% Area of Wall including windows - whole building
wall_A = 4 * sqrt(f_A) * f_h ;%m2 

% Area of window - whole building
window_A = W2W_ratio * wall_A; %m2 

% Area of opaque wall area - whole building
opaque_wall_area = wall_A - window_A; %m2

% area of roof and ceiling
roof_area = f_A; %m2 

%% Air Capacitance of building
cp_air = 1005;           % Specific heat of air (J/kg·K)
rho_air = 1.2;          % Air density (kg/m³)
C_air = rho_air * cp_air * V; %Air Capacitance of whole building together 
%% U values of building 

% Wall
w_U = 0.18;% wall U value  W/m²K
r_U = 0.08;% roof u value W/m²K <=0.13
win_U = 1.1;% window U value W/m²K

%% -------Material inputs-------- %%
%% wall 
% Gypsum, Wood ,mineral wool, CLT, Poly PE barrier 
d_wall = [0.0127, 0.0254, 0.20, 0.0127, 0.0000762];         % thickness of wall material m
lambda_wall = [0.9022,0.0150, 0.0346, 0.0837, 0.2885];       % conductivity of wall material W/mK
rho_wall = [800.925, 608.703, 56.06, 480.56, 929.073];      % density of wall material kg/m3
cp_wall = [1088.568, 1632.852, 837.36, 1549.12, 1884.06];  % specific heat capacity of wall material K/Kg.K

    %Calculating U value for Wall
    %Calculate R-value for each layer
    R_layers_wall = d_wall ./ lambda_wall;   % [m²·K/W]
    R_total_wall = sum(R_layers_wall);
    % U value = 1/R
    U_wall = 1 / R_total_wall; % if surface resistance are added 
    fprintf('U_wall     = %.4f W/m²·K\n', U_wall);

%% Same for floor, roof, window
% concrete, PIR, Finish

%% Floor except ground floor 
d_floor   = [0.02, 0.10, 0.15]; 
lambda_floor=[0.4,1.8,0.022];
rho_floor = [700, 2300, 35];           
cp_floor  = [1600, 880, 1400];
    %R-values for each layer
    R_floor = d_floor ./ lambda_floor;
    R_total_floor = sum(R_floor);
    %  U-value
    U_floor = 1 / R_total_floor;
    fprintf('U_floor    = %.4f W/m²·K\n', U_floor);

%% Roof 
%Membrane insulation metal deckoning
d_roof = [0.0095, 0.169291, 0.0015];
lambda_roof = [0.20, 0.022, 1771.845];
rho_roof = [1121.26, 265, 7679.91];
cp_roof = [1461.193, 837.36, 418.68];

    R_roof = d_roof ./ lambda_roof;
    U_roof=1/(sum(R_roof));
    fprintf('U_roof    = %.4f W/m²·K\n', U_roof);
%% ceiling-
% Gypsum board,Air gap, Concrete slab, Mineral Wool

d_ceiling = [0.0127, 0.10, 0.20, 0.25];
lambda_ceiling = [0.25, 0.20, 1.5, 0.0346];
rho_ceiling = [800, 1.2, 2200, 50];
cp_ceiling = [1090, 1000, 880, 837.36];
    R_ceiling = d_ceiling ./ lambda_ceiling;
    U_ceiling = 1 / sum(R_ceiling);
    fprintf('U_ceiling  = %.4f W/m²·K\n', U_ceiling);

%% Floor
d_floor_gnd     = [0.02, 0.10, 0.15];   % thickness of the ground floor material     
lambda_floor_gnd = [0.4, 1.8, 0.022];   % conductivity of the ground floor material 
rho_floor_gnd   = [700, 2300, 35];      % density of the of the ground floor material 
cp_floor_gnd    = [1600, 880, 1400];    % specific heat capacity of the ground floor material 
    R_ground=d_floor_gnd./lambda_floor_gnd;
    U_ground=1/sum(R_ground);
    fprintf('U_ground   = %.4f W/m²·K\n', U_ground);

 %% Window (double glazing with argon)
d_window = [0.004, 0.016, 0.004];
lambda_window = [1.0, 0.016, 1.0];
rho_window = [2240, 1.784, 2240];
cp_window = [837, 520, 837];
    R_window = d_window ./ lambda_window;
    U_window = 1 / sum(R_window);
    fprintf('U_window   = %.4f W/m²·K\n', U_window);
    
%% Resistance of the building
% ---- Per-Floor Area ----
    wall_A_f = wall_A / num_f;
    window_A_f = window_A / num_f;
    opaque_wall_A_f = opaque_wall_area / num_f;

    % ---- Ventilation loss (air change per hour) ----
    UA_vent = air_changes_per_hour * C_air / 3600; % whole building ventilation loss coefficient [W/K]
    UA_vent_f = UA_vent/num_f; % ventilation loss per floor
    R_zone = zeros(num_f, 1);  % Thermal resistance per zone
    % ---- Compute R_zone ----
    for f = 1:num_f
        if f == 1
            % Ground floor: includes floor loss (simplified as extra wall)
            UA = U_wall * opaque_wall_A_f + U_window * window_A_f + U_ground * f_A;
        elseif f == num_f
            % Roof floor: includes roof loss
            UA = U_wall * opaque_wall_A_f + U_window * window_A_f + U_roof * roof_area;
        else
            % Middle floors
            UA = U_wall * opaque_wall_A_f + U_window * window_A_f;
        end

        R_zone(f) = 1 / (UA + UA_vent_f);  % K/W
    end

    % ---- Use average R_total across all floors ----
    R_total = sum(R_zone);  % this is for the whole building  

%% ----- add internal gains -----
%%lets say 500W 
Q_int= 500; %per floor  % W


%% --Coefficient of Power 
COP =2.5;

%% -----Control parameters---
%[min,max]
T_set = [18 23]; % °C %% if changing here change in  no heat function also  

%% --Price ---
% %% e_price=
% price_preheat=0.2;% €/kWh %%below this it will charge
% price_peak=0.3;% €/kWh  %%above this it will discharge the thermal mass i.e the heating will be off


%% --- Thermal Storage Parameters --- %%
% let's consider the system has thermal storage, with each floor can
% utilise upto 1.5kWh
% totalStorage = 1.5 * num_f; in our base case the whole building has 15kWh
% of thermal storage 

storage_kWh=1.5*num_f;

% since most of our calculations are done in J and then converted to kWh,
% to avoid confusion, we will use storage in Joules and at the end we will
% convert it to kWh. The same heat pump charges the storage; hence COP=2.5


storage_J = storage_kWh *3.6e6;% storage in Joules 

% To avoid complete depletion of thermal storage and incase of critical
% situations, the storage can be discharged upto 20%, if required this 20%
% can be used by manual trigger 

SOC_min =0.2*storage_J; %J

% To avoid over heating and loss due to conventional heat exchange, the
% thermal storage can be heated up till 90 % only !

SOC_max=0.9*storage_J; %J

% Charging and discharging rate would be same; we assume it as 3e3 

charge_rate =3e3;  % Max charging power (W)
discharge_rate =3e3; % Max discharging power (W)

%%save it 
save('RCinput.mat',"T_set","T_initial",'-append');
% get max q_h
Q_neededORQmax;
%Price 
Price_data_filteration;
%opening Tin_Day file after execution of this file 
edit("Tin_Day.m");
