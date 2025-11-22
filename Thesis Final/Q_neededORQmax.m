% Put this before building constraints
dt = 3600;% hourly
load('JanFeb_Prices.mat');
Price = TT_full.price_values/1000; % convert MWh to kWh
load('weather_hourly.mat');
%weather_T=zeros(length(DateTime),1);
T_out = Temperature_clean;
n  = min(length(Price), length(T_out));
Price =Price(1:n);
%Price = price_values(1:n);
T_out = T_out(1:n);
T_min = T_set(1);
T_max= T_set(2);
Q_needed = zeros(n, num_f);
load("RCinput.mat");
Q_max= 300;
for f = 1:num_f
    Q_needed(:,f) = max(0, (T_min - T_out(1:n))/R_zone(f) - Q_int);  % W
end
Qmax_req = max(Q_needed, [], 'all');   % W needed to always hold T_min
fprintf('Q_max required to hold T_min everywhere: %.0f W\n', Qmax_req);

if Q_max < Qmax_req
    Q_max = Qmax_req;                 % <-- smallest change that guarantees feasibility
    fprintf('Q_max increased to %.0f W to make model feasible.\n', Q_max);
end

% Highest T_min you can enforce with current Q_max
Tmin_feasible = inf;
for f = 1:num_f
    Tmin_feasible = min(Tmin_feasible, min(T_out(1:n) + R_zone(f)*(Q_max + Q_int)));
end
if T_min > Tmin_feasible
    fprintf('T_min=%.2f°C is too high; reducing to %.2f°C for feasibility.\n', T_min, Tmin_feasible);
    T_min = Tmin_feasible;
end
