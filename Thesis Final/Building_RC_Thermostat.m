function [T_in_t,T_in_t_avg,Q_HVAC_t,Q_HVAC_avg,E_HVAC_kWh,E_HVAC_zonewise_kWh,E_HVAC_J] = Building_RC_Thermostat(weather_T, dt, tot_t,num_f,C_zone, R_total,R_zone, T_initial,T_set,Q_int,COP)

% Time discretization
n_hours = floor(tot_t / dt);
% indoor temp per floor
T_in_t = zeros(n_hours, num_f); %%initialising the array
%T_in at the start of the data 
T_in_t(1,:) = T_initial;
% energy use in joules
E_HVAC_J = zeros(n_hours,num_f); %Pre-allocate 
% heating power per floor
Q_HVAC_t= zeros(n_hours, num_f);
delta = 0.7; % 0.7 degree deadband   
% C_zone = TMF_vector .* C_air_vector; % capacitance per floor

% Time loop
for i = 1:n_hours
    T_out = weather_T(i);

    for f = 1:num_f
        % Thermostat control for each zone
        if T_in_t(i, f) < T_set(1)
            % Heating required
            Q_HVAC_t(i, f) = ((T_set(2) - T_in_t(i, f)) / (dt / C_zone(f))) ...
                         - ((T_out - T_in_t(i, f)) / R_total) - Q_int;
            % Prevent negative heating (no cooling)
            Q_HVAC_t(i, f) = max(Q_HVAC_t(i, f), 0);

        elseif T_in_t(i, f) > T_set(2)
            Q_HVAC_t(i, f) = 0; % heating off
        else
            % Within deadband → hold previous state
            if i == 1
                Q_HVAC_t(i, f) = 0;
            else
                Q_HVAC_t(i, f) = Q_HVAC_t(i - 1, f);
            end
        end

        % Euler integration to update temperature
        dTdt = (T_out - T_in_t(i, f)) / R_zone(f) + Q_HVAC_t(i, f) + Q_int;
        T_in_t(i + 1, f) = T_in_t(i, f) + (dt / C_zone(f)) * dTdt;

        % Energy used
        E_HVAC_J(i, f) = Q_HVAC_t(i, f) * dt/COP;
    end
end
% T_in_tavg = zeros(n_hours,1);                % preallocate
% for i = 1:n_hours
%     T_in_tavg(i) = mean(T_in_t(i, :));       % average across floors at time i
% end
%% --average value of the floors
T_in_t_avg=mean(T_in_t,2);
%Q_HVAC_tavg = zeros(n_hours,1);                % preallocate
% for i = 1:n_hours
%     Q_HVAC_tavg(i) = mean(Q_HVAC(i, :));       % average across floors at time i
% end
%% average value of the floors 
Q_HVAC_avg = mean(Q_HVAC_t,2);
%disp(E_HVAC_t_J)
% Total HVAC energy per time step [kWh] across all floors
E_HVAC_total_kWh = sum(E_HVAC_J, 2) / (3600 * 1000);
% Convert energy from J to kWh
E_HVAC_zonewise_kWh = sum(E_HVAC_J, 1) / (3600 * 1000);  % [1 × num_f]

%% Total HVAC energy over all time and floors [kWh]
E_HVAC_kWh = sum(E_HVAC_total_kWh);  % scalar
end

