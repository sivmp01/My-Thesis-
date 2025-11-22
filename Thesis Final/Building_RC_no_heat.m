function [T_in_noheat,T_in_nh_avg, Q_HVAC_req,E_HVAC_req_kWh,E_HVAC_totreq] = Building_RC_no_heat(weather_T, dt, tot_t, num_f, R_zone, C_zone, T_initial, Q_int,COP)
    steps = floor(tot_t/dt) +1;  % +1 to include initial time
    n_hours = floor(tot_t / dt);
    T_in_noheat = zeros(n_hours, num_f);
    T_in_noheat(1, :) = T_initial(:).';
    N = length(weather_T);
    T_out = weather_T(:);
    % ---- Time evolution (no HVAC heating, only internal gains) ----
    for i = 1:n_hours
        for f = 1:num_f
        %T_out(i) = weather_T(i);  % Outdoor temp at current step
        % Heat flow (simplified): only one zone used here
        dTdt = (T_out(i) - T_in_noheat(i,f)) / R_zone(f) + Q_int;
        T_in_noheat(i + 1,f) = T_in_noheat(i,f) + (dt / C_zone(f)) * dTdt;
        end
    end
    T_target =21;
    for i = 1:n_hours
        for f = 1:num_f
            delta_T = T_target - T_in_noheat(i,f);
            loss_to_outside = (T_out(i) - T_in_noheat(i,f)) / R_zone(f);
            required_gain = (C_zone(f)/dt) * delta_T;
    
            Q_HVAC_req(i,f) = required_gain - loss_to_outside - Q_int;
            E_HVAC_req(i, f) = Q_HVAC_req(i, f) * dt/COP;
            % Update temperature with only internal gains (your original model)
            % dTdt1 = loss_to_outside + Q_int;
            % T_in(i + 1,f) = T_in(i,f) + (dt / C_zone(f)) * dTdt1;
        end
    end
    %disp(E_HVAC_t_J)
    T_in_nh_avg = zeros(n_hours,1);% preallocate
for i = 1:n_hours
    T_in_nh_avg(i) = mean(T_in_noheat(i, :));% average across floors at time i
end
% Total HVAC energy per time step [kWh] across all floors
E_HVAC_totreq = sum(E_HVAC_req, 2) / (3600 * 1000);

%% Total HVAC energy over all time and floors [kWh]
E_HVAC_req_kWh = sum(E_HVAC_totreq);  % scalar
end

