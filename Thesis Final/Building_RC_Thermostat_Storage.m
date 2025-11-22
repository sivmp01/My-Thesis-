function [T_in_t_s, T_in_t_savg, Q_HVAC_ts, Q_HVAC_tsavg, ...
          E_HVACwithTs_kWh, E_HVAC_ts_zonewise_kWh, ...
          E_HVAC_ts_J, SOC_vector, Q_storage_in, Q_storage_out] = ...
    Building_RC_Thermostat_Storage(weather_T, dt, tot_t, num_f, ...
                                   C_zone, R_zone, ...
                                   T_initial, T_set, Q_int, ...
                                   COP, storage_J, ...
                                   SOC_min, SOC_max, ...
                                   charge_rate, discharge_rate)
% Initialisation
n_steps = floor (tot_t/dt);
T_in_t_s= zeros(n_steps, num_f); %input temperature during thermostat and storage 
T_in_t_s(1,:)= T_initial;
Q_HVAC_ts =zeros(n_steps,num_f); %W
E_HVAC_ts_J = zeros(n_steps,num_f); %J

% Soc value at each hour // there is only one storage in the whole building
% but the usage of each floor is propotionised.Thus SOC is (n_steps,1)

SOC_vector=zeros(n_steps,1);  %Wh
% initially the soc is in Joules then we convert it to W as per Q_hvac 
% which is done in the next equation as this is a preallocate;
SOC_vector(1)=0.5* storage_J/3.6e3;  %Wh 
% at the start of the simulation the storage is assumed to be 50%
% this SOC_vector is in Wh; not Ws so /3.6e3 and not dt ; 1 Wh =3600J;
% Power(W) = Energy(J) /Time(s); so if you want to conver it to Ws then *
% 3.6e3/dt

SOC_max = SOC_max / 3600;               % Wh
SOC_min = SOC_min / 3600;               % Wh

Q_storage_in = zeros(n_steps,1);
Q_storage_out= zeros(n_steps,1);

for i = 1:n_steps-1
    T_out = weather_T(i);
    for f = 1:num_f
        if T_in_t_s(i,f) < T_set(1)
            Q_HVAC_ts(i,f)=((T_set(2)-T_in_t_s(i,f))/(dt/C_zone(f))) - ((T_out-T_in_t_s(i,f))/R_zone(f))- Q_int;

            % if storage is required 
            Q_storageneeded = max(Q_HVAC_ts(i,f),0); %W
            E_availabile_Wh=SOC_vector(i)-SOC_min; %Wh
            Q_reqdischarge=min(discharge_rate,E_availabile_Wh*3600/dt); %% W ; Power (W) =Energy_available(Wh)*3600/dt(s)
            Q_discharged = min(Q_storageneeded,Q_reqdischarge); % W
            Q_storage_out(i)=Q_storage_out(i)+Q_discharged; % W
            SOC_vector(i+1)=SOC_vector(i)-Q_discharged*dt/3600; % W to Ws to Wh

            % remaining heat from HVAC 
            Q_HVAC_ts(i,f) = Q_storageneeded - Q_discharged; % W
            E_HVAC_ts_J(i,f)=Q_HVAC_ts(i,f)*dt/COP; %Electrical energy in J

        elseif T_in_t_s(i,f)>T_set(2)
            Q_HVAC_ts(i,f)=0;
            SOC_vector(i+1) = SOC_vector(i);
        else 
            Q_HVAC_ts(i,f) = Q_HVAC_ts(i-1,f);
            SOC_vector(i+1)=SOC_vector(i);
        end

        %Charge storage if spce available 
        % charge possible is Charge_p
        charge_p= (SOC_max) - SOC_vector(i+1); % Wh
        charge_p_W =charge_p*3600/dt; % wh to W 
        if charge_p >0
            Q_charge = min(charge_rate, charge_p_W);  %  W
            %Q_charge= max(0,min(charge_rate,charge_p_W));
            % how much energy can be stored at this time step max(0,..) is to only have non negative result
            
            Q_storage_in(i)=Q_charge ;%W
            SOC_vector(i+1) = SOC_vector(i+1) + Q_charge *dt/3600; %Wh
            E_HVAC_ts_J(i,f)= E_HVAC_ts_J(i,f)+Q_charge*dt/COP; % J
        end
        dTdt=(T_out - T_in_t_s(i,f))/R_zone(f) + Q_HVAC_ts(i,f) +Q_int;
        T_in_t_s(i+1,f)=T_in_t_s(i,f)+(dt/C_zone(f))*dTdt;
    end
end

T_in_t_savg=mean(T_in_t_s,2);
Q_HVAC_tsavg =mean(Q_HVAC_ts,2);
E_HVAC_totalwithstorage_kWh =sum(E_HVAC_ts_J,2)/3.6e6; % J to kWh
E_HVAC_ts_zonewise_kWh =sum(E_HVAC_ts_J,1)/3.6e6;
E_HVACwithTs_kWh = sum(E_HVAC_totalwithstorage_kWh);

% figure (1);
% time = (1:n_steps) * dt / 3600;
% 
% yyaxis left
% plot(time, SOC_vector/1000 , 'b-');  % Wh
% 
% ylabel('SOC (kWh)');
% yyaxis right
% plot(time, Q_storage_out / 1000, 'r--');  % W to kW
% 
% ylabel('Storage Discharge (kW)');
% xlabel('Time (hours)');
% title('Thermal Storage Operation');
% legend('SOC', 'Q_{out}');
% grid on;
% 
% figure (2);
% plot(time, cumsum(E_HVAC_ts_J)/3.6e6);
% xlabel('Time (hours)');
% ylabel('Cumulative Energy per Zone (kWh)');
% legend("Zone 1", "Zone 2"," ...","Zone N");
% title('Cumulative HVAC Energy Consumption per Zone');
% grid on;

end
