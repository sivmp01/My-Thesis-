% PARAMETERS
floor_area = 100;  % m²
lighting_density = 10;  % W/m²
lighting_power = lighting_density * floor_area;

num_people = 3;
sensible_heat_per_person = 75;  % W
equipment_power = 500;  % W

% TIME (24h hourly simulation)
dt = 1; % 1 hour
hours = 0:dt:23;

% SCHEDULES (simplified as 0 or 1)
lighting_schedule = double((hours >= 6) & (hours <= 22));
occupancy_schedule = double((hours >= 17) | (hours <= 8));
equipment_schedule = double((hours >= 8) & (hours <= 22));

% INTERNAL GAINS CALCULATION
Q_lighting = lighting_power * lighting_schedule;
Q_occupancy = num_people * sensible_heat_per_person * occupancy_schedule;
Q_equipment = equipment_power * equipment_schedule;

Q_internal = Q_lighting + Q_occupancy + Q_equipment;

% Plot to visualize
% plot(hours, Q_internal, '-o')
% xlabel('Hour of Day')
% ylabel('Internal Gains (W)')
% title('Internal Heat Gains Schedule')
% grid on
