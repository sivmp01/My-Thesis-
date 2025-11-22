%% Load SMHI Hourly Temperature Data

% File path (adjust if needed)
filename = 'C:\Users\sk090\MATLAB\Projects\Jan 24 till Feb 24.csv';

% Read data into table
data = readtable(filename);

% Extract columns
Date = string(data{:,1});         % Datum column
Time =string( data{:,2});         % Tid (UTC) column
Temperature = data{:,3};          % Lufttemperatur column
Quality = data{:,4};              % Kvalitet column

% Combine date and time into full timestamp
% DateTime = datetime(strcat(Date, {' '}, Time), 'InputFormat','yyyy-MM-dd HH:mm:ss');
DateTime = datetime(Date + " " + Time, 'InputFormat','yyyy-MM-dd HH:mm:ss');

% Optional: filter only good quality data
% good_idx = strcmp(Quality,'G');
% DateTime = DateTime(good_idx);
% Temperature = Temperature(good_idx);

% %% Plot to check data
% plot(DateTime, Temperature, '-')
% xlabel('Date Time')
% ylabel('Air Temperature (°C)')
% title('Hourly Temperature from SMHI')
% grid on

DateTime_expected = (DateTime(1):hours(1):DateTime(end))';

% Compare length:
fprintf('Expected length: %d\n', length(DateTime_expected));
fprintf('Actual length: %d\n', length(DateTime));

% Find missing timestamps
missing_idx = ~ismember(DateTime_expected, DateTime);
missing_times = DateTime_expected(missing_idx);

% Display missing timestamps
disp('Missing timestamps:');
disp(missing_times);

% Create timetable for easy filling
TT = timetable(DateTime, Temperature);

% Re-time to full hourly grid
TT_full = retime(TT, DateTime_expected, 'linear');  % forward-fill

% Check what was filled
filled_idx = ismissing(TT_full.Temperature);
fprintf('Number of missing filled: %d\n', sum(filled_idx));

% Now you have fully clean data
DateTime_clean = TT_full.DateTime;
Temperature_clean = TT_full.Temperature;

%%checking again
missing_idx = ~ismember(DateTime_expected, DateTime_clean);
missing_times_check = DateTime_expected(missing_idx);

disp('Missing timestamps:');
disp(missing_times);
weather_table_clean = table(DateTime_clean, Temperature_clean);

%% OPTIONAL: Save as MAT file for your model
save('weather_hourly.mat','DateTime_clean','Temperature_clean');
%% Extracting coldest & warmest days + hour indices
day = dateshift(DateTime_clean,'start','day');              % day bucket
T = table(DateTime_clean, day, Temperature_clean);

% Daily stats
G = groupsummary(T,'day',{'min','max','mean'},'Temperature_clean');

% Coldest day (by daily minimum temperature)
[~, iCold]   = min(G.min_Temperature_clean);
coldest_day  = G.day(iCold);
coldest_temp = G.min_Temperature_clean(iCold);
coldest_idx  = (day == coldest_day) & (Temperature_clean == coldest_temp);
coldest_hours_idx = find(coldest_idx);   % numeric indices (e.g., 200th hour)

% Hottest/Warmest day (by daily maximum temperature)
[~, iHot]    = max(G.max_Temperature_clean);
hottest_day  = G.day(iHot);
hottest_temp = G.max_Temperature_clean(iHot);
hottest_idx  = (day == hottest_day) & (Temperature_clean == hottest_temp);
hottest_hours_idx = find(hottest_idx);   % numeric indices (e.g., 200th hour)

% Print nicely
fprintf('\nColdest day: %s | Tmin = %.2f °C at hour(s): %s\n', ...
    datestr(coldest_day,'yyyy-mm-dd'), coldest_temp, mat2str(coldest_hours_idx));

fprintf('Warmest day: %s | Tmax = %.2f °C at hour(s): %s\n', ...
    datestr(hottest_day,'yyyy-mm-dd'), hottest_temp, mat2str(hottest_hours_idx));
