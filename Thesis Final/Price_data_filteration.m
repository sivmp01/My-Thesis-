%% this file optimises the Price data 
%clc; clear;

% Step 1: Load Excel Data
filename = "C:\Users\sk090\MATLAB\Projects\GUI_ENERGY_PRICES_202401010000-202501010000.csv";  % Replace with actual filename
data = readtable(filename);

% Step 2: Extract and combine Var1 (date) and Var2 (hour)
date_strings = string(data.Var1);      % e.g., "01/01/2024"
hour_strings = string(data.Var2);      % e.g., "23:00:00"
datetime_strings = date_strings + " " + hour_strings;

% Step 3: Convert to datetime
timestamps = datetime(datetime_strings, 'InputFormat', 'dd/MM/yyyy HH:mm:ss');

% Step 4: Extract prices from Var5
% Var5 format: '00:00:00,29.56' → extract the price after the comma
raw_prices = string(data.Var5);
price_parts = split(raw_prices, ',');
price_values = str2double(price_parts(:,2));  % This gives you numeric price values

% Step 5: Create timetable
TT = timetable(timestamps, price_values);

% Step 6: Generate expected hourly range for Jan & Feb 2024
expectedTimes = (datetime(2024,1,1,0,0,0):hours(1):datetime(2024,2,29,23,0,0))';

% Step 7: Retime to fill missing timestamps (use previous value)
TT_full = retime(TT, expectedTimes, 'previous');

% Step 8: Check for missing values
missingIdx = isnan(TT_full.price_values);
if any(missingIdx)
    fprintf('Missing timestamps:\n');
    disp(expectedTimes(missingIdx));
else
    fprintf('No missing data found . 1440 entries found for Jan–Feb.\n');
end
% Step 10: Extract price vector
price_vector = TT_full.price_values/1000; % MWh to kWh 
assert(length(price_vector) == 1440, ' Price vector is not 1440 entries long.');
fprintf("Cleaned price vector is ready for use in your simulation.\n");

% Step 9: Save outputs
writetimetable(TT_full, 'Cleaned_JanFeb_Prices.xlsx');
save('JanFeb_Prices.mat', 'TT_full','price_vector');

