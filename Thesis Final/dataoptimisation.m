%% This file optimises the Energy prices it finds the missing index ( Daylight savings ) 
% optimises so it can be used in the project 

% File path (adjust if needed)
filename = 'C:\Users\sk090\Downloads\Energyprices 2024.csv';

% Read the full table (assuming first row is header, separator is comma or semicolon depending on file)
opts = detectImportOptions(filename);
data = readtable(filename, opts);

% Extract MTU column as strings
MTU_strings = string(data{:,1});

% Extract start time from MTU (everything before " - ")
start_strings = extractBefore(MTU_strings, ' -');

% Convert to datetime
MTU_start = datetime(start_strings, 'InputFormat', 'dd/MM/yyyy HH:mm:ss');

% Generate full expected hourly time series
start_full = min(MTU_start);
end_full = max(MTU_start);
expected_time = (start_full:hours(1):end_full)';

% Find missing timestamps
missing_idx = ~ismember(expected_time, MTU_start);
missing_times = expected_time(missing_idx);

% Display results
disp('Missing timestamps:');
disp(missing_times);

fprintf('Total expected timestamps: %d\n', length(expected_time));
fprintf('Actual timestamps: %d\n', length(MTU_start));
fprintf('Missing hours: %d\n', sum(missing_idx));