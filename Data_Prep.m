%% Load and Organize EEG Data from CSV
% This script reads a CSV file containing EEG data, organizes it into a structured
% format, and saves it as a .mat file. It handles empty values and extracts
% relevant columns based on the study requirements.

%% Clear Workspace
clear;
clc;
close all;

%% File Information
filename = 'PT2YY.csv'; % Update with your filename
output_filename = 'PT2YY_Structured_Data.mat';

%% Read CSV File with Original Column Headers
% Read the CSV file into a table, preserving original column headers
data = readtable(filename, 'TextType', 'string', 'VariableNamingRule', 'preserve');

%% Identify Relevant Columns
% EEG Channels
eeg_channels = [
    "EEG.AF3", "EEG.F3", "EEG.F4", ...
    "EEG.FC5", "EEG.FC6", "EEG.AF4"];

% Power Bands (Alpha, High Beta, Theta)
power_bands = struct(...
    "Alpha", [
        "POW.AF3.Alpha", "POW.F3.Alpha", "POW.F4.Alpha", ...
        "POW.FC5.Alpha", "POW.FC6.Alpha", "POW.AF4.Alpha"], ...
    "BetaH", [
        "POW.AF3.BetaH", "POW.F3.BetaH", "POW.F4.BetaH", ...
        "POW.FC5.BetaH", "POW.FC6.BetaH", "POW.AF4.BetaH"], ...
    "Theta", [
        "POW.AF3.Theta", "POW.F3.Theta", "POW.F4.Theta", ...
        "POW.FC5.Theta", "POW.FC6.Theta", "POW.AF4.Theta"]);

% Marker Columns
marker_columns = ["MarkerType", "MarkerValueInt"];

%% Extract and Organize Data
% Initialize structure
eeg_struct = struct(...
    'Timestamp', [], ...
    'EEG_Data', struct(...
        'AF3', [], 'F3', [], 'F4', [], ...
        'FC5', [], 'FC6', [], 'AF4', []), ...
    'Power_Bands', struct(...
        'Alpha', struct(...
            'AF3', [], 'F3', [], 'F4', [], ...
            'FC5', [], 'FC6', [], 'AF4', []), ...
        'BetaH', struct(...
            'AF3', [], 'F3', [], 'F4', [], ...
            'FC5', [], 'FC6', [], 'AF4', []), ...
        'Theta', struct(...
            'AF3', [], 'F3', [], 'F4', [], ...
            'FC5', [], 'FC6', [], 'AF4', [])), ...
    'Markers', struct('Type', [], 'Value', []));

% Extract EEG Data
for i = 1:length(eeg_channels)
    channel_name = eeg_channels(i);
    field_name = channel_name(5:end); % Extract field name (e.g., "AF3")
    eeg_struct.EEG_Data.(field_name) = data{:, channel_name};
end

% Extract Power Bands
for band = fieldnames(power_bands)'
    band_name = band{1};
    band_channels = power_bands.(band_name);
    
    for i = 1:length(band_channels)
        channel_name = band_channels(i);
        field_name = channel_name(5:end); % Extract field name (e.g., "AF3")
        eeg_struct.Power_Bands.(band_name).(field_name) = data{:, channel_name};
    end
end

% Extract Markers
eeg_struct.Markers.Type = data{:, "MarkerType"};
eeg_struct.Markers.Value = data{:, "MarkerValueInt"};

% Extract Timestamps
eeg_struct.Timestamp = data{:, "Timestamp"};

%% Handle Empty Values
% Identify empty values and replace with NaN
eeg_struct = recursively_replace_empty_with_nan(eeg_struct);

%% Save to .mat File
save(output_filename, "eeg_struct");
disp(['Data saved to ', output_filename]);

%% Helper Function to Replace Empty Values with NaN
function struct_data = recursively_replace_empty_with_nan(struct_data)
    fields = fieldnames(struct_data);
    for i = 1:length(fields)
        field = fields{i};
        if isstruct(struct_data.(field))
            struct_data.(field) = recursively_replace_empty_with_nan(struct_data.(field));
        else
            data = struct_data.(field);
            empty_indices = isempty(data) | ismissing(data);
            if any(empty_indices)
                data(empty_indices) = NaN;
                struct_data.(field) = data;
            end
        end
    end
end