%% Clear Workspace
clear;
clc;
close all;

%% 1. Data Extraction and Preparation

% Define participant files
participants = {'PT1NN', 'PT2YY', 'PT3YY', 'PT4YY', 'PT7YN'};
num_participants = numel(participants);

% Initialize data structures
all_data = struct();
all_markers = struct();

% Load data for each participant
for i = 1:num_participants
    participant_id = participants{i};
    all_data(i).id = participant_id;
    all_markers(i).id = participant_id;
    all_data(i).data = readtable([participant_id '.csv'],'VariableNamingRule', 'preserve');
    % Replace Dots with Underscores in Column Headers
    % Get the original column names
    original_vars= all_data(i).data.Properties.VariableNames;
    % Replace dots with underscores
    new_vars = strrep(original_vars, '.', '_');
    % Update the table's variable names
    all_data(i).data.Properties.VariableNames = new_vars;
    all_markers(i).markers = readtable([participant_id '_Interval Marker.csv'],'VariableNamingRule', 'preserve');
    % populate data for other processes
    all_data(i).dataFiltered = all_data(i).data;
    all_data(i).cleaned_data = all_data(i).data;
    all_data(i).ica_components = all_data(i).data;
    all_data(i).epoched_data = all_data(i).data;
end

%% 2. Preprocessing
% Get the column names
all_vars= all_data(i).data.Properties.VariableNames;
% Idendify the channels of interest
channels = {'EEG_AF3', 'EEG_AF4', 'EEG_F3', 'EEG_F4', 'EEG_FC5', ...
    'EEG_FC6', 'EEG_P8', 'EEG_F8'}; % Target  channels here
%channels = {'EEG_AF3','EEG_AF4', 'EEG_F7', 'EEG_F3', 'EEG_F4', 'EEG_FC5',...
%    'EEG_FC6', 'EEG_P8', 'EEG_F8'};

% Find indices of target channels
[~, idx] = ismember(channels, all_vars);

% Load EEGLAB
addpath('eeglab2025.0.0'); % Update with your EEGLAB path
eeglab;

% Sampling frequency
fs = 256;

% Define filter parameters
high_cutoff = 0.8; % High-pass filter at 0.8 Hz
low_cutoff = 40;    % Low-pass filter at 40 Hz

% Apply high-pass and low-pass filters
for i = 1:num_participants
    for ch = 1:length(channels)
        % High-pass filter
        Hpf = eegfilt(all_data(i).data.(idx(ch))', fs, high_cutoff, 0);
        % Low-pass filter
        Lpf = eegfilt(Hpf, fs, 0, low_cutoff);
        all_data(i).dataFiltered.(idx(ch)) = Lpf';
    end
end

%% 3. Artifact Rejection
for i = 1:num_participants
    for ch = 1:length(channels)
        EEG = pop_importdata('dataformat', 'array', 'data', all_data(i).dataFiltered.(idx(ch))', 'srate', fs);
        
        % Perform joint probability artifact rejection (already implemented)
        EEG = pop_jointprob(EEG, 1, 1:EEG.nbchan, 5, 5, 0, 0);
        
        % Apply amplitude thresholding (±100 µV) for artifact rejection
        % Find epochs where the absolute value of any channel exceeds 100 µV
        threshold = 100;  % µV
        artifact_epochs = find(any(abs(EEG.data) > threshold, 1));  % Find epochs exceeding the threshold
        
        % Reject epochs based on amplitude threshold
        EEG = pop_rejepoch(EEG, artifact_epochs, 0);
        
        % Store cleaned data
        all_data(i).cleaned_data.(idx(ch)) = EEG.data';
    end 
end

%% 4. Component Inspection (Multi-Channel ICA)
% Initialize data structures
all_eeg = struct(); % New structure to store EEG data for each participant
for i = 1:num_participants
    participant_id = participants{i};
    %---------------------------------------------------
    % 1) Gather all target EEG channels into one matrix
    %    Rows = channels, Columns = time points
    %---------------------------------------------------
    data_matrix = [];
    for ch = 1:length(channels)
        % Make sure each channel is a row
        data_matrix(ch, :) = all_data(i).cleaned_data.(channels{ch});
    end

    %---------------------------------------------------
    % 2) Create EEGLAB EEG structure from multi-channel data
    %---------------------------------------------------
    EEG = pop_importdata( ...
        'dataformat', 'array', ...
        'data', data_matrix, ...  % multi-channel matrix
        'srate', fs, ...
        'nbchan', size(data_matrix, 1) ...  % number of channels
    );
    
    %---------------------------------------------------
    % Update channel labels to match the 10-20 system
    %---------------------------------------------------
    % Extract channel labels by removing the 'EEG_' prefix
    standard_labels = strrep(channels, 'EEG_', '');
    
    % Assign the updated labels to EEG.chanlocs
    EEG.chanlocs = struct('labels', standard_labels);

    %---------------------------------------------------
    % Load standard 10-20 channel locations
    %---------------------------------------------------
    EEG = pop_chanedit(EEG, 'lookup', 'standard_1020.elc');

    %---------------------------------------------------
    % Ensure EEG structure is correctly configured
    %---------------------------------------------------
    EEG = eeg_checkset(EEG);

    %---------------------------------------------------
    % 3) Run ICA decomposition
    %---------------------------------------------------
    EEG = pop_runica(EEG, 'extended', 1); % You can adjust parameters as needed
    EEG = eeg_checkset(EEG);

    %---------------------------------------------------
    % 4) Plot ICA components & scroll through activations
    %    (Manually identify artifact components)
    %---------------------------------------------------
    pop_topoplot(EEG, 0, 1:EEG.nbchan, ...
        ['Component activations - Participant: ', all_data(i).id]);
    pop_eegplot(EEG, 0, 1, 1);

    %---------------------------------------------------
    % 5) Manually remove artifact components
    %    (Requires user to specify the artifact component indices)
    %---------------------------------------------------
    % For example:
    % artifact_components = [1, 2]; % <--- user-defined
    % EEG = pop_subcomp(EEG, artifact_components, 0);

    %---------------------------------------------------
    % 7) Save EEG structure to a new data structure file
    %---------------------------------------------------
    % Create a new structure to store EEG data for the current participant
    all_eeg(i).id = participant_id;
    all_eeg(i).EEG = EEG;

    % Optionally, save the full EEG structure to a file
    pop_saveset(EEG, 'filename', [all_eeg(i).id '_ICAcleaned.set']);
end

%% 5. Epoching and Baseline Correction
for i = 1:num_participants
    %---------------------------------------------------
    % 1) Import ICA cleaned data as a continuous matrix
    %---------------------------------------------------
    EEG = all_eeg(i).EEG;
    
    %---------------------------------------------------
    % 2) Ensure EEG structure is correctly configured
    %---------------------------------------------------
    EEG = eeg_checkset(EEG);
    
    %---------------------------------------------------
    % 3) Build the event structure from markers, modifying 'W' events
    %---------------------------------------------------
    markers = all_markers(i).markers;
    EEG.event = [];
    
    for j = 1:height(markers)
        event = struct();
        event.latency = round(markers.latency(j) * fs);  % Convert seconds to samples
        event.marker_value = markers.marker_value(j);
        % Check if the event type is 'W' and remove marker_value
        if strncmp(markers.type{j}, 'W',1)
            event.type = 'W';  % Retain only the event type
        else
            event.type = markers.type{j};
        end
        
        EEG.event = [EEG.event, event];
        all_eeg(i).EEG.event = [EEG.event, event];
    end
    
    %---------------------------------------------------
    % 4) Create epochs from -200 ms to 1000 ms relative to events
    %---------------------------------------------------
    EEG = pop_epoch(EEG, {'W', 'T1', 'T2'}, [-0.2 1.0]); % Convert to seconds
    
    %---------------------------------------------------
    % 5) Apply baseline correction using -200 to 0 ms
    %---------------------------------------------------
    EEG = pop_rmbase(EEG, [-0.2 0]); % Convert to seconds
    
    %---------------------------------------------------
    % 6) Store the baseline corrected epochs
    %---------------------------------------------------
    all_eeg(i).epoched_data = EEG.data;
end 
%% 6. Data Preparation: Log-transform
for i = 1:num_participants
    %---------------------------------------------------
    % 1) Access epoched data from all_eeg structure
    %---------------------------------------------------
    epoched_data = all_eeg(i).epoched_data;
    
    %---------------------------------------------------
    % 2) Apply log-transform to each channel
    %---------------------------------------------------
    % Justification,
    %{
    Absolute Values (abs(epoched_data)), the data involves spectral power
    analysis). And we are dealing with event-related spectral perturbations
    (ERSP).   
    %}
    for ch = 1:size(epoched_data, 1)
        % Ensure data is positive to avoid log(0) issues
        epoched_data(ch, :, :) = log10(abs(epoched_data(ch, :, :)) + 1e-6);
    end
    
    %---------------------------------------------------
    % 3) Store the transformed data back into all_eeg
    %---------------------------------------------------
    all_eeg(i).logTransformed_data = epoched_data;
end

%% 7. Data Preparation: Z-score Normalization
for i = 1:num_participants
    %---------------------------------------------------
    % 1) Access log-transformed data from all_eeg structure
    %---------------------------------------------------
    log_data = all_eeg(i).logTransformed_data;
    
    %---------------------------------------------------
    % 2) Compute mean and standard deviation across trials for each channel
    %---------------------------------------------------
    [num_channels, num_timepoints, num_trials] = size(log_data);
    mean_data = mean(log_data, 3);  % Mean across trials
    std_data = std(log_data, 0, 3); % Standard deviation across trials
    
    %---------------------------------------------------
    % 3) Apply Z-score normalization to each trial
    %---------------------------------------------------
    normalized_data = zeros(num_channels, num_timepoints, num_trials);
    for tr = 1:num_trials
        for ch = 1:num_channels
            % Normalize each trial using the mean and std from all trials
            normalized_data(ch, :, tr) = (log_data(ch, :, tr) - mean_data(ch, :)) ./ (std_data(ch, :) + 1e-6);
        end
    end
    
    %---------------------------------------------------
    % 4) Store the normalized data back into all_eeg
    %---------------------------------------------------
    all_eeg(i).normalized_data = normalized_data;
end

%% 8. Data Preparation: Peak Extraction and Analysis
% Define time windows for peak extraction
time_windows = struct('P200', [150, 250], 'P300', [300, 500], 'N400', ...
    [350, 500], 'LPP', [500, 800], 'N200', [250 350]);

% Define EEG frequency bands with their corresponding frequency ranges (in Hz)
frequency_bands = struct(...
    'Theta', [4.8, 7.2], ...
    'Alpha', [9.3, 10.8], ...
    'Beta', [12, 16], ...
    'HighBeta', [18, 25] ...
    );

% Sampling rate (used for calculations, not passed to pop_eegfiltnew)
fs = 256;

% Initialize a structure to store peak data
peak_data = struct();
EEG_data = struct();
EEG_data.nbchan = all_eeg(1).EEG.nbchan; % Number of channels
EEG_data.srate = all_eeg(1).EEG.srate;   % Sampling rate (256 Hz)
% Note: pnts and trials added for compatibility with pop_eegfiltnew

% Iterate over each participant
for i = 1:num_participants
    participant_id = participants{i};
    
    % Load normalized data and set additional fields
    EEG_data.data = all_eeg(i).normalized_data; % Channels x Time x Trials
    EEG_data.pnts = size(all_eeg(i).normalized_data, 2); % Number of time points
    EEG_data.trials = size(all_eeg(i).normalized_data, 3); % Number of trials
    
    [num_channels, num_timepoints, num_trials] = size(EEG_data.data);

    % Initialize structure for each participant
    peak_data(i).id = participant_id;
    peak_data(i).peaks = struct();
    
    % Iterate through each frequency band
    for band_name = fieldnames(frequency_bands)'
        band_name = band_name{1};
        band_range = frequency_bands.(band_name);
        peak_data(i).peaks.(band_name) = struct();
        
        % Calculate transition bandwidth (df) and filter order
        df = max(band_range) * 0.25; % Transition bandwidth as 25% of the upper band edge
        filt_order = round(3.3 * fs / df); % Filter order based on transition bandwidth
        filt_order = filt_order + mod(filt_order, 2); % Ensure filter order is even
        
        % Apply band-pass filter to isolate the frequency band
        EEG_filtered = pop_eegfiltnew(EEG_data, band_range(1), band_range(2), filt_order, 0, 0, 0);
        
        
        % Iterate through each time window (P200, P300, etc.)
        for window_name = fieldnames(time_windows)'
            window_name = window_name{1};
            window = time_windows.(window_name);
            sample_range = round(window * fs / 1000); % Convert time window to sample indices
            
            % Initialize storage for peak data
            peak_data(i).peaks.(band_name).(window_name) = struct();
            
            % Iterate through each channel
            for ch = 1:num_channels
                % Extract data segment for the current channel and time window
                data_segment = EEG_filtered.data(ch, sample_range(1):sample_range(2), :);
                
                % Find positive (max) and negative (min) peaks across trials
                [max_amp, max_idx] = max(data_segment, [], 2);
                [min_amp, min_idx] = min(data_segment, [], 2);
                
                % Convert sample indices to time (in milliseconds)
                time_vector = (sample_range(1):sample_range(2)) / fs * 1000;
                peak_times_max = time_vector(max_idx);
                peak_times_min = time_vector(min_idx);
                
                % Store peak amplitudes and times
                peak_data(i).peaks.(band_name).(window_name).channel(ch).max_amplitude = max_amp;
                peak_data(i).peaks.(band_name).(window_name).channel(ch).max_time = peak_times_max;
                peak_data(i).peaks.(band_name).(window_name).channel(ch).min_amplitude = min_amp;
                peak_data(i).peaks.(band_name).(window_name).channel(ch).min_time = peak_times_min;
            end
        end
    end
end

% Save the extracted peak data for further analysis
save('PeakAnalysisResults.mat', 'peak_data');



%% Frontal Alpha Asymmetry (FAA) Calculation
% Define windows to process
windows = {'P200', 'P300', 'N200', 'LPP'};

% Define trial segments
ANW = 240;    % W: trials 1 to 240
ANT1 = 30;    % T1: trials 241 to 270
ANT2 = 30;    % T2: trials 271 to 300
segments = struct(...
    'W', 1:ANW, ...
    'T1', ANW+1:ANW+ANT1, ...
    'T2', ANW+ANT1+1:ANW+ANT1+ANT2, ...
    'T1T2', ANW+1:ANW+ANT1+ANT2);

% Initialize storage for FAA results
faa_results = struct();

% Step 1: Compute FAA across channels and segment-based metrics
for i = 1:num_participants
    participant_id = participants{i};
    
    % Process each window
    for w = 1:length(windows)
        win = windows{w};
        
        % Extract Alpha Peak Power from F3 and F4 channels (log-transformed per data prep)
        if strcmp(win, 'N200')
            f3_alpha = peak_data(i).peaks.Alpha.(win).channel(3).min_amplitude;
            f4_alpha = peak_data(i).peaks.Alpha.(win).channel(4).min_amplitude;
        else 
            f3_alpha = peak_data(i).peaks.Alpha.(win).channel(3).max_amplitude;
            f4_alpha = peak_data(i).peaks.Alpha.(win).channel(4).max_amplitude;
        end 
        
        num_trials = length(f3_alpha);
        
        % Check if enough trials exist
        if num_trials < ANW + ANT1 + ANT2
            warning(['Participant ', participant_id, ' in window ', win, ...
                     ' has only ', num2str(num_trials), ' trials, less than required ', ...
                     num2str(ANW + ANT1 + ANT2), '. Skipping.']);
            continue;
        end
        
        % Compute reference mean and std for W (control)
        mean_f3_w = mean(f3_alpha(segments.W), 'omitnan');
        std_f3_w = std(f3_alpha(segments.W), 'omitnan') + 1e-6; % Avoid division by zero
        mean_f4_w = mean(f4_alpha(segments.W), 'omitnan');
        std_f4_w = std(f4_alpha(segments.W), 'omitnan') + 1e-6;
        
        % Compute FAA_W (raw difference as baseline)
        FAA_W = mean_f4_w - mean_f3_w;
        if ~isfinite(FAA_W)
            FAA_W = NaN;
        end
        
        % Store participant and window info
        if w == 1
            faa_results(i).id = participant_id;
        end
        faa_results(i).windows.(win) = struct();
        
        % Compute metrics for each segment
        segment_names = fieldnames(segments);
        for s = 1:length(segment_names)
            seg_name = segment_names{s};
            seg_indices = segments.(seg_name);
            
            % Segment data
            f3_seg = f3_alpha(seg_indices);
            f4_seg = f4_alpha(seg_indices);
            
            % Compute trial-level z-scores relative to W
            z_f3_seg = (f3_seg - mean_f3_w) / std_f3_w;
            z_f4_seg = (f4_seg - mean_f4_w) / std_f4_w;
            
            % Compute per-trial FAA (F4 - F3)
            FAA_seg_j = z_f4_seg - z_f3_seg; % Should be 1×N
            FAA_seg_j = squeeze(FAA_seg_j); % Remove singleton dimensions
            if size(FAA_seg_j, 1) > 1 % Ensure row vector
                FAA_seg_j = FAA_seg_j';
            end
            
            % Compute mean FAA and other statistics
            if strcmp(seg_name, 'W')
                FAA_seg = FAA_W; % Use raw difference for W
            else
                FAA_seg = mean(FAA_seg_j, 'omitnan'); % Mean FAA per condition
            end
            sd_faa_seg = std(FAA_seg_j, 'omitnan');
            min_faa_seg = min(FAA_seg_j, [], 'omitnan');
            max_faa_seg = max(FAA_seg_j, [], 'omitnan');
            
            if ~isfinite(FAA_seg)
                FAA_seg = NaN;
                sd_faa_seg = NaN;
                min_faa_seg = NaN;
                max_faa_seg = NaN;
            end
            
            % Store results
            faa_results(i).windows.(win).(seg_name).mean_f3 = mean(f3_seg, 'omitnan');
            faa_results(i).windows.(win).(seg_name).sd_f3 = std(f3_seg, 'omitnan');
            faa_results(i).windows.(win).(seg_name).min_f3 = min(f3_seg, [], 'omitnan');
            faa_results(i).windows.(win).(seg_name).max_f3 = max(f3_seg, [], 'omitnan');
            faa_results(i).windows.(win).(seg_name).mean_f4 = mean(f4_seg, 'omitnan');
            faa_results(i).windows.(win).(seg_name).sd_f4 = std(f4_seg, 'omitnan');
            faa_results(i).windows.(win).(seg_name).min_f4 = min(f4_seg, [], 'omitnan');
            faa_results(i).windows.(win).(seg_name).max_f4 = max(f4_seg, [], 'omitnan');
            faa_results(i).windows.(win).(seg_name).z_f3 = mean(z_f3_seg, 'omitnan');
            faa_results(i).windows.(win).(seg_name).z_f4 = mean(z_f4_seg, 'omitnan');
            faa_results(i).windows.(win).(seg_name).FAA = FAA_seg; % Mean FAA
            faa_results(i).windows.(win).(seg_name).FAA_trials = FAA_seg_j; % Per-trial FAA
            faa_results(i).windows.(win).(seg_name).sd_faa = sd_faa_seg;
            faa_results(i).windows.(win).(seg_name).min_faa = min_faa_seg;
            faa_results(i).windows.(win).(seg_name).max_faa = max_faa_seg;
        end
    end
end

% Step 2: Create Segment-Based FAA Table
% Initialize table with Segment as third column
faa_table = cell(num_participants * length(windows) * length(segment_names) + 1, 15);
faa_table(1, :) = {'Participant', 'Window', 'Segment', ...
                   'Mean F3', 'SD F3', 'Min F3', 'Max F3', ...
                   'Mean F4', 'SD F4', 'Min F4', 'Max F4', ...
                   'Mean FAA(z)', 'SD FAA', 'Min FAA', 'Max FAA'};

row = 2;
for w = 1:length(windows)
    win = windows{w};
    for i = 1:num_participants
        if ~isfield(faa_results, 'id') || isempty(faa_results(i).id) || ...
           ~isfield(faa_results(i).windows, win)
            for s = 1:length(segment_names)
                seg_name = segment_names{s};
                faa_table(row, :) = {participants{i}, win, seg_name, ...
                                     NaN, NaN, NaN, NaN, ... % F3
                                     NaN, NaN, NaN, NaN, ... % F4
                                     NaN, NaN, NaN, NaN};    % FAA
                row = row + 1;
            end
            continue;
        end
        participant_id = faa_results(i).id;
        
        % Fill table for each segment
        for s = 1:length(segment_names)
            seg_name = segment_names{s};
            seg_data = faa_results(i).windows.(win).(seg_name);
            
            faa_table(row, :) = {participant_id, win, seg_name, ...
                                 seg_data.mean_f3, seg_data.sd_f3, seg_data.min_f3, seg_data.max_f3, ...
                                 seg_data.mean_f4, seg_data.sd_f4, seg_data.min_f4, seg_data.max_f4, ...
                                 seg_data.FAA, seg_data.sd_faa, seg_data.min_faa, seg_data.max_faa};
            row = row + 1;
        end
    end
end

% Convert to table
faa_table = cell2table(faa_table(2:end, :), 'VariableNames', faa_table(1, :));

% Display and save the table
disp('Segment-Based FAA Table:');
disp(faa_table);
writetable(faa_table, 'Segment_FAA_Table.csv');


%% Arousal (Beta/Alpha Ratio) Calculation
% Windows to process
windows = {'P300', 'N400', 'LPP'};

% Define trial segments
ANW = 240;    % W: trials 1 to 240
ANT1 = 30;    % T1: trials 241 to 270
ANT2 = 30;    % T2: trials 271 to 300
segments = struct(...
    'W', 1:ANW, ...
    'T1', ANW+1:ANW+ANT1, ...
    'T2', ANW+ANT1+1:ANW+ANT1+ANT2, ...
    'T1T2', ANW+1:ANW+ANT1+ANT2);

% Initialize storage for BAR results
segment_bar_results = struct();

% Step 1: Compute mean beta and alpha across channels and segment-based BAR
for i = 1:num_participants
    participant_id = participants{i};
    
    % Process each window
    for w = 1:length(windows)
        win = windows{w};
        
        % Extract Beta Peak Power from 4 frontal channels
        if strcmp(win, 'N400')
            beta_af3 = peak_data(i).peaks.HighBeta.(win).channel(1).min_amplitude;
            beta_af4 = peak_data(i).peaks.HighBeta.(win).channel(2).min_amplitude;
            beta_f3  = peak_data(i).peaks.HighBeta.(win).channel(3).min_amplitude;
            beta_f4  = peak_data(i).peaks.HighBeta.(win).channel(4).min_amplitude;
            % Extract Alpha Peak Power from 4 frontal channels
            alpha_af3 = peak_data(i).peaks.Alpha.(win).channel(1).min_amplitude;
            alpha_af4 = peak_data(i).peaks.Alpha.(win).channel(2).min_amplitude;
            alpha_f3  = peak_data(i).peaks.Alpha.(win).channel(3).min_amplitude;
            alpha_f4  = peak_data(i).peaks.Alpha.(win).channel(4).min_amplitude;
        else
            beta_af3 = peak_data(i).peaks.HighBeta.(win).channel(1).max_amplitude;
            beta_af4 = peak_data(i).peaks.HighBeta.(win).channel(2).max_amplitude;
            beta_f3  = peak_data(i).peaks.HighBeta.(win).channel(3).max_amplitude;
            beta_f4  = peak_data(i).peaks.HighBeta.(win).channel(4).max_amplitude;
            % Extract Alpha Peak Power from 4 frontal channels
            alpha_af3 = peak_data(i).peaks.Alpha.(win).channel(1).max_amplitude;
            alpha_af4 = peak_data(i).peaks.Alpha.(win).channel(2).max_amplitude;
            alpha_f3  = peak_data(i).peaks.Alpha.(win).channel(3).max_amplitude;
            alpha_f4  = peak_data(i).peaks.Alpha.(win).channel(4).max_amplitude;
        end 
        
        mean_beta = (beta_af3 + beta_af4 + beta_f3 + beta_f4) / 4;
        mean_alpha = (alpha_af3 + alpha_af4 + alpha_f3 + alpha_f4) / 4;
        
        num_trials = length(mean_beta);
        
        % Check if enough trials exist
        if num_trials < ANW + ANT1 + ANT2
            warning(['Participant ', participant_id, ' in window ', win, ...
                     ' has only ', num2str(num_trials), ' trials, less than required ', ...
                     num2str(ANW + ANT1 + ANT2), '. Skipping.']);
            continue;
        end
        
        % Compute reference mean and std for W (for z-scores)
        mean_beta_w = mean(mean_beta(segments.W), 'omitnan');
        std_beta_w = std(mean_beta(segments.W), 'omitnan') + 1e-6; % Avoid division by zero
        mean_alpha_w = mean(mean_alpha(segments.W), 'omitnan');
        std_alpha_w = std(mean_alpha(segments.W), 'omitnan') + 1e-6;
        
        % Compute BAR_W (dynamic baseline)
        BAR_W = mean_beta_w / mean_alpha_w;
        if ~isfinite(BAR_W)
            BAR_W = NaN; % Handle invalid baseline
        end
        
        % Store participant and window info
        if w == 1
            segment_bar_results(i).id = participant_id;
        end
        segment_bar_results(i).windows.(win) = struct();
        segment_bar_results(i).windows.(win).acrosschannel.mean_beta = mean_beta;
        segment_bar_results(i).windows.(win).acrosschannel.mean_alpha = mean_alpha;
        % Compute metrics for each segment
        segment_names = fieldnames(segments);
        for s = 1:length(segment_names)
            seg_name = segment_names{s};
            seg_indices = segments.(seg_name);
            
            % Segment data
            beta_seg = mean_beta(seg_indices);
            alpha_seg = mean_alpha(seg_indices);
            
            % Compute statistics
            mean_beta_seg = mean(beta_seg, 'omitnan');
            sd_beta_seg = std(beta_seg, 'omitnan');
            min_beta_seg = min(beta_seg, [], 'omitnan');
            max_beta_seg = max(beta_seg, [], 'omitnan');
            
            mean_alpha_seg = mean(alpha_seg, 'omitnan');
            sd_alpha_seg = std(alpha_seg, 'omitnan');
            min_alpha_seg = min(alpha_seg, [], 'omitnan');
            max_alpha_seg = max(alpha_seg, [], 'omitnan');
            
            % Compute trial-level z-scores relative to W
            z_beta_seg = (beta_seg - mean_beta_w) / std_beta_w;
            z_alpha_seg = (alpha_seg - mean_alpha_w) / std_alpha_w;
            
            % Compute BAR
            if strcmp(win, 'N400')
                if strcmp(seg_name, 'W')
                    BAR_seg = abs(BAR_W); % Use raw mean ratio for W
                    BAR_seg_j = abs(z_beta_seg ./ z_alpha_seg);
                    mean_z_beta_seg = mean(z_beta_seg, 'omitnan'); 
                    mean_z_alpha_seg = mean(z_alpha_seg, 'omitnan');
                else
                    BAR_seg_j = abs(z_beta_seg ./ z_alpha_seg);
                    mean_z_beta_seg = mean(z_beta_seg, 'omitnan'); 
                    mean_z_alpha_seg = mean(z_alpha_seg, 'omitnan');
                    BAR_seg = mean_z_beta_seg/mean_z_alpha_seg;
                    %BAR_seg = mean(BAR_seg_j, 'omitnan');
                    
                    if ~isfinite(BAR_seg)
                        BAR_seg = NaN; % Handle NaN/infinite BAR
                    end
                end
            else 
                if strcmp(seg_name, 'W')
                    BAR_seg = BAR_W; % Use raw mean ratio for W
                    BAR_seg_j = abs(z_beta_seg ./ z_alpha_seg);
                    mean_z_beta_seg = mean(z_beta_seg, 'omitnan'); 
                    mean_z_alpha_seg = mean(z_alpha_seg, 'omitnan');
                    
                else
                    BAR_seg_j = z_beta_seg ./ z_alpha_seg;
                    mean_z_beta_seg = mean(z_beta_seg, 'omitnan'); 
                    mean_z_alpha_seg = mean(z_alpha_seg, 'omitnan');
                    BAR_seg = mean_z_beta_seg/mean_z_alpha_seg;
                    %BAR_seg = mean(BAR_seg_j, 'omitnan');
                    if ~isfinite(BAR_seg)
                        BAR_seg = NaN; % Handle NaN/infinite BAR
                    end
                end
            end 
            
            % Store mean z-scores for T1 and T2 to compute T1T2 BAR_change
            if strcmp(seg_name, 'T1')
                z_beta_t1 = mean(z_beta_seg, 'omitnan');
                z_alpha_t1 = mean(z_alpha_seg, 'omitnan');
                BAR_seg_j = abs(z_beta_seg ./ z_alpha_seg);
                mean_z_beta_seg = mean(z_beta_seg, 'omitnan'); 
                mean_z_alpha_seg = mean(z_alpha_seg, 'omitnan');
                BAR_seg = mean_z_beta_seg/mean_z_alpha_seg;
            elseif strcmp(seg_name, 'T2')
                z_beta_t2 = mean(z_beta_seg, 'omitnan');
                z_alpha_t2 = mean(z_alpha_seg, 'omitnan');
                BAR_seg_j = abs(z_beta_seg ./ z_alpha_seg);
                mean_z_beta_seg = mean(z_beta_seg, 'omitnan'); 
                mean_z_alpha_seg = mean(z_alpha_seg, 'omitnan');
                BAR_seg = mean_z_beta_seg/mean_z_alpha_seg;
            end
            
            % Compute BAR_change
            if strcmp(seg_name, 'W')
                BAR_change = 0; % Baseline, no change
            elseif strcmp(seg_name, 'T1T2')
                % T1T2: Ratio of z-score differences
                BAR_change = (abs((abs(z_beta_t2 - z_alpha_t2)/abs(z_beta_t1 - z_alpha_t1 + 1e-6))-BAR_W)/BAR_W)*100;
                if ~isfinite(BAR_change)
                    BAR_change = NaN; % Handle invalid ratio
                end
            else
                % T1 and T2: Percent change relative to BAR_W
                BAR_change = (abs(BAR_seg - BAR_W) / BAR_W) * 100;
                if ~isfinite(BAR_change)
                    BAR_change = NaN; % Handle invalid percent change
                end
            end
            
            % Store results
            segment_bar_results(i).windows.(win).(seg_name).mean_beta = mean_beta_seg;
            segment_bar_results(i).windows.(win).(seg_name).sd_beta = sd_beta_seg;
            segment_bar_results(i).windows.(win).(seg_name).min_beta = min_beta_seg;
            segment_bar_results(i).windows.(win).(seg_name).max_beta = max_beta_seg;
            segment_bar_results(i).windows.(win).(seg_name).mean_alpha = mean_alpha_seg;
            segment_bar_results(i).windows.(win).(seg_name).sd_alpha = sd_alpha_seg;
            segment_bar_results(i).windows.(win).(seg_name).min_alpha = min_alpha_seg;
            segment_bar_results(i).windows.(win).(seg_name).max_alpha = max_alpha_seg;
            segment_bar_results(i).windows.(win).(seg_name).z_beta = mean_z_beta_seg;
            segment_bar_results(i).windows.(win).(seg_name).z_alpha =  mean_z_alpha_seg;
            segment_bar_results(i).windows.(win).(seg_name).BAR = BAR_seg;
            segment_bar_results(i).windows.(win).(seg_name).BAR_trials = BAR_seg_j;
            segment_bar_results(i).windows.(win).(seg_name).BAR_change = BAR_change;
        end
    end
end

% Step 2: Create Segment-Based Table for All Windows with Changes
% Initialize table
segment_table = cell(num_participants * length(windows) * length(segment_names) + 1, 15);
segment_table(1, :) = {'Participant', 'Window', 'Segment', ...
                       'Mean Beta', 'SD Beta', 'Min Beta', 'Max Beta', ...
                       'Mean Alpha', 'SD Alpha', 'Min Alpha', 'Max Alpha', ...
                       'MeanBeta(z)','MeanAlpha(z)', 'MeanBAR(z)', 'BAR Change'};

row = 2;
for w = 1:length(windows)
    win = windows{w};
    for i = 1:num_participants
        if ~isfield(segment_bar_results, 'id') || isempty(segment_bar_results(i).id) || ...
           ~isfield(segment_bar_results(i).windows, win)
            continue; % Skip if participant/window was excluded
        end
        participant_id = segment_bar_results(i).id;
        
        % Fill table for all segments
        for s = 1:length(segment_names)
            seg_name = segment_names{s};
            seg_data = segment_bar_results(i).windows.(win).(seg_name);
            
            segment_table(row, :) = {participant_id, win, seg_name, ...
                                     seg_data.mean_beta, seg_data.sd_beta, ...
                                     seg_data.min_beta, seg_data.max_beta, ...
                                     seg_data.mean_alpha, seg_data.sd_alpha, ...
                                     seg_data.min_alpha, seg_data.max_alpha, ...
                                     seg_data.z_beta, seg_data.z_alpha,seg_data.BAR, ...
                                     seg_data.BAR_change};
            row = row + 1;
        end
    end
end

% Convert to table, excluding empty rows if any
segment_table = cell2table(segment_table(2:row-1, :), 'VariableNames', segment_table(1, :));

% Display and save the table
disp('Segment-Based Beta/Alpha Ratio (BAR) Table with Changes for All Windows:');
disp(segment_table);
writetable(segment_table, 'Segment_BAR_Table_with_Changes.csv');

%% High Beta Power (Stress) Calculation
% Define windows to process
windows = {'N200', 'P200', 'P300', 'N400', 'LPP'};

% Define trial segments
ANW = 240;    % W: trials 1 to 240
ANT1 = 30;    % T1: trials 241 to 270
ANT2 = 30;    % T2: trials 271 to 300
segments = struct(...
    'W', 1:ANW, ...
    'T1', ANW+1:ANW+ANT1, ...
    'T2', ANW+ANT1+1:ANW+ANT1+ANT2, ...
    'T1T2', ANW+1:ANW+ANT1+ANT2);

% Initialize storage for High Beta Power results
hb_results = struct();

% Step 1: Compute High Beta Power across channels and segment-based metrics
for i = 1:num_participants
    participant_id = participants{i};
    
    % Process each window
    for w = 1:length(windows)
        win = windows{w};
        
        % Extract High Beta Peak Power from 4 frontal channels (AF3, AF4, F3, F4)
        if strcmp(win, 'N200') || strcmp(win, 'N400')
            hb_af3 = peak_data(i).peaks.HighBeta.(win).channel(1).min_amplitude;
            hb_af4 = peak_data(i).peaks.HighBeta.(win).channel(2).min_amplitude;
            hb_f3 = peak_data(i).peaks.HighBeta.(win).channel(3).min_amplitude;
            hb_f4 = peak_data(i).peaks.HighBeta.(win).channel(4).min_amplitude;
        else
            hb_af3 = peak_data(i).peaks.HighBeta.(win).channel(1).max_amplitude;
            hb_af4 = peak_data(i).peaks.HighBeta.(win).channel(2).max_amplitude;
            hb_f3 = peak_data(i).peaks.HighBeta.(win).channel(3).max_amplitude;
            hb_f4 = peak_data(i).peaks.HighBeta.(win).channel(4).max_amplitude;
        end 
        
        % Compute mean High Beta Power across 4 channels per trial
        mean_hb = (hb_af3 + hb_af4 + hb_f3 + hb_f4) / 4;
        
        num_trials = length(mean_hb);
        
        % Check if enough trials exist
        if num_trials < ANW + ANT1 + ANT2
            warning(['Participant ', participant_id, ' in window ', win, ...
                     ' has only ', num2str(num_trials), ' trials, less than required ', ...
                     num2str(ANW + ANT1 + ANT2), '. Skipping.']);
            continue;
        end
        
        % Compute reference mean and std for W (control)
        mean_hb_w = mean(mean_hb(segments.W), 'omitnan');
        std_hb_w = std(mean_hb(segments.W), 'omitnan') + 1e-6; % Avoid division by zero
        
        % Store participant and window info
        if w == 1
            hb_results(i).id = participant_id;
        end
        hb_results(i).windows.(win) = struct();
        hb_results(i).windows.(win).acrosschannel.mean_hb = mean_hb;
        % Compute metrics for each segment
        segment_names = fieldnames(segments);
        for s = 1:length(segment_names)
            seg_name = segment_names{s};
            seg_indices = segments.(seg_name);
            
            % Segment data
            hb_seg = mean_hb(seg_indices);
            f3_seg = hb_f3(seg_indices);
            f4_seg = hb_f4(seg_indices);
            af3_seg = hb_af3(seg_indices);
            af4_seg = hb_af4(seg_indices);
            % Compute statistics for overall HB
            mean_hb_seg = mean(hb_seg, 'omitnan');
            sd_hb_seg = std(hb_seg, 'omitnan');
            min_hb_seg = min(hb_seg, [], 'omitnan');
            max_hb_seg = max(hb_seg, [], 'omitnan');
            
            % Compute statistics for F3
            mean_f3_seg = mean(f3_seg, 'omitnan');
            sd_f3_seg = std(f3_seg, 'omitnan');
            min_f3_seg = min(f3_seg, [], 'omitnan');
            max_f3_seg = max(f3_seg, [], 'omitnan');
            
            % Compute statistics for F4
            mean_f4_seg = mean(f4_seg, 'omitnan');
            sd_f4_seg = std(f4_seg, 'omitnan');
            min_f4_seg = min(f4_seg, [], 'omitnan');
            max_f4_seg = max(f4_seg, [], 'omitnan');
            
            % Compute statistics for AF3
            mean_af3_seg = mean(af3_seg, 'omitnan');
            sd_af3_seg = std(af3_seg, 'omitnan');
            min_af3_seg = min(af3_seg, [], 'omitnan');
            max_af3_seg = max(af3_seg, [], 'omitnan');
            
            % Compute statistics for AF4
            mean_af4_seg = mean(af4_seg, 'omitnan');
            sd_af4_seg = std(af4_seg, 'omitnan');
            min_af4_seg = min(af4_seg, [], 'omitnan');
            max_af4_seg = max(af4_seg, [], 'omitnan');
            
            % Compute trial-level z-scores relative to W for overall HB
            z_hb_seg = (hb_seg - mean_hb_w) / std_hb_w;
            z_hb_seg = squeeze(z_hb_seg); % Remove singleton dimensions
            if size(z_hb_seg, 1) > 1 % Ensure row vector
                z_hb_seg = z_hb_seg';
            end
            
            % Compute mean z-score for the segment
            mean_z_hb_seg = mean(z_hb_seg, 'omitnan');
            if ~isfinite(mean_z_hb_seg)
                mean_z_hb_seg = NaN;
                sd_hb_seg = NaN;
                min_hb_seg = NaN;
                max_hb_seg = NaN;
                mean_f3_seg = NaN; sd_f3_seg = NaN; min_f3_seg = NaN; max_f3_seg = NaN;
                mean_f4_seg = NaN; sd_f4_seg = NaN; min_f4_seg = NaN; max_f4_seg = NaN;
                mean_af3_seg = NaN; sd_af3_seg = NaN; min_af3_seg = NaN; max_af3_seg = NaN;
                mean_af4_seg = NaN; sd_af4_seg = NaN; min_af4_seg = NaN; max_af4_seg = NaN;
            end
            
            % Store mean z-scores for T1 and T2 to compute T1T2 HB_change
            if strcmp(seg_name, 'T1')
                z_hb_t1 = mean_z_hb_seg;
            elseif strcmp(seg_name, 'T2')
                z_hb_t2 = mean_z_hb_seg;
            end
            
            % Compute HB_change
            if strcmp(seg_name, 'W')
                HB_change = 0; % Baseline, no change
            elseif strcmp(seg_name, 'T1T2')
                % T1T2: Ratio of z-score differences
                HB_change = z_hb_t2 / (z_hb_t1 + 1e-6);
                if ~isfinite(HB_change)
                    HB_change = NaN;
                end
            else
                % T1 and T2: Percent change relative to mean_hb_w
                HB_change = (abs(mean_hb_seg - mean_hb_w) / mean_hb_w) * 100;
                if ~isfinite(HB_change)
                    HB_change = NaN;
                end
            end
            
            % Store results
            hb_results(i).windows.(win).(seg_name).mean_hb = mean_hb_seg;
            hb_results(i).windows.(win).(seg_name).sd_hb = sd_hb_seg;
            hb_results(i).windows.(win).(seg_name).min_hb = min_hb_seg;
            hb_results(i).windows.(win).(seg_name).max_hb = max_hb_seg;
            hb_results(i).windows.(win).(seg_name).mean_f3 = mean_f3_seg;
            hb_results(i).windows.(win).(seg_name).sd_f3 = sd_f3_seg;
            hb_results(i).windows.(win).(seg_name).min_f3 = min_f3_seg;
            hb_results(i).windows.(win).(seg_name).max_f3 = max_f3_seg;
            hb_results(i).windows.(win).(seg_name).mean_f4 = mean_f4_seg;
            hb_results(i).windows.(win).(seg_name).sd_f4 = sd_f4_seg;
            hb_results(i).windows.(win).(seg_name).min_f4 = min_f4_seg;
            hb_results(i).windows.(win).(seg_name).max_f4 = max_f4_seg;
            hb_results(i).windows.(win).(seg_name).mean_af3 = mean_af3_seg;
            hb_results(i).windows.(win).(seg_name).sd_af3 = sd_af3_seg;
            hb_results(i).windows.(win).(seg_name).min_af3 = min_af3_seg;
            hb_results(i).windows.(win).(seg_name).max_af3 = max_af3_seg;
            hb_results(i).windows.(win).(seg_name).mean_af4 = mean_af4_seg;
            hb_results(i).windows.(win).(seg_name).sd_af4 = sd_af4_seg;
            hb_results(i).windows.(win).(seg_name).min_af4 = min_af4_seg;
            hb_results(i).windows.(win).(seg_name).max_af4 = max_af4_seg;
            hb_results(i).windows.(win).(seg_name).z_hb = mean_z_hb_seg;
            hb_results(i).windows.(win).(seg_name).hb_trials = z_hb_seg; % Per-trial z-scores
            hb_results(i).windows.(win).(seg_name).HB_change = HB_change;
        end
    end
end

% Step 2: Create Segment-Based High Beta Power Table
% Initialize table with Segment as third column
hb_table = cell(num_participants * length(windows) * length(segment_names) + 1, 23);
hb_table(1, :) = {'Participant', 'Window', 'Segment', ...
                  'MeanHB(z)', 'SD HB', 'Min HB', 'Max HB', ...
                  'Mean F3', 'SD F3', 'Min F3', 'Max F3', ...
                  'Mean F4', 'SD F4', 'Min F4', 'Max F4',...
                  'Mean AF3', 'SD AF3', 'Min AF3', 'Max AF3', ...
                  'Mean AF4', 'SD AF4', 'Min AF4', 'Max AF4'};

row = 2;
for w = 1:length(windows)
    win = windows{w};
    for i = 1:num_participants
        if ~isfield(hb_results, 'id') || isempty(hb_results(i).id) || ...
           ~isfield(hb_results(i).windows, win)
            for s = 1:length(segment_names)
                seg_name = segment_names{s};
                hb_table(row, :) = {participants{i}, win, seg_name, ...
                                    NaN, NaN, NaN, NaN, ... % HB
                                    NaN, NaN, NaN, NaN, ... % F3
                                    NaN, NaN, NaN, NaN, ... % F4
                                    NaN, NaN, NaN, NaN, ... % AF3
                                    NaN, NaN, NaN, NaN};    % AF4
                row = row + 1;
            end
            continue;
        end
        participant_id = hb_results(i).id;
        
        % Fill table for each segment
        for s = 1:length(segment_names)
            seg_name = segment_names{s};
            seg_data = hb_results(i).windows.(win).(seg_name);
            
            hb_table(row, :) = {participant_id, win, seg_name, ...
                                seg_data.mean_hb, seg_data.sd_hb, seg_data.min_hb, seg_data.max_hb, ...
                                seg_data.mean_f3, seg_data.sd_f3, seg_data.min_f3, seg_data.max_f3, ...
                                seg_data.mean_f4, seg_data.sd_f4, seg_data.min_f4, seg_data.max_f4, ...
                                seg_data.mean_af3, seg_data.sd_af3, seg_data.min_af3, seg_data.max_af3,...
                                seg_data.mean_af4, seg_data.sd_af4, seg_data.min_af4, seg_data.max_af4};
            row = row + 1;
        end
    end
end

% Convert to table
hb_table = cell2table(hb_table(2:end, :), 'VariableNames', hb_table(1, :));

% Display and save the table
disp('Segment-Based High Beta Power (Stress) Table:');
disp(hb_table);
writetable(hb_table, 'Segment_High_Beta_Power_Table.csv');

%% Theta/Beta Ratio (TBR) Calculation for Target Trials (P300, N400, LPP)
% Define windows to process for TBR
tbr_windows = {'P300', 'P200'};

% Define trial segments
ANW = 240;    % W: trials 1 to 240
ANT1 = 30;    % T1: trials 241 to 270
ANT2 = 30;    % T2: trials 271 to 300
segments = struct(...
    'W', 1:ANW, ...
    'T1', ANW+1:ANW+ANT1, ...
    'T2', ANW+ANT1+1:ANW+ANT1+ANT2, ...
    'T1T2', ANW+1:ANW+ANT1+ANT2);

% Initialize storage for TBR results
tbr_results = struct();

% Step 1: Compute mean Theta and Beta across channels and segment-based TBR
for i = 1:num_participants
    participant_id = participants{i};
    
    % Process each window
    for w = 1:length(tbr_windows)
        win = tbr_windows{w};
        if strcmp(win, 'N400')
            % Extract Theta Peak Power from 4 frontal channels
            theta_af3 = peak_data(i).peaks.Theta.(win).channel(1).min_amplitude;
            theta_af4 = peak_data(i).peaks.Theta.(win).channel(2).min_amplitude;
            theta_f3  = peak_data(i).peaks.Theta.(win).channel(3).min_amplitude;
            theta_f4  = peak_data(i).peaks.Theta.(win).channel(4).min_amplitude;     
            % Extract Beta Peak Power from 4 frontal channels
            beta_af3 = peak_data(i).peaks.HighBeta.(win).channel(1).min_amplitude;
            beta_af4 = peak_data(i).peaks.HighBeta.(win).channel(2).min_amplitude;
            beta_f3  = peak_data(i).peaks.HighBeta.(win).channel(3).min_amplitude;
            beta_f4  = peak_data(i).peaks.HighBeta.(win).channel(4).min_amplitude;
        else
            % Extract Theta Peak Power from 4 frontal channels
            theta_af3 = peak_data(i).peaks.Theta.(win).channel(1).max_amplitude;
            theta_af4 = peak_data(i).peaks.Theta.(win).channel(2).max_amplitude;
            theta_f3  = peak_data(i).peaks.Theta.(win).channel(3).max_amplitude;
            theta_f4  = peak_data(i).peaks.Theta.(win).channel(4).max_amplitude;     
            % Extract Beta Peak Power from 4 frontal channels
            beta_af3 = peak_data(i).peaks.HighBeta.(win).channel(1).max_amplitude;
            beta_af4 = peak_data(i).peaks.HighBeta.(win).channel(2).max_amplitude;
            beta_f3  = peak_data(i).peaks.HighBeta.(win).channel(3).max_amplitude;
            beta_f4  = peak_data(i).peaks.HighBeta.(win).channel(4).max_amplitude;
            
        end
        
        mean_theta = (theta_af3 + theta_af4 + theta_f3 + theta_f4) / 4;
        mean_beta = (beta_af3 + beta_af4 + beta_f3 + beta_f4) / 4;
        
        num_trials = length(mean_theta);
        
        % Check if enough trials exist
        if num_trials < ANW + ANT1 + ANT2
            warning(['Participant ', participant_id, ' in window ', win, ...
                     ' has only ', num2str(num_trials), ' trials, less than required ', ...
                     num2str(ANW + ANT1 + ANT2), '. Skipping.']);
            continue;
        end
        
        % Compute reference mean and std for W (for z-scores)
        mean_theta_w = mean(mean_theta(segments.W), 'omitnan');
        std_theta_w = std(mean_theta(segments.W), 'omitnan') + 1e-6; % Avoid division by zero
        mean_beta_w = mean(mean_beta(segments.W), 'omitnan');
        std_beta_w = std(mean_beta(segments.W), 'omitnan') + 1e-6;
        
        % Compute TBR_W (dynamic baseline)
        TBR_W = mean_theta_w / mean_beta_w;
        if ~isfinite(TBR_W)
            TBR_W = NaN; % Handle invalid baseline
        end
        
        % Store participant and window info
        if w == 1
            tbr_results(i).id = participant_id;
        end
        tbr_results(i).windows.(win) = struct();
        
        % Compute metrics for each segment
        segment_names = fieldnames(segments);
        for s = 1:length(segment_names)
            seg_name = segment_names{s};
            seg_indices = segments.(seg_name);
            
            % Segment data
            theta_seg = mean_theta(seg_indices);
            beta_seg = mean_beta(seg_indices);
            
            % Compute statistics
            mean_theta_seg = mean(theta_seg, 'omitnan');
            sd_theta_seg = std(theta_seg, 'omitnan');
            min_theta_seg = min(theta_seg, [], 'omitnan');
            max_theta_seg = max(theta_seg, [], 'omitnan');
            
            mean_beta_seg = mean(beta_seg, 'omitnan');
            sd_beta_seg = std(beta_seg, 'omitnan');
            min_beta_seg = min(beta_seg, [], 'omitnan');
            max_beta_seg = max(beta_seg, [], 'omitnan');
            
            % Compute trial-level z-scores relative to W
            z_theta_seg = (theta_seg - mean_theta_w) / std_theta_w;
            z_beta_seg = (beta_seg - mean_beta_w) / std_beta_w;
            
            % Compute TBR
            if strcmp(win, 'N400')
                if strcmp(seg_name, 'W')
                    TBR_seg = TBR_W; % Use raw mean ratio for W
                    mean_z_theta_seg = mean(z_theta_seg, 'omitnan');
                    mean_z_beta_seg = mean(z_beta_seg, 'omitnan');
                    TBR_seg_j = z_theta_seg ./ abs(z_beta_seg);
                else
                    TBR_seg_j = z_theta_seg ./ abs(z_beta_seg);
                    %TBR_seg = mean(TBR_seg_j, 'omitnan');
                    mean_z_theta_seg = mean(z_theta_seg, 'omitnan');
                    mean_z_beta_seg = mean(z_beta_seg, 'omitnan');
                    TBR_seg = mean_z_theta_seg/ mean_z_beta_seg ;
                    if ~isfinite(TBR_seg)
                        TBR_seg = NaN; % Handle NaN/infinite TBR
                    end
                end
            else 
                if strcmp(seg_name, 'W')
                    TBR_seg = abs(TBR_W); % Use raw mean ratio for W
                    TBR_seg_j = z_theta_seg ./ abs(z_beta_seg);
                    mean_z_theta_seg = mean(z_theta_seg, 'omitnan');
                    mean_z_beta_seg = mean(z_beta_seg, 'omitnan');
                    
                else
                    TBR_seg_j = (z_theta_seg ./ z_beta_seg);
                    %TBR_seg = mean(TBR_seg_j, 'omitnan');
                    mean_z_theta_seg = mean(z_theta_seg, 'omitnan');
                    mean_z_beta_seg = mean(z_beta_seg, 'omitnan');
                    TBR_seg = mean_z_theta_seg/ mean_z_beta_seg ;
                    if ~isfinite(TBR_seg)
                        TBR_seg = NaN; % Handle NaN/infinite TBR
                    end
                end
            end 
            
            % Store mean z-scores for T1 and T2 to compute T1T2 TBR_change
            if strcmp(seg_name, 'T1')
                z_theta_t1 = mean(z_theta_seg, 'omitnan');
                z_beta_t1 = mean(z_beta_seg, 'omitnan');
                TBR_seg_j = z_theta_seg ./ abs(z_beta_seg);
                mean_z_theta_seg = mean(z_theta_seg, 'omitnan');
                mean_z_beta_seg = mean(z_beta_seg, 'omitnan');
                TBR_seg = mean_z_theta_seg/ mean_z_beta_seg ;
            elseif strcmp(seg_name, 'T2')
                z_theta_t2 = mean(z_theta_seg, 'omitnan');
                z_beta_t2 = mean(z_beta_seg, 'omitnan');
                TBR_seg_j = z_theta_seg ./ abs(z_beta_seg);
                mean_z_theta_seg = mean(z_theta_seg, 'omitnan');
                mean_z_beta_seg = mean(z_beta_seg, 'omitnan');
                TBR_seg = mean_z_theta_seg/ mean_z_beta_seg ;
            end
            
            % Compute TBR_change
            if strcmp(seg_name, 'W')
                TBR_change = 0; % Baseline, no change
            elseif strcmp(seg_name, 'T1T2')
                % T1T2: Ratio of z-score differences
                TBR_change = abs((z_theta_t2 - z_beta_t2) / (z_theta_t1 - z_beta_t1 + 1e-6));
                if ~isfinite(TBR_change)
                    TBR_change = NaN; % Handle invalid ratio
                end
            else
                % T1 and T2: Percent change relative to TBR_W
                TBR_change = (abs(TBR_seg - TBR_W) / TBR_W) * 100;
                if ~isfinite(TBR_change)
                    TBR_change = NaN; % Handle invalid percent change
                end
            end
            
            % Store results
            tbr_results(i).windows.(win).(seg_name).mean_theta = mean_theta_seg;
            tbr_results(i).windows.(win).(seg_name).sd_theta = sd_theta_seg;
            tbr_results(i).windows.(win).(seg_name).min_theta = min_theta_seg;
            tbr_results(i).windows.(win).(seg_name).max_theta = max_theta_seg;
            tbr_results(i).windows.(win).(seg_name).mean_beta = mean_beta_seg;
            tbr_results(i).windows.(win).(seg_name).sd_beta = sd_beta_seg;
            tbr_results(i).windows.(win).(seg_name).min_beta = min_beta_seg;
            tbr_results(i).windows.(win).(seg_name).max_beta = max_beta_seg;
            tbr_results(i).windows.(win).(seg_name).z_theta = mean_z_theta_seg;
            tbr_results(i).windows.(win).(seg_name).z_beta = mean_z_beta_seg;
            tbr_results(i).windows.(win).(seg_name).TBR = TBR_seg;
            tbr_results(i).windows.(win).(seg_name).TBR_change = TBR_change;
            tbr_results(i).windows.(win).(seg_name).TBR_trials = TBR_seg_j;
        end
    end
end

% Step 2: Create Segment-Based Table for All Windows with Changes
% Initialize table
tbr_table = cell(num_participants * length(tbr_windows) * length(segment_names) + 1, 14);
tbr_table(1, :) = {'Participant', 'Window', 'Segment', ...
                   'Mean Theta', 'SD Theta', 'Min Theta', 'Max Theta', ...
                   'Mean Beta', 'SD Beta', 'Min Beta', 'Max Beta', ...
                   'meanTheta(z)','meanBeta(z)' 'TBR'};

row = 2;
for w = 1:length(tbr_windows)
    win = tbr_windows{w};
    for i = 1:num_participants
        if ~isfield(tbr_results, 'id') || isempty(tbr_results(i).id) || ...
           ~isfield(tbr_results(i).windows, win)
            continue; % Skip if participant/window was excluded
        end
        participant_id = tbr_results(i).id;
        
        % Fill table for all segments
        for s = 1:length(segment_names)
            seg_name = segment_names{s};
            seg_data = tbr_results(i).windows.(win).(seg_name);
            
            tbr_table(row, :) = {participant_id, win, seg_name, ...
                                 seg_data.mean_theta, seg_data.sd_theta, ...
                                 seg_data.min_theta, seg_data.max_theta, ...
                                 seg_data.mean_beta, seg_data.sd_beta, ...
                                 seg_data.min_beta, seg_data.max_beta, ...
                                 seg_data.z_theta,seg_data.z_beta seg_data.TBR}; %...
                                 %seg_data.TBR_change};
            row = row + 1;
        end
    end
end

% Convert to table, excluding empty rows if any
tbr_table = cell2table(tbr_table(2:row-1, :), 'VariableNames', tbr_table(1, :));

% Display and save the table
disp('Segment-Based Theta/Beta Ratio (TBR) Table with Changes for All Windows:');
disp(tbr_table);
writetable(tbr_table, 'Segment_TBR_Table_with_Changes.csv');

%% Theta/Beta Ratio (TBR) Calculation for Target Trials (N200)
% Define windows to process for TBR
tbr_windows = {'N200'};

% Define trial segments
ANW = 240;    % W: trials 1 to 240
ANT1 = 30;    % T1: trials 241 to 270
ANT2 = 30;    % T2: trials 271 to 300
segments = struct(...
    'W', 1:ANW, ...
    'T1', ANW+1:ANW+ANT1, ...
    'T2', ANW+ANT1+1:ANW+ANT1+ANT2, ...
    'T1T2', ANW+1:ANW+ANT1+ANT2);

% Initialize storage for TBR results
tbr_results = struct();

% Step 1: Compute mean Theta and Beta across channels and segment-based TBR
for i = 1:num_participants
    participant_id = participants{i};
    
    % Process each window
    for w = 1:length(tbr_windows)
        win = tbr_windows{w};
        if strcmp(win, 'N200')
            % Extract Theta Peak Power from 4 frontal channels
            theta_af3 = peak_data(i).peaks.Theta.(win).channel(1).min_amplitude;
            theta_af4 = peak_data(i).peaks.Theta.(win).channel(2).min_amplitude;
            theta_f3  = peak_data(i).peaks.Theta.(win).channel(3).min_amplitude;
            theta_f4  = peak_data(i).peaks.Theta.(win).channel(4).min_amplitude;     
            % Extract Beta Peak Power from 4 frontal channels
            beta_af3 = peak_data(i).peaks.HighBeta.(win).channel(1).min_amplitude;
            beta_af4 = peak_data(i).peaks.HighBeta.(win).channel(2).min_amplitude;
            beta_f3  = peak_data(i).peaks.HighBeta.(win).channel(3).min_amplitude;
            beta_f4  = peak_data(i).peaks.HighBeta.(win).channel(4).min_amplitude;
        else
            % Extract Theta Peak Power from 4 frontal channels
            theta_af3 = peak_data(i).peaks.Theta.(win).channel(1).max_amplitude;
            theta_af4 = peak_data(i).peaks.Theta.(win).channel(2).max_amplitude;
            theta_f3  = peak_data(i).peaks.Theta.(win).channel(3).max_amplitude;
            theta_f4  = peak_data(i).peaks.Theta.(win).channel(4).max_amplitude;     
            % Extract Beta Peak Power from 4 frontal channels
            beta_af3 = peak_data(i).peaks.HighBeta.(win).channel(1).max_amplitude;
            beta_af4 = peak_data(i).peaks.HighBeta.(win).channel(2).max_amplitude;
            beta_f3  = peak_data(i).peaks.HighBeta.(win).channel(3).max_amplitude;
            beta_f4  = peak_data(i).peaks.HighBeta.(win).channel(4).max_amplitude;
            
        end
        
        mean_theta = (theta_af3 + theta_af4 + theta_f3 + theta_f4) / 4;
        mean_beta = (beta_af3 + beta_af4 + beta_f3 + beta_f4) / 4;
        
        num_trials = length(mean_theta);
        
        % Check if enough trials exist
        if num_trials < ANW + ANT1 + ANT2
            warning(['Participant ', participant_id, ' in window ', win, ...
                     ' has only ', num2str(num_trials), ' trials, less than required ', ...
                     num2str(ANW + ANT1 + ANT2), '. Skipping.']);
            continue;
        end
        
        % Compute reference mean and std for W (for z-scores)
        mean_theta_w = mean(mean_theta(segments.W), 'omitnan');
        std_theta_w = std(mean_theta(segments.W), 'omitnan') + 1e-6; % Avoid division by zero
        mean_beta_w = mean(mean_beta(segments.W), 'omitnan');
        std_beta_w = std(mean_beta(segments.W), 'omitnan') + 1e-6;
        
        % Compute TBR_W (dynamic baseline)
        TBR_W = mean_theta_w / mean_beta_w;
        if ~isfinite(TBR_W)
            TBR_W = NaN; % Handle invalid baseline
        end
        
        % Store participant and window info
        if w == 1
            tbr_results(i).id = participant_id;
        end
        tbr_results(i).windows.(win) = struct();
        
        % Compute metrics for each segment
        segment_names = fieldnames(segments);
        for s = 1:length(segment_names)
            seg_name = segment_names{s};
            seg_indices = segments.(seg_name);
            
            % Segment data
            theta_seg = mean_theta(seg_indices);
            beta_seg = mean_beta(seg_indices);
            
            % Compute statistics
            mean_theta_seg = mean(theta_seg, 'omitnan');
            sd_theta_seg = std(theta_seg, 'omitnan');
            min_theta_seg = min(theta_seg, [], 'omitnan');
            max_theta_seg = max(theta_seg, [], 'omitnan');
            
            mean_beta_seg = mean(beta_seg, 'omitnan');
            sd_beta_seg = std(beta_seg, 'omitnan');
            min_beta_seg = min(beta_seg, [], 'omitnan');
            max_beta_seg = max(beta_seg, [], 'omitnan');
            
            % Compute trial-level z-scores relative to W
            z_theta_seg = (theta_seg - mean_theta_w) / std_theta_w;
            z_beta_seg = (beta_seg - mean_beta_w) / std_beta_w;
            
            % Compute TBR
            if strcmp(win, 'N200')
                if strcmp(seg_name, 'W')
                    TBR_seg = TBR_W; % Use raw mean ratio for W
                    mean_z_theta_seg = mean(z_theta_seg, 'omitnan');
                    mean_z_beta_seg = mean(z_beta_seg, 'omitnan');
                    TBR_seg_j = z_theta_seg ./ abs(z_beta_seg);
                else
                    TBR_seg_j = z_theta_seg ./ abs(z_beta_seg);
                    %TBR_seg = mean(TBR_seg_j, 'omitnan');
                    mean_z_theta_seg = mean(z_theta_seg, 'omitnan');
                    mean_z_beta_seg = mean(z_beta_seg, 'omitnan');
                    TBR_seg = mean_z_theta_seg/ mean_z_beta_seg ;
                    if ~isfinite(TBR_seg)
                        TBR_seg = NaN; % Handle NaN/infinite TBR
                    end
                end
            else 
                if strcmp(seg_name, 'W')
                    TBR_seg = abs(TBR_W); % Use raw mean ratio for W
                    TBR_seg_j = z_theta_seg ./ abs(z_beta_seg);
                    mean_z_theta_seg = mean(z_theta_seg, 'omitnan');
                    mean_z_beta_seg = mean(z_beta_seg, 'omitnan');
                    
                else
                    TBR_seg_j = (z_theta_seg ./ z_beta_seg);
                    %TBR_seg = mean(TBR_seg_j, 'omitnan');
                    mean_z_theta_seg = mean(z_theta_seg, 'omitnan');
                    mean_z_beta_seg = mean(z_beta_seg, 'omitnan');
                    TBR_seg = mean_z_theta_seg/ mean_z_beta_seg ;
                    if ~isfinite(TBR_seg)
                        TBR_seg = NaN; % Handle NaN/infinite TBR
                    end
                end
            end 
            
            % Store mean z-scores for T1 and T2 to compute T1T2 TBR_change
            if strcmp(seg_name, 'T1')
                z_theta_t1 = mean(z_theta_seg, 'omitnan');
                z_beta_t1 = mean(z_beta_seg, 'omitnan');
                TBR_seg_j = z_theta_seg ./ abs(z_beta_seg);
                mean_z_theta_seg = mean(z_theta_seg, 'omitnan');
                mean_z_beta_seg = mean(z_beta_seg, 'omitnan');
                TBR_seg = mean_z_theta_seg/ mean_z_beta_seg ;
            elseif strcmp(seg_name, 'T2')
                z_theta_t2 = mean(z_theta_seg, 'omitnan');
                z_beta_t2 = mean(z_beta_seg, 'omitnan');
                TBR_seg_j = z_theta_seg ./ abs(z_beta_seg);
                mean_z_theta_seg = mean(z_theta_seg, 'omitnan');
                mean_z_beta_seg = mean(z_beta_seg, 'omitnan');
                TBR_seg = mean_z_theta_seg/ mean_z_beta_seg ;
            end
            
            % Compute TBR_change
            if strcmp(seg_name, 'W')
                TBR_change = 0; % Baseline, no change
            elseif strcmp(seg_name, 'T1T2')
                % T1T2: Ratio of z-score differences
                TBR_change = abs((z_theta_t2 - z_beta_t2) / (z_theta_t1 - z_beta_t1 + 1e-6));
                if ~isfinite(TBR_change)
                    TBR_change = NaN; % Handle invalid ratio
                end
            else
                % T1 and T2: Percent change relative to TBR_W
                TBR_change = (abs(TBR_seg - TBR_W) / TBR_W) * 100;
                if ~isfinite(TBR_change)
                    TBR_change = NaN; % Handle invalid percent change
                end
            end
            
            % Store results
            tbr_results(i).windows.(win).(seg_name).mean_theta = mean_theta_seg;
            tbr_results(i).windows.(win).(seg_name).sd_theta = sd_theta_seg;
            tbr_results(i).windows.(win).(seg_name).min_theta = min_theta_seg;
            tbr_results(i).windows.(win).(seg_name).max_theta = max_theta_seg;
            tbr_results(i).windows.(win).(seg_name).mean_beta = mean_beta_seg;
            tbr_results(i).windows.(win).(seg_name).sd_beta = sd_beta_seg;
            tbr_results(i).windows.(win).(seg_name).min_beta = min_beta_seg;
            tbr_results(i).windows.(win).(seg_name).max_beta = max_beta_seg;
            tbr_results(i).windows.(win).(seg_name).z_theta = mean_z_theta_seg;
            tbr_results(i).windows.(win).(seg_name).z_beta = mean_z_beta_seg;
            tbr_results(i).windows.(win).(seg_name).TBR = TBR_seg;
            tbr_results(i).windows.(win).(seg_name).TBR_change = TBR_change;
            tbr_results(i).windows.(win).(seg_name).TBR_trials = TBR_seg_j;
        end
    end
end

% Step 2: Create Segment-Based Table for All Windows with Changes
% Initialize table
tbr_table = cell(num_participants * length(tbr_windows) * length(segment_names) + 1, 14);
tbr_table(1, :) = {'Participant', 'Window', 'Segment', ...
                   'Mean Theta', 'SD Theta', 'Min Theta', 'Max Theta', ...
                   'Mean Beta', 'SD Beta', 'Min Beta', 'Max Beta', ...
                   'meanTheta(z)','meanBeta(z)' 'TBR'};

row = 2;
for w = 1:length(tbr_windows)
    win = tbr_windows{w};
    for i = 1:num_participants
        if ~isfield(tbr_results, 'id') || isempty(tbr_results(i).id) || ...
           ~isfield(tbr_results(i).windows, win)
            continue; % Skip if participant/window was excluded
        end
        participant_id = tbr_results(i).id;
        
        % Fill table for all segments
        for s = 1:length(segment_names)
            seg_name = segment_names{s};
            seg_data = tbr_results(i).windows.(win).(seg_name);
            
            tbr_table(row, :) = {participant_id, win, seg_name, ...
                                 seg_data.mean_theta, seg_data.sd_theta, ...
                                 seg_data.min_theta, seg_data.max_theta, ...
                                 seg_data.mean_beta, seg_data.sd_beta, ...
                                 seg_data.min_beta, seg_data.max_beta, ...
                                 seg_data.z_theta,seg_data.z_beta seg_data.TBR}; %...
                                 %seg_data.TBR_change};
            row = row + 1;
        end
    end
end

% Convert to table, excluding empty rows if any
tbr_table = cell2table(tbr_table(2:row-1, :), 'VariableNames', tbr_table(1, :));

% Display and save the table
disp('Segment-Based Theta/Beta Ratio (TBR) Table with Changes for All Windows:');
disp(tbr_table);
writetable(tbr_table, 'Segment_TBR_N200_Table_with_Changes.csv');

%% P200 - Exploratory Table
% Define trial segments (consistent with FAA)
ANW = 240;    % W: trials 1 to 240
ANT1 = 30;    % T1: trials 241 to 270
ANT2 = 30;    % T2: trials 271 to 300
segments = struct(...
    'W', 1:ANW, ...
    'T1', ANW+1:ANW+ANT1, ...
    'T2', ANW+ANT1+1:ANW+ANT1+ANT2, ...
    'T1T2', ANW+1:ANW+ANT1+ANT2);

% Initialize storage for P200 exploratory results
p200_exploratory_results = struct();

% Step 1: Compute Theta power and z-scores for P200 window
for i = 1:num_participants
    participant_id = participants{i};
    
    % Extract Theta Peak Power for P200 window from 4 frontal channels (FC5, FC6)
    theta_p200 = (peak_data(i).peaks.Theta.P200.channel(5).max_amplitude + ... % FC5
                  peak_data(i).peaks.Theta.P200.channel(6).max_amplitude ... % FC6
                  )/2; % Mean across channels
                  %peak_data(i).peaks.Theta.P200.channel(3).max_amplitude + ... % F3
                  %peak_data(i).peaks.Theta.P200.channel(4).max_amplitude) / 4; % Mean across channels
    
    num_trials = length(theta_p200);
    
    % Check if enough trials exist
    if num_trials < ANW + ANT1 + ANT2
        warning(['Participant ', participant_id, ' in P200 window has only ', ...
                 num2str(num_trials), ' trials, less than required ', ...
                 num2str(ANW + ANT1 + ANT2), '. Skipping.']);
        continue;
    end
    
    % Compute reference mean and std for W (control)
    mean_theta_w = mean(theta_p200(segments.W), 'omitnan');
    std_theta_w = std(theta_p200(segments.W), 'omitnan') + 1e-6; % Avoid division by zero
    
    % Store participant info
    p200_exploratory_results(i).id = participant_id;
    
    % Compute metrics for each segment
    segment_names = fieldnames(segments);
    for s = 1:length(segment_names)
        seg_name = segment_names{s};
        seg_indices = segments.(seg_name);
        
        % Segment data
        theta_seg = theta_p200(seg_indices);
        
        % Compute trial-level z-scores relative to W
        z_theta_seg = (theta_seg - mean_theta_w) / std_theta_w;
        z_theta_seg = squeeze(z_theta_seg); % Remove singleton dimensions
        if size(z_theta_seg, 1) > 1 % Ensure row vector
            z_theta_seg = z_theta_seg';
        end
        
        % Compute statistics
        mean_z_theta = mean(z_theta_seg, 'omitnan');
        sd_z_theta = std(z_theta_seg, 'omitnan');
        min_z_theta = min(z_theta_seg, [], 'omitnan');
        max_z_theta = max(z_theta_seg, [], 'omitnan');
        
        if ~isfinite(mean_z_theta)
            mean_z_theta = NaN;
            sd_z_theta = NaN;
            min_z_theta = NaN;
            max_z_theta = NaN;
        end
        
        % Store results
        p200_exploratory_results(i).(seg_name).mean_z_theta = mean_z_theta;
        p200_exploratory_results(i).(seg_name).sd_z_theta = sd_z_theta;
        p200_exploratory_results(i).(seg_name).min_z_theta = min_z_theta;
        p200_exploratory_results(i).(seg_name).max_z_theta = max_z_theta;
        p200_exploratory_results(i).(seg_name).z_theta_trials = z_theta_seg; % Per-trial z-scores
    end
end

% Step 2: Create P200 Exploratory Table
% Initialize table with Segment renamed to Metric
p200_table = cell(num_participants * length(segment_names) + 1, 6);
p200_table(1, :) = {'Participant', 'Segment', 'Mean Theta(z)', 'SD', 'Min', 'Max'};

segment_names = {'W', 'T1', 'T2', 'T1T2'}; % Match your requested labels (T1+T2 as T1T2)
row = 2;
for i = 1:num_participants
    if ~isfield(p200_exploratory_results, 'id') || isempty(p200_exploratory_results(i).id)
        for s = 1:length(segment_names)
            seg_name = segment_names{s};
            p200_table(row, :) = {participants{i}, seg_name, NaN, NaN, NaN, NaN};
            row = row + 1;
        end
        continue;
    end
    participant_id = p200_exploratory_results(i).id;
    
    % Fill table for each segment
    for s = 1:length(segment_names)
        seg_name = segment_names{s};
        seg_data = p200_exploratory_results(i).(seg_name);
        
        p200_table(row, :) = {participant_id, seg_name, ...
                              seg_data.mean_z_theta, seg_data.sd_z_theta, ...
                              seg_data.min_z_theta, seg_data.max_z_theta};
        row = row + 1;
    end
end

% Convert to table
p200_table = cell2table(p200_table(2:end, :), 'VariableNames', p200_table(1, :));

% Display and save the table
disp('P200 - Exploratory Table:');
disp(p200_table);
writetable(p200_table, 'P200_Exploratory_Table.csv');

%% N400 - Exploratory Table
% Define trial segments (consistent with FAA and High Beta)
ANW = 240;    % W: trials 1 to 240
ANT1 = 30;    % T1: trials 241 to 270
ANT2 = 30;    % T2: trials 271 to 300
segments = struct(...
    'W', 1:ANW, ...
    'T1', ANW+1:ANW+ANT1, ...
    'T2', ANW+ANT1+1:ANW+ANT1+ANT2, ...
    'T1T2', ANW+1:ANW+ANT1+ANT2);

% Initialize storage for N400 exploratory results
n400_exploratory_results = struct();

% Step 1: Compute Theta power and z-scores for N400 window
for i = 1:num_participants
    participant_id = participants{i};
    
    % Extract Theta Peak Power for N400 window from 4 frontal channels (FC5, FC6)
    theta_n400 = (peak_data(i).peaks.Theta.N400.channel(5).min_amplitude + ... % FC53
                  peak_data(i).peaks.Theta.N400.channel(6).min_amplitude  ... % FC6
                  )/2; % Mean across channels
                  %peak_data(i).peaks.Theta.N400.channel(3).max_amplitude + ... % F3
                  %peak_data(i).peaks.Theta.N400.channel(4).max_amplitude) / 4; % Mean across channels
    
    num_trials = length(theta_n400);
    
    % Check if enough trials exist
    if num_trials < ANW + ANT1 + ANT2
        warning(['Participant ', participant_id, ' in N400 window has only ', ...
                 num2str(num_trials), ' trials, less than required ', ...
                 num2str(ANW + ANT1 + ANT2), '. Skipping.']);
        continue;
    end
    
    % Compute reference mean and std for W (control)
    mean_theta_w = mean(theta_n400(segments.W), 'omitnan');
    std_theta_w = std(theta_n400(segments.W), 'omitnan') + 1e-6; % Avoid division by zero
    
    % Store participant info
    n400_exploratory_results(i).id = participant_id;
    
    % Compute metrics for each segment
    segment_names = fieldnames(segments);
    for s = 1:length(segment_names)
        seg_name = segment_names{s};
        seg_indices = segments.(seg_name);
        
        % Segment data
        theta_seg = theta_n400(seg_indices);
        
        % Compute trial-level z-scores relative to W
        z_theta_seg = (theta_seg - mean_theta_w) / std_theta_w;
        z_theta_seg = squeeze(z_theta_seg); % Remove singleton dimensions
        if size(z_theta_seg, 1) > 1 % Ensure row vector
            z_theta_seg = z_theta_seg';
        end
        
        % Compute statistics
        mean_z_theta = mean(z_theta_seg, 'omitnan');
        sd_z_theta = std(z_theta_seg, 'omitnan');
        min_z_theta = min(z_theta_seg, [], 'omitnan');
        max_z_theta = max(z_theta_seg, [], 'omitnan');
        
        if ~isfinite(mean_z_theta)
            mean_z_theta = NaN;
            sd_z_theta = NaN;
            min_z_theta = NaN;
            max_z_theta = NaN;
        end
        
        % Store results
        n400_exploratory_results(i).(seg_name).mean_z_theta = mean_z_theta;
        n400_exploratory_results(i).(seg_name).sd_z_theta = sd_z_theta;
        n400_exploratory_results(i).(seg_name).min_z_theta = min_z_theta;
        n400_exploratory_results(i).(seg_name).max_z_theta = max_z_theta;
        n400_exploratory_results(i).(seg_name).z_theta_trials = z_theta_seg; % Per-trial z-scores
    end
end

% Step 2: Create N400 Exploratory Table
% Initialize table with Metric column
n400_table = cell(num_participants * length(segment_names) + 1, 6);
n400_table(1, :) = {'Participant', 'Segment', 'Mean Theta(z)', 'SD', 'Min', 'Max'};

segment_names = {'W', 'T1', 'T2', 'T1T2'}; % Match your requested labels (T1+T2 as T1T2)
row = 2;
for i = 1:num_participants
    if ~isfield(n400_exploratory_results, 'id') || isempty(n400_exploratory_results(i).id)
        for s = 1:length(segment_names)
            seg_name = segment_names{s};
            n400_table(row, :) = {participants{i}, seg_name, NaN, NaN, NaN, NaN};
            row = row + 1;
        end
        continue;
    end
    participant_id = n400_exploratory_results(i).id;
    
    % Fill table for each segment
    for s = 1:length(segment_names)
        seg_name = segment_names{s};
        seg_data = n400_exploratory_results(i).(seg_name);
        
        n400_table(row, :) = {participant_id, seg_name, ...
                              seg_data.mean_z_theta, seg_data.sd_z_theta, ...
                              seg_data.min_z_theta, seg_data.max_z_theta};
        row = row + 1;
    end
end

% Convert to table
n400_table = cell2table(n400_table(2:end, :), 'VariableNames', n400_table(1, :));

% Display and save the table
disp('N400 - Exploratory Table:');
disp(n400_table);
writetable(n400_table, 'N400_Exploratory_Table.csv');

%% H1 - T1 vs. T2 Analysis and Visuals (P300 Window Only)
% Define trial segments (consistent with previous codes)
ANW = 240;    % W: trials 1 to 240
ANT1 = 30;    % T1: trials 241 to 270
ANT2 = 30;    % T2: trials 271 to 300
segments = struct(...
    'W', 1:ANW, ...
    'T1', ANW+1:ANW+ANT1, ...
    'T2', ANW+ANT1+1:ANW+ANT1+ANT2, ...
    'T1T2', ANW+1:ANW+ANT1+ANT2);

% Step 1: Extract T1 and T2 mean z-scores and SDs from precalculated results
t1_t2_results = struct();
for i = 1:num_participants
    participant_id = participants{i};
    
    % Extract BAR data for P300 window
    if isfield(segment_bar_results(i), 'windows') && isfield(segment_bar_results(i).windows, 'P300')
        % T1
        bar_t1 = segment_bar_results(i).windows.P300.T1.BAR;
        beta_t1 = segment_bar_results(i).windows.P300.T1.mean_beta;
        alpha_t1 = segment_bar_results(i).windows.P300.T1.mean_alpha;
        sd_beta_t1 = segment_bar_results(i).windows.P300.T1.sd_beta;
        sd_alpha_t1 = segment_bar_results(i).windows.P300.T1.sd_alpha;
        bar_t1_sd = abs(bar_t1) * sqrt((sd_beta_t1 / (beta_t1 + 1e-6))^2 + (sd_alpha_t1 / (alpha_t1 + 1e-6))^2);
        
        % T2
        bar_t2 = segment_bar_results(i).windows.P300.T2.BAR;
        beta_t2 = segment_bar_results(i).windows.P300.T2.mean_beta;
        alpha_t2 = segment_bar_results(i).windows.P300.T2.mean_alpha;
        sd_beta_t2 = segment_bar_results(i).windows.P300.T2.sd_beta;
        sd_alpha_t2 = segment_bar_results(i).windows.P300.T2.sd_alpha;
        bar_t2_sd = abs(bar_t2) * sqrt((sd_beta_t2 / (beta_t2 + 1e-6))^2 + (sd_alpha_t2 / (alpha_t2 + 1e-6))^2);
    else
        bar_t1 = NaN; bar_t2 = NaN;
        bar_t1_sd = NaN; bar_t2_sd = NaN;
    end
    
    % Extract High Beta data for P300 window
    if isfield(hb_results(i), 'windows') && isfield(hb_results(i).windows, 'P300')
        hb_t1 = hb_results(i).windows.P300.T1.z_hb; % Mean z-score for T1
        hb_t1_sd = hb_results(i).windows.P300.T1.sd_hb; % Provided SD for T1
        hb_t2 = hb_results(i).windows.P300.T2.z_hb; % Mean z-score for T2
        hb_t2_sd = hb_results(i).windows.P300.T2.sd_hb; % Provided SD for T2
    else
        hb_t1 = NaN; hb_t2 = NaN;
        hb_t1_sd = NaN; hb_t2_sd = NaN;
    end
    
    % Store results
    t1_t2_results(i).id = participant_id;
    t1_t2_results(i).bar_t1 = bar_t1;
    t1_t2_results(i).bar_t2 = bar_t2;
    t1_t2_results(i).bar_t1_sd = bar_t1_sd;
    t1_t2_results(i).bar_t2_sd = bar_t2_sd;
    t1_t2_results(i).hb_t1 = hb_t1;
    t1_t2_results(i).hb_t2 = hb_t2;
    t1_t2_results(i).hb_t1_sd = hb_t1_sd;
    t1_t2_results(i).hb_t2_sd = hb_t2_sd;
end

% Step 2: Compute Aggregates for Visuals
bar_t1_means = [t1_t2_results(:).bar_t1];
bar_t2_means = [t1_t2_results(:).bar_t2];
bar_t1_sds = [t1_t2_results(:).bar_t1_sd];
bar_t2_sds = [t1_t2_results(:).bar_t2_sd];

hb_t1_means = [t1_t2_results(:).hb_t1];
hb_t2_means = [t1_t2_results(:).hb_t2];
hb_t1_sds = [t1_t2_results(:).hb_t1_sd];
hb_t2_sds = [t1_t2_results(:).hb_t2_sd];

% Step 3: Side-by-Side Bar Chart (All Participants)
figure('Name', 'H1 Side-by-Side Bar Chart (P300)');
x = 1:num_participants; % 1:5 for 5 participants
bar_width = 3;

% BAR group (5 participants, T1 and T2)
bar_positions_bar = [x-0.2; x+0.2]'; % Shift T1 left, T2 right for BAR
b_bar = bar(bar_positions_bar, [bar_t1_means', bar_t2_means'], bar_width);
b_bar(1).FaceColor = [0.4660, 0.6740, 0.1880]; % T1: Green
b_bar(2).FaceColor = [0.9290, 0.6940, 0.1250]; % T2: Yellow

hold on;

% High Beta group (5 participants, T1 and T2, offset by num_participants + gap)
bar_positions_hb = [x-0.2 + num_participants + 1; x+0.2 + num_participants + 1]';
b_hb = bar(bar_positions_hb, [hb_t1_means', hb_t2_means'], bar_width);
b_hb(1).FaceColor = [0.4660, 0.6740, 0.1880]; % T1: Green
b_hb(2).FaceColor = [0.9290, 0.6940, 0.1250]; % T2: Yellow

% Error bars for BAR
errorbar(bar_positions_bar(:, 1), bar_t1_means, bar_t1_sds, 'k',...
    'LineStyle', 'none', 'LineWidth', 1);
errorbar(bar_positions_bar(:, 2), bar_t2_means, bar_t2_sds, 'k',...
    'LineStyle', 'none', 'LineWidth', 1);

% Error bars for High Beta
errorbar(bar_positions_hb(:, 1), hb_t1_means, hb_t1_sds, 'k',...
    'LineStyle', 'none', 'LineWidth', 1);
errorbar(bar_positions_hb(:, 2), hb_t2_means, hb_t2_sds, 'k',...
    'LineStyle', 'none', 'LineWidth', 1);

xticks([1:num_participants, (num_participants+2):(2*num_participants+1)]);
xticklabels([arrayfun(@(i) sprintf('P%d BAR', i), 1:num_participants, 'UniformOutput', false), ...
             arrayfun(@(i) sprintf('P%d HB', i), 1:num_participants, 'UniformOutput', false)]);
xtickangle(45);
ylabel('Mean Z-Score');
title('T1 vs. T2 Mean Z-Scores by Participant (P300 Window) ± 1 SD');
legend({'T1', 'T2'}, 'Location', 'best');
grid on;
hold off;
saveas(gcf, 'H1_Side_by_Side_Bar_Chart_All_Participants.png');

% Step 4: Paired Dot Plot
figure('Name', 'H1 Paired Dot Plot (P300)');
colors = lines(num_participants); % Distinct colors for participants

% BAR Subplot
subplot(2, 1, 1);
hold on;
for i = 1:num_participants
    plot([1, 2], [bar_t1_means(i), bar_t2_means(i)], 'o-', ...
         'Color', colors(i,:), 'MarkerFaceColor', colors(i,:), 'LineWidth', 1);
end
xticks([1, 2]);
xticklabels({'T1', 'T2'});
ylabel('Z-Score');
title('BAR - T1 vs. T2 (P300)');
legend(participants, 'Location', 'bestoutside', 'Box', 'off');
grid on;

% High Beta Subplot
subplot(2, 1, 2);
hold on;
for i = 1:num_participants
    plot([1, 2], [hb_t1_means(i), hb_t2_means(i)], 'o-', ...
         'Color', colors(i,:), 'MarkerFaceColor', colors(i,:), 'LineWidth', 1);
end
xticks([1, 2]);
xticklabels({'T1', 'T2'});
ylabel('Z-Score');
title('High Beta - T1 vs. T2 (P300)');
legend(participants, 'Location', 'bestoutside', 'Box', 'off');
grid on;
hold off;
saveas(gcf, 'H1_Paired_Dot_Plot.png');


%% Wilcoxon Signed-Rank Tests
% Step 1: Initialize variables and parameters
num_participants = 5;
window = 'P300';
T1_range = 241:270; % 30 values per participant for T1
T2_range = 271:300; % 30 values per participant for T2

% Initialize storage for T1 and T2 data
hb_zscoredT1 = zeros(num_participants, 30);
hb_zscoredT2 = zeros(num_participants, 30);
bar_dataT1 = zeros(num_participants, 30);
bar_dataT2 = zeros(num_participants, 30);

% Step 2: Process each participant's data
for i = 1:num_participants
    % High Beta Processing
    % Extract raw data
    hb_data = hb_results(i).windows.P300.acrosschannel.mean_hb;
    hb_data_T1 = hb_data(T1_range);
    hb_data_T2 = hb_data(T2_range);
    
    % Z-score using W stimuli mean and SD
    w_mean = hb_results(i).windows.P300.W.mean_hb;
    w_sd = hb_results(i).windows.P300.W.sd_hb;
    hb_zscoredT1(i,:) = (hb_data_T1 - w_mean) / w_sd;
    hb_zscoredT2(i,:) = (hb_data_T2 - w_mean) / w_sd;
    
    % BAR Processing
    % Extract Alpha and Beta data
    alpha_data = segment_bar_results(i).windows.P300.acrosschannel.mean_alpha;
    beta_data = segment_bar_results(i).windows.P300.acrosschannel.mean_beta;
    alpha_data_T1 = alpha_data(T1_range);
    alpha_data_T2 = alpha_data(T2_range);
    beta_data_T1 = beta_data(T1_range);
    beta_data_T2 = beta_data(T2_range);
    
    % Z-score using W stimuli mean and SD for both Alpha and Beta
    w_alpha_mean = segment_bar_results(i).windows.P300.W.mean_alpha;
    w_alpha_sd = segment_bar_results(i).windows.P300.W.sd_alpha;
    alpha_zscoredT1 = (alpha_data_T1 - w_alpha_mean) / w_alpha_sd;
    alpha_zscoredT2 = (alpha_data_T2 - w_alpha_mean) / w_alpha_sd;
    
    w_beta_mean = segment_bar_results(i).windows.P300.W.mean_beta;
    w_beta_sd = segment_bar_results(i).windows.P300.W.sd_beta;
    beta_zscoredT1 = (beta_data_T1 - w_beta_mean) / w_beta_sd;
    beta_zscoredT2 = (beta_data_T2 - w_beta_mean) / w_beta_sd;
    
    % Calculate BAR as Alpha/Beta
    bar_dataT1(i,:) = alpha_zscoredT1 ./ beta_zscoredT1;
    bar_dataT2(i,:) = alpha_zscoredT2 ./ beta_zscoredT2;
end


% Step 3: Perform Wilcoxon Signed-Rank Tests and collect results

  

results = cell(2, 7); % [Metric, T1 Mean, T2 Mean, Z, p-value, Effect Size, Interpretation]
row = 1;

[bar_p, ~, bar_stats] = signrank(bar_dataT1(:), bar_dataT2(:));
    if isfield(bar_stats, 'zval')
        bar_z = bar_stats.zval;
    else
        bar_z = NaN;
    end
    bar_r = abs(bar_z) / sqrt(150);
    bar_interp = 'Small';
    if bar_r >= 0.5
        bar_interp = 'Large';
    elseif bar_r >= 0.3
        bar_interp = 'Medium';
    end
    bar_t1_mean = mean(bar_dataT1(:));
    bar_t2_mean = mean(bar_dataT2(:));
 % High Beta
[hb_p, ~, hb_stats] = signrank(hb_zscoredT1(:), hb_zscoredT2(:));
    if isfield(hb_stats, 'zval')
        hb_z = hb_stats.zval;
    else
        hb_z = NaN;
    end
    hb_r = abs(hb_z) / sqrt(150);
    hb_interp = 'Small';
    if hb_r >= 0.5
        hb_interp = 'Large';
    elseif hb_r >= 0.3
        hb_interp = 'Medium';
    end
    hb_t1_mean = mean(hb_zscoredT1(:));
    hb_t2_mean = mean(hb_zscoredT2(:));
    % Store results
    results(row, :) = {'BAR', bar_t1_mean, bar_t2_mean, bar_z, bar_p, bar_r, bar_interp};
    row = row + 1;
    
    results(row, :) = {'High Beta', hb_t1_mean, hb_t2_mean, hb_z, hb_p, hb_r, hb_interp};
    row = row + 1;


% Step 4: Create and display the results table
results_table = cell2table(results, ...
    'VariableNames', {'Metric', 'T1 Mean', 'T2 Mean', 'Z-score', 'p-value', 'Effect Size (r)', 'Effect Size Interpretation'});

disp('T1 vs. T2 for H1:');
disp(results_table);

% Step 5: Save the table to a CSV file
writetable(results_table, 'H1_T1_vs_T2_Results.csv');

% Compute averages and SDs across participants

mean_bar_T1 = mean(bar_t1_means);
mean_bar_T2 = mean(bar_t2_means);
mean_hb_T1 = mean(hb_t1_means);
mean_hb_T2 = mean(hb_t2_means);
Sample_size = 30*ones(1,5);
T1_BAR_pooledSTD = sum((Sample_size-1).*bar_t1_sds.^2)./sum(Sample_size-1);
T2_BAR_pooledSTD = sum((Sample_size-1).*bar_t2_sds.^2)./sum(Sample_size-1);
T1_HB_pooledSTD = sum((Sample_size-1).*hb_t1_sds.^2)./sum(Sample_size-1);
T2_HB_pooledSTD = sum((Sample_size-1).*hb_t2_sds.^2)./sum(Sample_size-1);

sd_bar_T1 = sqrt(T1_BAR_pooledSTD);
sd_bar_T2 = sqrt(T2_BAR_pooledSTD);
sd_hb_T1 = sqrt(T1_HB_pooledSTD);
sd_hb_T2 = sqrt(T2_HB_pooledSTD);

% Create figure and plot averaged bar chart
figure('Name', 'Averaged T1 and T2 Bar Plots Grouped by BAR and High Beta');
bar_data = [mean_bar_T1, mean_bar_T2; mean_hb_T1, mean_hb_T2];
b = bar(bar_data);

% Set colors
b(1).FaceColor = [0.4660, 0.6740, 0.1880]; % T1: Green
b(2).FaceColor = [0.9290, 0.6940, 0.1250]; % T2: Yellow

hold on;
% Add error bars
errorbar(b(1).XEndPoints, [mean_bar_T1, mean_hb_T1], [sd_bar_T1, sd_hb_T1], ...
    'k', 'LineStyle', 'none', 'LineWidth', 1);
errorbar(b(2).XEndPoints, [mean_bar_T2, mean_hb_T2], [sd_bar_T2, sd_hb_T2], ...
    'k', 'LineStyle', 'none', 'LineWidth', 1);

% Customize plot
xticks([1, 2]);
xticklabels({'BAR', 'High Beta'});
ylabel('Mean Z-Score');
title('Averaged T1 and T2 Mean Z-Scores ± 1 SD (P300 Window)');
legend({'T1', 'T2'}, 'Location', 'best');
grid on;
hold off;

% Save figure
saveas(gcf, 'Averaged_T1_T2_Bar_Plot.png');


%% H2 - Peak and Z-Score Graphs for All Participants and Y* Group with Baseline, Statistics, and Cohen's d
% Define parameters
T1_trials = segments.T1;      % Trial indices for T1 condition (241:270)
W_trials = segments.W;        % Trial indices for W condition (baseline)
chanlocs = all_eeg(1).EEG.chanlocs;  % Channel locations
num_channels = length(chanlocs);     % Number of EEG channels (assumed 8: AF3, AF4, F3, F4, FC5, FC6, P8, F8)
num_participants = length(participants);  % Number of participants (5)
all_windows = {'W', 'P200', 'N200', 'P300', 'LPP'};  % ERP windows including W as baseline
metrics = {'FAA', 'BAR', 'TBR', 'HB'};
y_group_indices = [2, 3, 4, 5];  % Indices for Y* group: PT2YY, PT3YY, PT4YY, PT7YN

% Find F3 and F4 channel indices (for FAA)
f3_idx = find(strcmp({chanlocs.labels}, 'F3'));
f4_idx = find(strcmp({chanlocs.labels}, 'F4'));
if isempty(f3_idx) || isempty(f4_idx)
    error('F3 or F4 channel not found in chanlocs.');
end

% Define colors for channels (for peak graphs with multiple channels)
channel_colors = lines(num_channels);  % MATLAB's 'lines' colormap for distinct colors

% Define time mapping for ERP windows (ms)
% Baseline (W): assume 0–100 ms with midpoint = 50 ms
% P200: 150–250 ms, midpoint = 200 ms
% N200: 180–240 ms, midpoint = 210 ms
% P300: 300–400 ms, midpoint = 350 ms
% LPP: 400–800 ms, midpoint = 600 ms
timeCenters = [50, 175, 250, 350, 510];  % one value per window (W, P200, N200, P300, LPP)
% For shading, we use the ERP windows (excluding baseline, index 1)
windowBounds = [150 200;   % P200 bounds
                210 290;   % N200 bounds
                310 390;   % P300 bounds
                420 600];  % LPP bounds
yMin_default = -1;
yMax_default = 2;

% Initialize arrays to store z-scores for statistical analysis
faa_zscores_all = zeros(num_participants, length(all_windows));  % [participants, windows]
bar_zscores_all = zeros(num_participants, length(all_windows));
tbr_zscores_all = zeros(num_participants, length(all_windows));
hb_zscores_all = zeros(num_participants, length(all_windows));

%% --- Peak Graphs for Individuals ---
for i = 1:num_participants
    participant_id = participants{i};
    
    %% FAA Peak Graph (F3 vs F4) with W baseline
    figure('Name', ['FAA Peaks - ' participant_id], 'Position', [100, 100, 800, 600]);
    hold on;
    f3_peaks = zeros(1, length(all_windows));
    f4_peaks = zeros(1, length(all_windows));
    % W segment (baseline)
    alpha_peaks_w = zeros(num_channels, length(peak_data(i).peaks.Alpha.P200.channel(1).max_amplitude));
    for ch = 1:num_channels
        alpha_peaks_w(ch, :) = peak_data(i).peaks.Alpha.P200.channel(ch).max_amplitude;
    end
    f3_peaks(1) = mean(alpha_peaks_w(f3_idx, W_trials), 'omitnan');
    f4_peaks(1) = mean(alpha_peaks_w(f4_idx, W_trials), 'omitnan');
    % T1 segments
    for w = 2:length(all_windows)
        win = all_windows{w};
        alpha_peaks = zeros(num_channels, length(peak_data(i).peaks.Alpha.(win).channel(1).max_amplitude));
        for ch = 1:num_channels
            alpha_peaks(ch, :) = peak_data(i).peaks.Alpha.(win).channel(ch).max_amplitude;
        end
        f3_peaks(w) = mean(alpha_peaks(f3_idx, T1_trials), 'omitnan');
        f4_peaks(w) = mean(alpha_peaks(f4_idx, T1_trials), 'omitnan');
    end
    plot(timeCenters, f3_peaks, '-o', 'Color', [1, 0.5, 0], 'LineWidth', 2, 'DisplayName', 'F3');
    plot(timeCenters, f4_peaks, '-o', 'Color', [0, 0.5, 1], 'LineWidth', 2, 'DisplayName', 'F4');
    % Shade ERP windows only (skip baseline: index 1)
    for w = 2:length(all_windows)
        boundsIndex = w - 1;
        tStart = windowBounds(boundsIndex, 1);
        tEnd   = windowBounds(boundsIndex, 2);
        h = fill([tStart, tEnd, tEnd, tStart], [yMin_default, yMin_default, yMax_default, yMax_default], [0.8, 0.9, 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
    end
    xticks(timeCenters);
    xticklabels(all_windows);
    xlabel('Time (ms)');
    ylabel('Peak Z-Score (μV)');
    title(['FAA Peaks Across ERP Components (F3 vs F4) - ' participant_id]);
    xlim([0, 800]);
    ylim([yMin_default, yMax_default]);
    grid on;
    legend('Location','best');
    hold off;
    saveas(gcf, ['H2_PeakGraph_FAA_' participant_id '.png']);
    close(gcf);
    
    %% BAR Peak Graph (per channel) with W baseline
    figure('Name', ['BAR Peaks - ' participant_id], 'Position', [100, 100, 800, 600]);
    hold on;
    bar_peaks = zeros(num_channels, length(all_windows));
    % W segment
    alpha_peaks_w = zeros(num_channels, length(peak_data(i).peaks.Alpha.P200.channel(1).max_amplitude));
    beta_peaks_w = zeros(num_channels, length(peak_data(i).peaks.Beta.P200.channel(1).max_amplitude));
    for ch = 1:num_channels
        alpha_peaks_w(ch, :) = peak_data(i).peaks.Alpha.P200.channel(ch).max_amplitude;
        beta_peaks_w(ch, :) = peak_data(i).peaks.Beta.P200.channel(ch).max_amplitude;
    end
    bar_per_trial_w = beta_peaks_w(:, W_trials) ./ alpha_peaks_w(:, W_trials);
    bar_peaks(:, 1) = mean(bar_per_trial_w, 2, 'omitnan');
    % T1 segments
    for w = 2:length(all_windows)
        win = all_windows{w};
        alpha_peaks = zeros(num_channels, length(peak_data(i).peaks.Alpha.(win).channel(1).max_amplitude));
        beta_peaks  = zeros(num_channels, length(peak_data(i).peaks.Beta.(win).channel(1).max_amplitude));
        for ch = 1:num_channels
            alpha_peaks(ch, :) = peak_data(i).peaks.Alpha.(win).channel(ch).max_amplitude;
            beta_peaks(ch, :)  = peak_data(i).peaks.Beta.(win).channel(ch).max_amplitude;
        end
        bar_per_trial = beta_peaks(:, T1_trials) ./ alpha_peaks(:, T1_trials);
        bar_peaks(:, w) = mean(bar_per_trial, 2, 'omitnan');
    end
    for ch = 1:num_channels
        plot(timeCenters, bar_peaks(ch, :), '-o', 'Color', channel_colors(ch, :), ...
            'LineWidth', 2, 'DisplayName', chanlocs(ch).labels);
    end
    % Shade ERP windows (skip baseline)
    for w = 2:length(all_windows)
        boundsIndex = w - 1;
        tStart = windowBounds(boundsIndex, 1);
        tEnd   = windowBounds(boundsIndex, 2);
        h = fill([tStart, tEnd, tEnd, tStart], [min(bar_peaks(:))-0.1, min(bar_peaks(:))-0.1, max(bar_peaks(:))+0.1, max(bar_peaks(:))+0.1], ...
            [0.8, 0.9, 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
    end
    xticks(timeCenters);
    xticklabels(all_windows);
    xlabel('Time (ms)');
    ylabel('Peak BAR');
    title(['BAR Peaks Across ERP Windows - ' participant_id]);
    grid on;
    legend('Location','best');
    hold off;
    saveas(gcf, ['H2_PeakGraph_BAR_' participant_id '.png']);
    close(gcf);
    
    %% TBR Peak Graph (per channel) with W baseline
    figure('Name', ['TBR Peaks - ' participant_id], 'Position', [100, 100, 800, 600]);
    hold on;
    tbr_peaks = zeros(num_channels, length(all_windows));
    % W segment
    theta_peaks_w = zeros(num_channels, length(peak_data(i).peaks.Theta.P200.channel(1).max_amplitude));
    beta_peaks_w  = zeros(num_channels, length(peak_data(i).peaks.Beta.P200.channel(1).max_amplitude));
    for ch = 1:num_channels
        theta_peaks_w(ch, :) = peak_data(i).peaks.Theta.P200.channel(ch).max_amplitude;
        beta_peaks_w(ch, :)  = peak_data(i).peaks.Beta.P200.channel(ch).max_amplitude;
    end
    tbr_per_trial_w = theta_peaks_w(:, W_trials) ./ beta_peaks_w(:, W_trials);
    tbr_peaks(:, 1) = mean(tbr_per_trial_w, 2, 'omitnan');
    % T1 segments
    for w = 2:length(all_windows)
        win = all_windows{w};
        theta_peaks = zeros(num_channels, length(peak_data(i).peaks.Theta.(win).channel(1).max_amplitude));
        beta_peaks  = zeros(num_channels, length(peak_data(i).peaks.Beta.(win).channel(1).max_amplitude));
        for ch = 1:num_channels
            theta_peaks(ch, :) = peak_data(i).peaks.Theta.(win).channel(ch).max_amplitude;
            beta_peaks(ch, :)  = peak_data(i).peaks.Beta.(win).channel(ch).max_amplitude;
        end
        tbr_per_trial = theta_peaks(:, T1_trials) ./ beta_peaks(:, T1_trials);
        tbr_peaks(:, w) = mean(tbr_per_trial, 2, 'omitnan');
    end
    for ch = 1:num_channels
        plot(timeCenters, tbr_peaks(ch, :), '-o', 'Color', channel_colors(ch, :), ...
            'LineWidth', 2, 'DisplayName', chanlocs(ch).labels);
    end
    % Shade ERP windows (skip baseline)
    for w = 2:length(all_windows)
        boundsIndex = w - 1;
        tStart = windowBounds(boundsIndex, 1);
        tEnd   = windowBounds(boundsIndex, 2);
        h = fill([tStart, tEnd, tEnd, tStart], [min(tbr_peaks(:))-0.1, min(tbr_peaks(:))-0.1, max(tbr_peaks(:))+0.1, max(tbr_peaks(:))+0.1], ...
            [0.8, 0.9, 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
    end
    xticks(timeCenters);
    xticklabels(all_windows);
    xlabel('Time (ms)');
    ylabel('Peak TBR');
    title(['TBR Peaks Across ERP Windows - ' participant_id]);
    grid on;
    legend('Location','best');
    hold off;
    saveas(gcf, ['H2_PeakGraph_TBR_' participant_id '.png']);
    close(gcf);
    
    %% HB Peak Graph (per channel) with W baseline
    figure('Name', ['HB Peaks - ' participant_id], 'Position', [100, 100, 800, 600]);
    hold on;
    hb_peaks = zeros(num_channels, length(all_windows));
    % W segment
    hb_peaks_w = zeros(num_channels, length(peak_data(i).peaks.HighBeta.P200.channel(1).max_amplitude));
    for ch = 1:num_channels
        hb_peaks_w(ch, :) = peak_data(i).peaks.HighBeta.P200.channel(ch).max_amplitude;
    end
    hb_peaks(:, 1) = mean(hb_peaks_w(:, W_trials), 2, 'omitnan');
    % T1 segments
    for w = 2:length(all_windows)
        win = all_windows{w};
        hb_peaks_all = zeros(num_channels, length(peak_data(i).peaks.HighBeta.(win).channel(1).max_amplitude));
        for ch = 1:num_channels
            hb_peaks_all(ch, :) = peak_data(i).peaks.HighBeta.(win).channel(ch).max_amplitude;
        end
        hb_peaks(:, w) = mean(hb_peaks_all(:, T1_trials), 2, 'omitnan');
    end
    for ch = 1:num_channels
        plot(timeCenters, hb_peaks(ch, :), '-o', 'Color', channel_colors(ch, :), ...
            'LineWidth', 2, 'DisplayName', chanlocs(ch).labels);
    end
    % Shade ERP windows (skip baseline)
    for w = 2:length(all_windows)
        boundsIndex = w - 1;
        tStart = windowBounds(boundsIndex, 1);
        tEnd   = windowBounds(boundsIndex, 2);
        h = fill([tStart, tEnd, tEnd, tStart], [min(hb_peaks(:))-0.1, min(hb_peaks(:))-0.1, max(hb_peaks(:))+0.1, max(hb_peaks(:))+0.1], ...
            [0.8, 0.9, 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
    end
    xticks(timeCenters);
    xticklabels(all_windows);
    xlabel('Time (ms)');
    ylabel('Peak Z-Score (μV)');
    title(['HB Peaks Across ERP Windows - ' participant_id]);
    grid on;
    legend('Location','best');
    hold off;
    saveas(gcf, ['H2_PeakGraph_HB_' participant_id '.png']);
    close(gcf);
end

%% --- Z-Score Graphs for Individuals ---
for i = 1:num_participants
    participant_id = participants{i};
    
    %% FAA Z-Score Graph with W baseline
    figure('Name', ['FAA Z-Scores - ' participant_id], 'Position', [100, 100, 800, 600]);
    hold on;
    faa_zscores = zeros(1, length(all_windows));
    % W segment
    alpha_peaks_w = zeros(num_channels, length(peak_data(i).peaks.Alpha.P200.channel(1).max_amplitude));
    for ch = 1:num_channels
        alpha_peaks_w(ch, :) = peak_data(i).peaks.Alpha.P200.channel(ch).max_amplitude;
    end
    f3_z_w = mean(alpha_peaks_w(f3_idx, W_trials), 'omitnan');
    f4_z_w = mean(alpha_peaks_w(f4_idx, W_trials), 'omitnan');
    faa_zscores(1) = f4_z_w - f3_z_w;  % FAA Z-Score = F4 - F3
    % T1 segments
    for w = 2:length(all_windows)
        win = all_windows{w};
        alpha_peaks = zeros(num_channels, length(peak_data(i).peaks.Alpha.(win).channel(1).max_amplitude));
        for ch = 1:num_channels
            alpha_peaks(ch, :) = peak_data(i).peaks.Alpha.(win).channel(ch).max_amplitude;
        end
        f3_z = mean(alpha_peaks(f3_idx, T1_trials), 'omitnan');
        f4_z = mean(alpha_peaks(f4_idx, T1_trials), 'omitnan');
        faa_zscores(w) = f4_z - f3_z;
    end
    plot(timeCenters, faa_zscores, '-o', 'Color', [1, 0.5, 0], 'LineWidth', 2);
    % Shade ERP windows (skip baseline)
    for w = 2:length(all_windows)
        boundsIndex = w - 1;
        tStart = windowBounds(boundsIndex, 1);
        tEnd   = windowBounds(boundsIndex, 2);
        h = fill([tStart, tEnd, tEnd, tStart], [min(faa_zscores)-0.1, min(faa_zscores)-0.1, max(faa_zscores)+0.1, max(faa_zscores)+0.1], ...
            [0.8, 0.9, 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
    end
    xticks(timeCenters);
    xticklabels(all_windows);
    xlabel('Time (ms)');
    ylabel('FAA Z-Score');
    title(['FAA Z-Score Across ERP Windows - ' participant_id]);
    grid on;
    hold off;
    saveas(gcf, ['H2_ZScoreGraph_FAA_' participant_id '.png']);
    close(gcf);
    faa_zscores_all(i, :) = faa_zscores;  % Store for stats
    
    %% BAR Z-Score Graph (average across channels) with W baseline
    figure('Name', ['BAR Z-Scores - ' participant_id], 'Position', [100, 100, 800, 600]);
    hold on;
    bar_zscores = zeros(1, length(all_windows));
    % W segment
    alpha_peaks_w = zeros(num_channels, length(peak_data(i).peaks.Alpha.P200.channel(1).max_amplitude));
    beta_peaks_w  = zeros(num_channels, length(peak_data(i).peaks.Beta.P200.channel(1).max_amplitude));
    for ch = 1:num_channels
        alpha_peaks_w(ch, :) = peak_data(i).peaks.Alpha.P200.channel(ch).max_amplitude;
        beta_peaks_w(ch, :)  = peak_data(i).peaks.Beta.P200.channel(ch).max_amplitude;
    end
    bar_per_trial_w = beta_peaks_w(:, W_trials) ./ alpha_peaks_w(:, W_trials);
    bar_T1_w = mean(bar_per_trial_w, 2, 'omitnan');
    bar_zscores(1) = mean(bar_T1_w, 'omitnan');
    % T1 segments
    for w = 2:length(all_windows)
        win = all_windows{w};
        alpha_peaks = zeros(num_channels, length(peak_data(i).peaks.Alpha.(win).channel(1).max_amplitude));
        beta_peaks  = zeros(num_channels, length(peak_data(i).peaks.Beta.(win).channel(1).max_amplitude));
        for ch = 1:num_channels
            alpha_peaks(ch, :) = peak_data(i).peaks.Alpha.(win).channel(ch).max_amplitude;
            beta_peaks(ch, :)  = peak_data(i).peaks.Beta.(win).channel(ch).max_amplitude;
        end
        bar_per_trial = beta_peaks(:, T1_trials) ./ alpha_peaks(:, T1_trials);
        bar_T1 = mean(bar_per_trial, 2, 'omitnan');
        bar_zscores(w) = mean(bar_T1, 'omitnan');
    end
    plot(timeCenters, bar_zscores, '-o', 'Color', [1, 0.5, 0], 'LineWidth', 2);
    % Shade ERP windows (skip baseline)
    for w = 2:length(all_windows)
        boundsIndex = w - 1;
        tStart = windowBounds(boundsIndex, 1);
        tEnd   = windowBounds(boundsIndex, 2);
        h = fill([tStart, tEnd, tEnd, tStart], [min(bar_zscores)-0.1, min(bar_zscores)-0.1, max(bar_zscores)+0.1, max(bar_zscores)+0.1], ...
            [0.8, 0.9, 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
    end
    xticks(timeCenters);
    xticklabels(all_windows);
    xlabel('Time (ms)');
    ylabel('BAR Z-Score');
    title(['BAR Z-Score Across ERP Windows - ' participant_id]);
    grid on;
    hold off;
    saveas(gcf, ['H2_ZScoreGraph_BAR_' participant_id '.png']);
    close(gcf);
    bar_zscores_all(i, :) = bar_zscores;
    
    %% TBR Z-Score Graph (average across channels) with W baseline
    figure('Name', ['TBR Z-Scores - ' participant_id], 'Position', [100, 100, 800, 600]);
    hold on;
    tbr_zscores = zeros(1, length(all_windows));
    % W segment
    theta_peaks_w = zeros(num_channels, length(peak_data(i).peaks.Theta.P200.channel(1).max_amplitude));
    beta_peaks_w  = zeros(num_channels, length(peak_data(i).peaks.Beta.P200.channel(1).max_amplitude));
    for ch = 1:num_channels
        theta_peaks_w(ch, :) = peak_data(i).peaks.Theta.P200.channel(ch).max_amplitude;
        beta_peaks_w(ch, :)  = peak_data(i).peaks.Beta.P200.channel(ch).max_amplitude;
    end
    tbr_per_trial_w = theta_peaks_w(:, W_trials) ./ beta_peaks_w(:, W_trials);
    tbr_T1_w = mean(tbr_per_trial_w, 2, 'omitnan');
    tbr_zscores(1) = mean(tbr_T1_w, 'omitnan');
    % T1 segments
    for w = 2:length(all_windows)
        win = all_windows{w};
        theta_peaks = zeros(num_channels, length(peak_data(i).peaks.Theta.(win).channel(1).max_amplitude));
        beta_peaks  = zeros(num_channels, length(peak_data(i).peaks.Beta.(win).channel(1).max_amplitude));
        for ch = 1:num_channels
            theta_peaks(ch, :) = peak_data(i).peaks.Theta.(win).channel(ch).max_amplitude;
            beta_peaks(ch, :)  = peak_data(i).peaks.Beta.(win).channel(ch).max_amplitude;
        end
        tbr_per_trial = theta_peaks(:, T1_trials) ./ beta_peaks(:, T1_trials);
        tbr_T1 = mean(tbr_per_trial, 2, 'omitnan');
        tbr_zscores(w) = mean(tbr_T1, 'omitnan');
    end
    plot(timeCenters, tbr_zscores, '-o', 'Color', [1, 0.5, 0], 'LineWidth', 2);
    % Shade ERP windows (skip baseline)
    for w = 2:length(all_windows)
        boundsIndex = w - 1;
        tStart = windowBounds(boundsIndex, 1);
        tEnd   = windowBounds(boundsIndex, 2);
        h = fill([tStart, tEnd, tEnd, tStart], [min(tbr_zscores)-0.1, min(tbr_zscores)-0.1, max(tbr_zscores)+0.1, max(tbr_zscores)+0.1], ...
            [0.8, 0.9, 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
    end
    xticks(timeCenters);
    xticklabels(all_windows);
    xlabel('Time (ms)');
    ylabel('TBR Z-Score');
    title(['TBR Z-Score Across ERP Windows - ' participant_id]);
    grid on;
    hold off;
    saveas(gcf, ['H2_ZScoreGraph_TBR_' participant_id '.png']);
    close(gcf);
    tbr_zscores_all(i, :) = tbr_zscores;
    
    %% HB Z-Score Graph (average across channels) with W baseline
    figure('Name', ['HB Z-Scores - ' participant_id], 'Position', [100, 100, 800, 600]);
    hold on;
    hb_zscores = zeros(1, length(all_windows));
    % W segment
    hb_peaks_w = zeros(num_channels, length(peak_data(i).peaks.HighBeta.P200.channel(1).max_amplitude));
    for ch = 1:num_channels
        hb_peaks_w(ch, :) = peak_data(i).peaks.HighBeta.P200.channel(ch).max_amplitude;
    end
    hb_T1_w = mean(hb_peaks_w(:, W_trials), 2, 'omitnan');
    hb_zscores(1) = mean(hb_T1_w, 'omitnan');
    % T1 segments
    for w = 2:length(all_windows)
        win = all_windows{w};
        hb_peaks_all = zeros(num_channels, length(peak_data(i).peaks.HighBeta.(win).channel(1).max_amplitude));
        for ch = 1:num_channels
            hb_peaks_all(ch, :) = peak_data(i).peaks.HighBeta.(win).channel(ch).max_amplitude;
        end
        hb_T1 = mean(hb_peaks_all(:, T1_trials), 2, 'omitnan');
        hb_zscores(w) = mean(hb_T1, 'omitnan');
    end
    plot(timeCenters, hb_zscores, '-o', 'Color', [1, 0.5, 0], 'LineWidth', 2);
    % Shade ERP windows (skip baseline)
    for w = 2:length(all_windows)
        boundsIndex = w - 1;
        tStart = windowBounds(boundsIndex, 1);
        tEnd   = windowBounds(boundsIndex, 2);
        h = fill([tStart, tEnd, tEnd, tStart], [min(hb_zscores)-0.1, min(hb_zscores)-0.1, max(hb_zscores)+0.1, max(hb_zscores)+0.1], ...
            [0.8, 0.9, 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
    end
    xticks(timeCenters);
    xticklabels(all_windows);
    xlabel('Time (ms)');
    ylabel('HB Z-Score');
    title(['HB Z-Score Across ERP Windows - ' participant_id]);
    grid on;
    hold off;
    saveas(gcf, ['H2_ZScoreGraph_HB_' participant_id '.png']);
    close(gcf);
    hb_zscores_all(i, :) = hb_zscores;
end

%% --- Peak Graphs for Y* Group (PT2YY, PT3YY, PT4YY, PT7YN) ---
% FAA Peak Graph (Y* Group)
figure('Name', 'FAA Peaks - Y* Group', 'Position', [100, 100, 800, 600]);
hold on;
f3_peaks_y = zeros(length(y_group_indices), length(all_windows));
f4_peaks_y = zeros(length(y_group_indices), length(all_windows));
for p = 1:length(y_group_indices)
    i = y_group_indices(p);
    % W segment
    alpha_peaks_w = zeros(num_channels, length(peak_data(i).peaks.Alpha.P200.channel(1).max_amplitude));
    for ch = 1:num_channels
        alpha_peaks_w(ch, :) = peak_data(i).peaks.Alpha.P200.channel(ch).max_amplitude;
    end
    f3_peaks_y(p, 1) = mean(alpha_peaks_w(f3_idx, W_trials), 'omitnan');
    f4_peaks_y(p, 1) = mean(alpha_peaks_w(f4_idx, W_trials), 'omitnan');
    % T1 segments
    for w = 2:length(all_windows)
        win = all_windows{w};
        alpha_peaks = zeros(num_channels, length(peak_data(i).peaks.Alpha.(win).channel(1).max_amplitude));
        for ch = 1:num_channels
            alpha_peaks(ch, :) = peak_data(i).peaks.Alpha.(win).channel(ch).max_amplitude;
        end
        f3_peaks_y(p, w) = mean(alpha_peaks(f3_idx, T1_trials), 'omitnan');
        f4_peaks_y(p, w) = mean(alpha_peaks(f4_idx, T1_trials), 'omitnan');
    end
end
f3_peaks_y_mean = mean(f3_peaks_y, 1, 'omitnan');
f4_peaks_y_mean = mean(f4_peaks_y, 1, 'omitnan');
plot(timeCenters, f3_peaks_y_mean, '-o', 'Color', [1, 0.5, 0], 'LineWidth', 2, 'DisplayName', 'F3');
plot(timeCenters, f4_peaks_y_mean, '-o', 'Color', [0, 0.5, 1], 'LineWidth', 2, 'DisplayName', 'F4');
% Shade ERP windows (skip baseline)
for w = 2:length(all_windows)
    boundsIndex = w - 1;
    tStart = windowBounds(boundsIndex, 1);
    tEnd   = windowBounds(boundsIndex, 2);
    h = fill([tStart, tEnd, tEnd, tStart], [yMin_default, yMin_default, yMax_default, yMax_default], [0.8, 0.9, 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
end
xticks(timeCenters);
xticklabels(all_windows);
xlabel('Time (ms)');
ylabel('Peak Z-Score (μV)');
title('FAA Peaks Across ERP Components (F3 vs F4) - Y* Group (PT2YY, PT3YY, PT4YY, PT7YN)');
xlim([0, 800]);
ylim([yMin_default, yMax_default]);
grid on;
legend('Location','best');
hold off;
saveas(gcf, 'H2_PeakGraph_FAA_YStar.png');
close(gcf);

%% BAR Peak Graph (Y* Group)
figure('Name', 'BAR Peaks - Y* Group', 'Position', [100, 100, 800, 600]);
hold on;
bar_peaks_y = zeros(length(y_group_indices), num_channels, length(all_windows));
for p = 1:length(y_group_indices)
    i = y_group_indices(p);
    % W segment
    alpha_peaks_w = zeros(num_channels, length(peak_data(i).peaks.Alpha.P200.channel(1).max_amplitude));
    beta_peaks_w  = zeros(num_channels, length(peak_data(i).peaks.Beta.P200.channel(1).max_amplitude));
    for ch = 1:num_channels
        alpha_peaks_w(ch, :) = peak_data(i).peaks.Alpha.P200.channel(ch).max_amplitude;
        beta_peaks_w(ch, :)  = peak_data(i).peaks.Beta.P200.channel(ch).max_amplitude;
    end
    bar_per_trial_w = beta_peaks_w(:, W_trials) ./ alpha_peaks_w(:, W_trials);
    bar_peaks_y(p, :, 1) = mean(bar_per_trial_w, 2, 'omitnan');
    % T1 segments
    for w = 2:length(all_windows)
        win = all_windows{w};
        alpha_peaks = zeros(num_channels, length(peak_data(i).peaks.Alpha.(win).channel(1).max_amplitude));
        beta_peaks  = zeros(num_channels, length(peak_data(i).peaks.Beta.(win).channel(1).max_amplitude));
        for ch = 1:num_channels
            alpha_peaks(ch, :) = peak_data(i).peaks.Alpha.(win).channel(ch).max_amplitude;
            beta_peaks(ch, :)  = peak_data(i).peaks.Beta.(win).channel(ch).max_amplitude;
        end
        bar_per_trial = beta_peaks(:, T1_trials) ./ alpha_peaks(:, T1_trials);
        bar_peaks_y(p, :, w) = mean(bar_per_trial, 2, 'omitnan');
    end
end
bar_peaks_y_mean = mean(bar_peaks_y, 1, 'omitnan');
bar_peaks_y_mean = squeeze(bar_peaks_y_mean);
for ch = 1:num_channels
    plot(timeCenters, bar_peaks_y_mean(ch, :), '-o', 'Color', channel_colors(ch, :), ...
        'LineWidth', 2, 'DisplayName', chanlocs(ch).labels);
end
% Shade ERP windows (skip baseline)
for w = 2:length(all_windows)
    boundsIndex = w - 1;
    tStart = windowBounds(boundsIndex, 1);
    tEnd   = windowBounds(boundsIndex, 2);
    h = fill([tStart, tEnd, tEnd, tStart], [min(bar_peaks_y_mean(:))-0.1, min(bar_peaks_y_mean(:))-0.1, max(bar_peaks_y_mean(:))+0.1, max(bar_peaks_y_mean(:))+0.1], ...
        [0.8, 0.9, 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
end
xticks(timeCenters);
xticklabels(all_windows);
xlabel('Time (ms)');
ylabel('Peak BAR');
title('BAR Peaks Across ERP Windows - Y* Group (PT2YY, PT3YY, PT4YY, PT7YN)');
grid on;
legend('Location','best');
hold off;
saveas(gcf, 'H2_PeakGraph_BAR_YStar.png');
close(gcf);

%% TBR Peak Graph (Y* Group)
figure('Name', 'TBR Peaks - Y* Group', 'Position', [100, 100, 800, 600]);
hold on;
tbr_peaks_y = zeros(length(y_group_indices), num_channels, length(all_windows));
for p = 1:length(y_group_indices)
    i = y_group_indices(p);
    % W segment
    theta_peaks_w = zeros(num_channels, length(peak_data(i).peaks.Theta.P200.channel(1).max_amplitude));
    beta_peaks_w  = zeros(num_channels, length(peak_data(i).peaks.Beta.P200.channel(1).max_amplitude));
    for ch = 1:num_channels
        theta_peaks_w(ch, :) = peak_data(i).peaks.Theta.P200.channel(ch).max_amplitude;
        beta_peaks_w(ch, :)  = peak_data(i).peaks.Beta.P200.channel(ch).max_amplitude;
    end
    tbr_per_trial_w = theta_peaks_w(:, W_trials) ./ beta_peaks_w(:, W_trials);
    tbr_peaks_y(p, :, 1) = mean(tbr_per_trial_w, 2, 'omitnan');
    % T1 segments
    for w = 2:length(all_windows)
        win = all_windows{w};
        theta_peaks = zeros(num_channels, length(peak_data(i).peaks.Theta.(win).channel(1).max_amplitude));
        beta_peaks  = zeros(num_channels, length(peak_data(i).peaks.Beta.(win).channel(1).max_amplitude));
        for ch = 1:num_channels
            theta_peaks(ch, :) = peak_data(i).peaks.Theta.(win).channel(ch).max_amplitude;
            beta_peaks(ch, :)  = peak_data(i).peaks.Beta.(win).channel(ch).max_amplitude;
        end
        tbr_per_trial = theta_peaks(:, T1_trials) ./ beta_peaks(:, T1_trials);
        tbr_peaks_y(p, :, w) = mean(tbr_per_trial, 2, 'omitnan');
    end
end
tbr_peaks_y_mean = mean(tbr_peaks_y, 1, 'omitnan');
tbr_peaks_y_mean = squeeze(tbr_peaks_y_mean);
for ch = 1:num_channels
    plot(timeCenters, tbr_peaks_y_mean(ch, :), '-o', 'Color', channel_colors(ch, :), ...
        'LineWidth', 2, 'DisplayName', chanlocs(ch).labels);
end
% Shade ERP windows (skip baseline)
for w = 2:length(all_windows)
    boundsIndex = w - 1;
    tStart = windowBounds(boundsIndex, 1);
    tEnd   = windowBounds(boundsIndex, 2);
    h = fill([tStart, tEnd, tEnd, tStart], [min(tbr_peaks_y_mean(:))-0.1, min(tbr_peaks_y_mean(:))-0.1, max(tbr_peaks_y_mean(:))+0.1, max(tbr_peaks_y_mean(:))+0.1], ...
        [0.8, 0.9, 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
end
xticks(timeCenters);
xticklabels(all_windows);
xlabel('Time (ms)');
ylabel('Peak TBR');
title('TBR Peaks Across ERP Windows - Y* Group (PT2YY, PT3YY, PT4YY, PT7YN)');
grid on;
legend('Location','best');
hold off;
saveas(gcf, 'H2_PeakGraph_TBR_YStar.png');
close(gcf);

%% HB Peak Graph (Y* Group)
figure('Name', 'HB Peaks - Y* Group', 'Position', [100, 100, 800, 600]);
hold on;
hb_peaks_y = zeros(length(y_group_indices), num_channels, length(all_windows));
for p = 1:length(y_group_indices)
    i = y_group_indices(p);
    % W segment
    hb_peaks_w = zeros(num_channels, length(peak_data(i).peaks.HighBeta.P200.channel(1).max_amplitude));
    for ch = 1:num_channels
        hb_peaks_w(ch, :) = peak_data(i).peaks.HighBeta.P200.channel(ch).max_amplitude;
    end
    hb_peaks_y(p, :, 1) = mean(hb_peaks_w(:, W_trials), 2, 'omitnan');
    % T1 segments
    for w = 2:length(all_windows)
        win = all_windows{w};
        hb_peaks_all = zeros(num_channels, length(peak_data(i).peaks.HighBeta.(win).channel(1).max_amplitude));
        for ch = 1:num_channels
            hb_peaks_all(ch, :) = peak_data(i).peaks.HighBeta.(win).channel(ch).max_amplitude;
        end
        hb_peaks_y(p, :, w) = mean(hb_peaks_all(:, T1_trials), 2, 'omitnan');
    end
end
hb_peaks_y_mean = mean(hb_peaks_y, 1, 'omitnan');
hb_peaks_y_mean = squeeze(hb_peaks_y_mean);
for ch = 1:num_channels
    plot(timeCenters, hb_peaks_y_mean(ch, :), '-o', 'Color', channel_colors(ch, :), ...
        'LineWidth', 2, 'DisplayName', chanlocs(ch).labels);
end
% Shade ERP windows (skip baseline)
for w = 2:length(all_windows)
    boundsIndex = w - 1;
    tStart = windowBounds(boundsIndex, 1);
    tEnd   = windowBounds(boundsIndex, 2);
    h = fill([tStart, tEnd, tEnd, tStart], [min(hb_peaks_y_mean(:))-0.1, min(hb_peaks_y_mean(:))-0.1, max(hb_peaks_y_mean(:))+0.1, max(hb_peaks_y_mean(:))+0.1], ...
        [0.8, 0.9, 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
end
xticks(timeCenters);
xticklabels(all_windows);
xlabel('Time (ms)');
ylabel('Peak Z-Score (μV)');
title('HB Peaks Across ERP Windows - Y* Group (PT2YY, PT3YY, PT4YY, PT7YN)');
grid on;
legend('Location','best');
hold off;
saveas(gcf, 'H2_PeakGraph_HB_YStar.png');
close(gcf);

%% --- Z-Score Graphs for Y* Group ---
% FAA Z-Score Graph (Y* Group)
figure('Name', 'FAA Z-Scores - Y* Group', 'Position', [100, 100, 800, 600]);
hold on;
faa_zscores_y = zeros(length(y_group_indices), length(all_windows));
for p = 1:length(y_group_indices)
    i = y_group_indices(p);
    % W segment
    alpha_peaks_w = zeros(num_channels, length(peak_data(i).peaks.Alpha.P200.channel(1).max_amplitude));
    for ch = 1:num_channels
        alpha_peaks_w(ch, :) = peak_data(i).peaks.Alpha.P200.channel(ch).max_amplitude;
    end
    f3_z_w = mean(alpha_peaks_w(f3_idx, W_trials), 'omitnan');
    f4_z_w = mean(alpha_peaks_w(f4_idx, W_trials), 'omitnan');
    faa_zscores_y(p, 1) = f4_z_w - f3_z_w;
    % T1 segments
    for w = 2:length(all_windows)
        win = all_windows{w};
        alpha_peaks = zeros(num_channels, length(peak_data(i).peaks.Alpha.(win).channel(1).max_amplitude));
        for ch = 1:num_channels
            alpha_peaks(ch, :) = peak_data(i).peaks.Alpha.(win).channel(ch).max_amplitude;
        end
        f3_z = mean(alpha_peaks(f3_idx, T1_trials), 'omitnan');
        f4_z = mean(alpha_peaks(f4_idx, T1_trials), 'omitnan');
        faa_zscores_y(p, w) = f4_z - f3_z;
    end
end
faa_zscores_y_mean = mean(faa_zscores_y, 1, 'omitnan');
plot(timeCenters, faa_zscores_y_mean, '-o', 'Color', [1, 0.5, 0], 'LineWidth', 2);
% Shade ERP windows (skip baseline)
for w = 2:length(all_windows)
    boundsIndex = w - 1;
    tStart = windowBounds(boundsIndex, 1);
    tEnd   = windowBounds(boundsIndex, 2);
    h = fill([tStart, tEnd, tEnd, tStart], [min(faa_zscores_y_mean)-0.1, min(faa_zscores_y_mean)-0.1, max(faa_zscores_y_mean)+0.1, max(faa_zscores_y_mean)+0.1], ...
        [0.8, 0.9, 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
end
xticks(timeCenters);
xticklabels(all_windows);
xlabel('Time (ms)');
ylabel('FAA Z-Score');
title('FAA Z-Score Across ERP Windows - Y* Group (PT2YY, PT3YY, PT4YY, PT7YN)');
grid on;
hold off;
saveas(gcf, 'H2_ZScoreGraph_FAA_YStar.png');
close(gcf);

%% BAR Z-Score Graph (Y* Group)
figure('Name', 'BAR Z-Scores - Y* Group', 'Position', [100, 100, 800, 600]);
hold on;
bar_zscores_y = zeros(length(y_group_indices), length(all_windows));
for p = 1:length(y_group_indices)
    i = y_group_indices(p);
    % W segment
    alpha_peaks_w = zeros(num_channels, length(peak_data(i).peaks.Alpha.P200.channel(1).max_amplitude));
    beta_peaks_w  = zeros(num_channels, length(peak_data(i).peaks.Beta.P200.channel(1).max_amplitude));
    for ch = 1:num_channels
        alpha_peaks_w(ch, :) = peak_data(i).peaks.Alpha.P200.channel(ch).max_amplitude;
        beta_peaks_w(ch, :)  = peak_data(i).peaks.Beta.P200.channel(ch).max_amplitude;
    end
    bar_per_trial_w = beta_peaks_w(:, W_trials) ./ alpha_peaks_w(:, W_trials);
    bar_T1_w = mean(bar_per_trial_w, 2, 'omitnan');
    bar_zscores_y(p, 1) = mean(bar_T1_w, 'omitnan');
    % T1 segments
    for w = 2:length(all_windows)
        win = all_windows{w};
        alpha_peaks = zeros(num_channels, length(peak_data(i).peaks.Alpha.(win).channel(1).max_amplitude));
        beta_peaks  = zeros(num_channels, length(peak_data(i).peaks.Beta.(win).channel(1).max_amplitude));
        for ch = 1:num_channels
            alpha_peaks(ch, :) = peak_data(i).peaks.Alpha.(win).channel(ch).max_amplitude;
            beta_peaks(ch, :)  = peak_data(i).peaks.Beta.(win).channel(ch).max_amplitude;
        end
        bar_per_trial = beta_peaks(:, T1_trials) ./ alpha_peaks(:, T1_trials);
        bar_T1 = mean(bar_per_trial, 2, 'omitnan');
        bar_zscores_y(p, w) = mean(bar_T1, 'omitnan');
    end
end
bar_zscores_y_mean = mean(bar_zscores_y, 1, 'omitnan');
plot(timeCenters, bar_zscores_y_mean, '-o', 'Color', [1, 0.5, 0], 'LineWidth', 2);
% Shade ERP windows (skip baseline)
for w = 2:length(all_windows)
    boundsIndex = w - 1;
    tStart = windowBounds(boundsIndex, 1);
    tEnd   = windowBounds(boundsIndex, 2);
    h = fill([tStart, tEnd, tEnd, tStart], [min(bar_zscores_y_mean)-0.1, min(bar_zscores_y_mean)-0.1, max(bar_zscores_y_mean)+0.1, max(bar_zscores_y_mean)+0.1], ...
        [0.8, 0.9, 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
end
xticks(timeCenters);
xticklabels(all_windows);
xlabel('Time (ms)');
ylabel('BAR Z-Score');
title('BAR Z-Score Across ERP Windows - Y* Group (PT2YY, PT3YY, PT4YY, PT7YN)');
grid on;
hold off;
saveas(gcf, 'H2_ZScoreGraph_BAR_YStar.png');
close(gcf);

%% TBR Z-Score Graph (Y* Group)
figure('Name', 'TBR Z-Scores - Y* Group', 'Position', [100, 100, 800, 600]);
hold on;
tbr_zscores_y = zeros(length(y_group_indices), length(all_windows));
for p = 1:length(y_group_indices)
    i = y_group_indices(p);
    % W segment
    theta_peaks_w = zeros(num_channels, length(peak_data(i).peaks.Theta.P200.channel(1).max_amplitude));
    beta_peaks_w  = zeros(num_channels, length(peak_data(i).peaks.Beta.P200.channel(1).max_amplitude));
    for ch = 1:num_channels
        theta_peaks_w(ch, :) = peak_data(i).peaks.Theta.P200.channel(ch).max_amplitude;
        beta_peaks_w(ch, :)  = peak_data(i).peaks.Beta.P200.channel(ch).max_amplitude;
    end
    tbr_per_trial_w = theta_peaks_w(:, W_trials) ./ beta_peaks_w(:, W_trials);
    temp = mean(tbr_per_trial_w, 2, 'omitnan');
    tbr_zscores_y(p, 1) = mean(temp, 'omitnan');

    %tbr_zscores_y(p, 1) = mean(tbr_per_trial_w, 'omitnan');
    % T1 segments
    for w = 2:length(all_windows)
        win = all_windows{w};
        theta_peaks = zeros(num_channels, length(peak_data(i).peaks.Theta.(win).channel(1).max_amplitude));
        beta_peaks  = zeros(num_channels, length(peak_data(i).peaks.Beta.(win).channel(1).max_amplitude));
        for ch = 1:num_channels
            theta_peaks(ch, :) = peak_data(i).peaks.Theta.(win).channel(ch).max_amplitude;
            beta_peaks(ch, :)  = peak_data(i).peaks.Beta.(win).channel(ch).max_amplitude;
        end
        tbr_per_trial = theta_peaks(:, T1_trials) ./ beta_peaks(:, T1_trials);
        temp = mean(tbr_per_trial, 2, 'omitnan');
        tbr_zscores_y(p, w) = mean(temp, 'omitnan');
    end
end
tbr_zscores_y_mean = mean(tbr_zscores_y, 1, 'omitnan');
plot(timeCenters, tbr_zscores_y_mean, '-o', 'Color', [1, 0.5, 0], 'LineWidth', 2);
% Shade ERP windows (skip baseline)
for w = 2:length(all_windows)
    boundsIndex = w - 1;
    tStart = windowBounds(boundsIndex, 1);
    tEnd   = windowBounds(boundsIndex, 2);
    h = fill([tStart, tEnd, tEnd, tStart], [min(tbr_zscores_y_mean)-0.1, min(tbr_zscores_y_mean)-0.1, max(tbr_zscores_y_mean)+0.1, max(tbr_zscores_y_mean)+0.1], ...
        [0.8, 0.9, 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
end
xticks(timeCenters);
xticklabels(all_windows);
xlabel('Time (ms)');
ylabel('TBR Z-Score');
title('TBR Z-Score Across ERP Windows - Y* Group (PT2YY, PT3YY, PT4YY, PT7YN)');
grid on;
hold off;
saveas(gcf, 'H2_ZScoreGraph_TBR_YStar.png');
close(gcf);

%% HB Z-Score Graph (Y* Group)
figure('Name', 'HB Z-Scores - Y* Group', 'Position', [100, 100, 800, 600]);
hold on;
hb_zscores_y = zeros(length(y_group_indices), length(all_windows));
for p = 1:length(y_group_indices)
    i = y_group_indices(p);
    % W segment
    hb_peaks_w = zeros(num_channels, length(peak_data(i).peaks.HighBeta.P200.channel(1).max_amplitude));
    for ch = 1:num_channels
        hb_peaks_w(ch, :) = peak_data(i).peaks.HighBeta.P200.channel(ch).max_amplitude;
    end
    temp = mean(hb_peaks_w(:, W_trials), 'omitnan');
    hb_zscores_y(p, 1) =  mean(temp, 'omitnan');
    % T1 segments
    for w = 2:length(all_windows)
        win = all_windows{w};
        hb_peaks_all = zeros(num_channels, length(peak_data(i).peaks.HighBeta.(win).channel(1).max_amplitude));
        for ch = 1:num_channels
            hb_peaks_all(ch, :) = peak_data(i).peaks.HighBeta.(win).channel(ch).max_amplitude;
        end
        temp = mean(hb_peaks_all(:, T1_trials), 'omitnan')
        hb_zscores_y(p, w) = mean(temp, 'omitnan') ;
    end
end
hb_zscores_y_mean = mean(hb_zscores_y, 1, 'omitnan');
plot(timeCenters, hb_zscores_y_mean, '-o', 'Color', [1, 0.5, 0], 'LineWidth', 2);
% Shade ERP windows (skip baseline)
for w = 2:length(all_windows)
    boundsIndex = w - 1;
    tStart = windowBounds(boundsIndex, 1);
    tEnd   = windowBounds(boundsIndex, 2);
    h = fill([tStart, tEnd, tEnd, tStart], [min(hb_zscores_y_mean)-0.1, min(hb_zscores_y_mean)-0.1, max(hb_zscores_y_mean)+0.1, max(hb_zscores_y_mean)+0.1], ...
        [0.8, 0.9, 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
end
xticks(timeCenters);
xticklabels(all_windows);
xlabel('Time (ms)');
ylabel('HB Z-Score');
title('HB Z-Score Across ERP Windows - Y* Group (PT2YY, PT3YY, PT4YY, PT7YN)');
grid on;
hold off;
saveas(gcf, 'H2_ZScoreGraph_HB_YStar.png');
close(gcf);

%% --- Descriptive Statistics and Spearman’s Correlation ---
% Compute mean z-scores for T1 (average across ERP windows P200, N200, P300, LPP) and W
t1_windows = 2:length(all_windows);  % P200, N200, P300, LPP
faa_zscores_t1 = mean(faa_zscores_all(:, t1_windows), 2, 'omitnan');  % Mean across T1 windows
faa_zscores_w  = faa_zscores_all(:, 1);  % W segment
bar_zscores_t1 = mean(bar_zscores_all(:, t1_windows), 2, 'omitnan');
bar_zscores_w  = bar_zscores_all(:, 1);
tbr_zscores_t1 = mean(tbr_zscores_all(:, t1_windows), 2, 'omitnan');
tbr_zscores_w  = tbr_zscores_all(:, 1);
hb_zscores_t1  = mean(hb_zscores_all(:, t1_windows), 2, 'omitnan');
hb_zscores_w   = hb_zscores_all(:, 1);

% Descriptive Statistics
desc_stats = table();
desc_stats.Metric = {'FAA'; 'FAA'; 'BAR'; 'BAR'; 'TBR'; 'TBR'; 'HB'; 'HB'};
desc_stats.Condition = {'T1'; 'W'; 'T1'; 'W'; 'T1'; 'W'; 'T1'; 'W'};
desc_stats.Mean = [
    mean(faa_zscores_t1, 'omitnan');
    mean(faa_zscores_w, 'omitnan');
    mean(bar_zscores_t1, 'omitnan');
    mean(bar_zscores_w, 'omitnan');
    mean(tbr_zscores_t1, 'omitnan');
    mean(tbr_zscores_w, 'omitnan');
    mean(hb_zscores_t1, 'omitnan');
    mean(hb_zscores_w, 'omitnan')
];
desc_stats.SD = [
    std(faa_zscores_t1, 'omitnan');
    std(faa_zscores_w, 'omitnan');
    std(bar_zscores_t1, 'omitnan');
    std(bar_zscores_w, 'omitnan');
    std(tbr_zscores_t1, 'omitnan');
    std(tbr_zscores_w, 'omitnan');
    std(hb_zscores_t1, 'omitnan');
    std(hb_zscores_w, 'omitnan')
];
disp('Descriptive Statistics (All Participants):');
disp(desc_stats);
writetable(desc_stats, 'H2_Descriptive_Statistics.csv');

% Spearman’s Correlation (T1 vs W)
[rho_faa, p_faa] = corr(faa_zscores_t1, faa_zscores_w, 'Type', 'Spearman', 'Rows', 'complete');
[rho_bar, p_bar] = corr(bar_zscores_t1, bar_zscores_w, 'Type', 'Spearman', 'Rows', 'complete');
[rho_tbr, p_tbr] = corr(tbr_zscores_t1, tbr_zscores_w, 'Type', 'Spearman', 'Rows', 'complete');
[rho_hb, p_hb]   = corr(hb_zscores_t1, hb_zscores_w, 'Type', 'Spearman', 'Rows', 'complete');

spearman_stats = table();
spearman_stats.Metric = {'FAA'; 'BAR'; 'TBR'; 'HB'};
spearman_stats.Rho = [rho_faa; rho_bar; rho_tbr; rho_hb];
spearman_stats.PValue = [p_faa; p_bar; p_tbr; p_hb];
disp('Spearman’s Correlation (T1 vs W, All Participants):');
disp(spearman_stats);
writetable(spearman_stats, 'H2_Spearman_Correlation.csv');

%% --- Cohen’s d for Y* Group (and PT1NN for discussion) ---
% Y* Group (PT2YY, PT3YY, PT4YY, PT7YN)
faa_zscores_y_t1 = mean(faa_zscores_all(y_group_indices, t1_windows), 2, 'omitnan');
faa_zscores_y_w  = faa_zscores_all(y_group_indices, 1);
bar_zscores_y_t1 = mean(bar_zscores_all(y_group_indices, t1_windows), 2, 'omitnan');
bar_zscores_y_w  = bar_zscores_all(y_group_indices, 1);
tbr_zscores_y_t1 = mean(tbr_zscores_all(y_group_indices, t1_windows), 2, 'omitnan');
tbr_zscores_y_w  = tbr_zscores_all(y_group_indices, 1);
hb_zscores_y_t1  = mean(hb_zscores_all(y_group_indices, t1_windows), 2, 'omitnan');
hb_zscores_y_w   = hb_zscores_all(y_group_indices, 1);

% PT1NN (for discussion, N=1)
faa_zscores_pt1_t1 = mean(faa_zscores_all(1, t1_windows), 'omitnan');
faa_zscores_pt1_w  = faa_zscores_all(1, 1);
bar_zscores_pt1_t1 = mean(bar_zscores_all(1, t1_windows), 'omitnan');
bar_zscores_pt1_w  = bar_zscores_all(1, 1);
tbr_zscores_pt1_t1 = mean(tbr_zscores_all(1, t1_windows), 'omitnan');
tbr_zscores_pt1_w  = tbr_zscores_all(1, 1);
hb_zscores_pt1_t1  = mean(hb_zscores_all(1, t1_windows), 'omitnan');
hb_zscores_pt1_w   = hb_zscores_all(1, 1);

% Cohen’s d for Y* Group
cohens_d_y = zeros(4, 1);
% FAA
mean_diff_faa_y = mean(faa_zscores_y_t1) - mean(faa_zscores_y_w);
pooled_sd_faa_y = sqrt((std(faa_zscores_y_t1)^2 + std(faa_zscores_y_w)^2) / 2);
cohens_d_y(1) = mean_diff_faa_y / pooled_sd_faa_y;
% BAR
mean_diff_bar_y = mean(bar_zscores_y_t1) - mean(bar_zscores_y_w);
pooled_sd_bar_y = sqrt((std(bar_zscores_y_t1)^2 + std(bar_zscores_y_w)^2) / 2);
cohens_d_y(2) = mean_diff_bar_y / pooled_sd_bar_y;
% TBR
mean_diff_tbr_y = mean(tbr_zscores_y_t1) - mean(tbr_zscores_y_w);
pooled_sd_tbr_y = sqrt((std(tbr_zscores_y_t1)^2 + std(tbr_zscores_y_w)^2) / 2);
cohens_d_y(3) = mean_diff_tbr_y / pooled_sd_tbr_y;
% HB
mean_diff_hb_y = mean(hb_zscores_y_t1) - mean(hb_zscores_y_w);
pooled_sd_hb_y = sqrt((std(hb_zscores_y_t1)^2 + std(hb_zscores_y_w)^2) / 2);
cohens_d_y(4) = mean_diff_hb_y / pooled_sd_hb_y;

% Cohen’s d for PT1NN (N=1, for discussion only)
cohens_d_pt1 = zeros(4, 1);
% FAA
mean_diff_faa_pt1 = faa_zscores_pt1_t1 - faa_zscores_pt1_w;
pooled_sd_faa_pt1 = pooled_sd_faa_y;  % Using Y* group SD as proxy
cohens_d_pt1(1) = mean_diff_faa_pt1 / pooled_sd_faa_pt1;
% BAR
mean_diff_bar_pt1 = bar_zscores_pt1_t1 - bar_zscores_pt1_w;
pooled_sd_bar_pt1 = pooled_sd_bar_y;
cohens_d_pt1(2) = mean_diff_bar_pt1 / pooled_sd_bar_pt1;
% TBR
mean_diff_tbr_pt1 = tbr_zscores_pt1_t1 - tbr_zscores_pt1_w;
pooled_sd_tbr_pt1 = pooled_sd_tbr_y;
cohens_d_pt1(3) = mean_diff_tbr_pt1 / pooled_sd_tbr_pt1;
% HB
mean_diff_hb_pt1 = hb_zscores_pt1_t1 - hb_zscores_pt1_w;
pooled_sd_hb_pt1 = pooled_sd_hb_y;
cohens_d_pt1(4) = mean_diff_hb_pt1 / pooled_sd_hb_pt1;

% Cohen’s d Table
cohens_d_table = table();
cohens_d_table.Metric = {'FAA'; 'BAR'; 'TBR'; 'HB'; 'FAA_PT1NN'; 'BAR_PT1NN'; 'TBR_PT1NN'; 'HB_PT1NN'};
cohens_d_table.Group  = [cohens_d_y; cohens_d_pt1];
disp('Cohens (T1 vs W, All Participants):');
disp(cohens_d_table);
writetable(cohens_d_table, 'H2_Cohens_Correlation.csv');


%% H3 - Peak and Z-Score Graphs for All Participants and Y* Group with Baseline, Statistics, and Cohen's d
% Define parameters
T2_trials = segments.T2;      % Trial indices for T2 condition (271:300)
W_trials = segments.W;        % Trial indices for W condition (baseline)
chanlocs = all_eeg(1).EEG.chanlocs;  % Channel locations
num_channels = length(chanlocs);     % Number of EEG channels (assumed 8: AF3, AF4, F3, F4, FC5, FC6, P8, F8)
num_participants = length(participants);  % Number of participants (5)
all_windows = {'W', 'P200', 'N200', 'P300', 'LPP'};  % ERP windows including W as baseline
metrics = {'FAA', 'BAR', 'TBR', 'HB'};
y_group_indices = [2, 3, 4, 5];  % Indices for Y* group: PT2YY, PT3YY, PT4YY, PT7YN

% Find F3 and F4 channel indices (for FAA)
f3_idx = find(strcmp({chanlocs.labels}, 'F3'));
f4_idx = find(strcmp({chanlocs.labels}, 'F4'));
if isempty(f3_idx) || isempty(f4_idx)
    error('F3 or F4 channel not found in chanlocs.');
end

% Define colors for channels (for peak graphs with multiple channels)
channel_colors = lines(num_channels);  % MATLAB's 'lines' colormap for distinct colors

% Define time mapping for ERP windows (ms)
% Baseline (W): assume 0–100 ms with midpoint = 50 ms
% P200: 150–250 ms, midpoint = 200 ms
% N200: 180–240 ms, midpoint = 210 ms
% P300: 300–400 ms, midpoint = 350 ms
% LPP: 400–800 ms, midpoint = 600 ms
timeCenters = [50, 175, 250, 350, 510];  % one value per window (W, P200, N200, P300, LPP)
% For shading, we use the ERP windows (excluding baseline, index 1)
windowBounds = [150 200;   % P200 bounds
                210 290;   % N200 bounds
                310 390;   % P300 bounds
                420 600];  % LPP bounds
yMin_default = -1;
yMax_default = 2;

% Initialize arrays to store z-scores for statistical analysis
faa_zscores_all = zeros(num_participants, length(all_windows));  % [participants, windows]
bar_zscores_all = zeros(num_participants, length(all_windows));
tbr_zscores_all = zeros(num_participants, length(all_windows));
hb_zscores_all = zeros(num_participants, length(all_windows));

%% --- Peak Graphs for Individuals ---
for i = 1:num_participants
    participant_id = participants{i};
    
    %% FAA Peak Graph (F3 vs F4) with W baseline
    figure('Name', ['FAA Peaks - ' participant_id], 'Position', [100, 100, 800, 600]);
    hold on;
    f3_peaks = zeros(1, length(all_windows));
    f4_peaks = zeros(1, length(all_windows));
    % W segment (baseline)
    alpha_peaks_w = zeros(num_channels, length(peak_data(i).peaks.Alpha.P200.channel(1).max_amplitude));
    for ch = 1:num_channels
        alpha_peaks_w(ch, :) = peak_data(i).peaks.Alpha.P200.channel(ch).max_amplitude;
    end
    f3_peaks(1) = mean(alpha_peaks_w(f3_idx, W_trials), 'omitnan');
    f4_peaks(1) = mean(alpha_peaks_w(f4_idx, W_trials), 'omitnan');
    % T2 segments
    for w = 2:length(all_windows)
        win = all_windows{w};
        alpha_peaks = zeros(num_channels, length(peak_data(i).peaks.Alpha.(win).channel(1).max_amplitude));
        for ch = 1:num_channels
            alpha_peaks(ch, :) = peak_data(i).peaks.Alpha.(win).channel(ch).max_amplitude;
        end
        f3_peaks(w) = mean(alpha_peaks(f3_idx, T2_trials), 'omitnan');
        f4_peaks(w) = mean(alpha_peaks(f4_idx, T2_trials), 'omitnan');
    end
    plot(timeCenters, f3_peaks, '-o', 'Color', [1, 0.5, 0], 'LineWidth', 2, 'DisplayName', 'F3');
    plot(timeCenters, f4_peaks, '-o', 'Color', [0, 0.5, 1], 'LineWidth', 2, 'DisplayName', 'F4');
    % Shade ERP windows only (skip baseline: index 1)
    for w = 2:length(all_windows)
        boundsIndex = w - 1;
        tStart = windowBounds(boundsIndex, 1);
        tEnd   = windowBounds(boundsIndex, 2);
        h = fill([tStart, tEnd, tEnd, tStart], [yMin_default, yMin_default, yMax_default, yMax_default], [0.8, 0.9, 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
    end
    xticks(timeCenters);
    xticklabels(all_windows);
    xlabel('Time (ms)');
    ylabel('Peak Z-Score (μV)');
    title(['FAA Peaks Across ERP Components (F3 vs F4) - ' participant_id]);
    xlim([0, 800]);
    ylim([yMin_default, yMax_default]);
    grid on;
    legend('Location','best');
    hold off;
    saveas(gcf, ['H3_PeakGraph_FAA_' participant_id '.png']);
    close(gcf);
    
    %% BAR Peak Graph (per channel) with W baseline
    figure('Name', ['BAR Peaks - ' participant_id], 'Position', [100, 100, 800, 600]);
    hold on;
    bar_peaks = zeros(num_channels, length(all_windows));
    % W segment
    alpha_peaks_w = zeros(num_channels, length(peak_data(i).peaks.Alpha.P200.channel(1).max_amplitude));
    beta_peaks_w = zeros(num_channels, length(peak_data(i).peaks.Beta.P200.channel(1).max_amplitude));
    for ch = 1:num_channels
        alpha_peaks_w(ch, :) = peak_data(i).peaks.Alpha.P200.channel(ch).max_amplitude;
        beta_peaks_w(ch, :) = peak_data(i).peaks.Beta.P200.channel(ch).max_amplitude;
    end
    bar_per_trial_w = beta_peaks_w(:, W_trials) ./ alpha_peaks_w(:, W_trials);
    bar_peaks(:, 1) = mean(bar_per_trial_w, 2, 'omitnan');
    % T2 segments
    for w = 2:length(all_windows)
        win = all_windows{w};
        alpha_peaks = zeros(num_channels, length(peak_data(i).peaks.Alpha.(win).channel(1).max_amplitude));
        beta_peaks  = zeros(num_channels, length(peak_data(i).peaks.Beta.(win).channel(1).max_amplitude));
        for ch = 1:num_channels
            alpha_peaks(ch, :) = peak_data(i).peaks.Alpha.(win).channel(ch).max_amplitude;
            beta_peaks(ch, :)  = peak_data(i).peaks.Beta.(win).channel(ch).max_amplitude;
        end
        bar_per_trial = beta_peaks(:, T2_trials) ./ alpha_peaks(:, T2_trials);
        bar_peaks(:, w) = mean(bar_per_trial, 2, 'omitnan');
    end
    for ch = 1:num_channels
        plot(timeCenters, bar_peaks(ch, :), '-o', 'Color', channel_colors(ch, :), ...
            'LineWidth', 2, 'DisplayName', chanlocs(ch).labels);
    end
    % Shade ERP windows (skip baseline)
    for w = 2:length(all_windows)
        boundsIndex = w - 1;
        tStart = windowBounds(boundsIndex, 1);
        tEnd   = windowBounds(boundsIndex, 2);
        h = fill([tStart, tEnd, tEnd, tStart], [min(bar_peaks(:))-0.1, min(bar_peaks(:))-0.1, max(bar_peaks(:))+0.1, max(bar_peaks(:))+0.1], ...
            [0.8, 0.9, 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
    end
    xticks(timeCenters);
    xticklabels(all_windows);
    xlabel('Time (ms)');
    ylabel('Peak BAR');
    title(['BAR Peaks Across ERP Windows - ' participant_id]);
    grid on;
    legend('Location','best');
    hold off;
    saveas(gcf, ['H3_PeakGraph_BAR_' participant_id '.png']);
    close(gcf);
    
    %% TBR Peak Graph (per channel) with W baseline
    figure('Name', ['TBR Peaks - ' participant_id], 'Position', [100, 100, 800, 600]);
    hold on;
    tbr_peaks = zeros(num_channels, length(all_windows));
    % W segment
    theta_peaks_w = zeros(num_channels, length(peak_data(i).peaks.Theta.P200.channel(1).max_amplitude));
    beta_peaks_w  = zeros(num_channels, length(peak_data(i).peaks.Beta.P200.channel(1).max_amplitude));
    for ch = 1:num_channels
        theta_peaks_w(ch, :) = peak_data(i).peaks.Theta.P200.channel(ch).max_amplitude;
        beta_peaks_w(ch, :)  = peak_data(i).peaks.Beta.P200.channel(ch).max_amplitude;
    end
    tbr_per_trial_w = theta_peaks_w(:, W_trials) ./ beta_peaks_w(:, W_trials);
    tbr_peaks(:, 1) = mean(tbr_per_trial_w, 2, 'omitnan');
    % T2 segments
    for w = 2:length(all_windows)
        win = all_windows{w};
        theta_peaks = zeros(num_channels, length(peak_data(i).peaks.Theta.(win).channel(1).max_amplitude));
        beta_peaks  = zeros(num_channels, length(peak_data(i).peaks.Beta.(win).channel(1).max_amplitude));
        for ch = 1:num_channels
            theta_peaks(ch, :) = peak_data(i).peaks.Theta.(win).channel(ch).max_amplitude;
            beta_peaks(ch, :)  = peak_data(i).peaks.Beta.(win).channel(ch).max_amplitude;
        end
        tbr_per_trial = theta_peaks(:, T2_trials) ./ beta_peaks(:, T2_trials);
        tbr_peaks(:, w) = mean(tbr_per_trial, 2, 'omitnan');
    end
    for ch = 1:num_channels
        plot(timeCenters, tbr_peaks(ch, :), '-o', 'Color', channel_colors(ch, :), ...
            'LineWidth', 2, 'DisplayName', chanlocs(ch).labels);
    end
    % Shade ERP windows (skip baseline)
    for w = 2:length(all_windows)
        boundsIndex = w - 1;
        tStart = windowBounds(boundsIndex, 1);
        tEnd   = windowBounds(boundsIndex, 2);
        h = fill([tStart, tEnd, tEnd, tStart], [min(tbr_peaks(:))-0.1, min(tbr_peaks(:))-0.1, max(tbr_peaks(:))+0.1, max(tbr_peaks(:))+0.1], ...
            [0.8, 0.9, 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
    end
    xticks(timeCenters);
    xticklabels(all_windows);
    xlabel('Time (ms)');
    ylabel('Peak TBR');
    title(['TBR Peaks Across ERP Windows - ' participant_id]);
    grid on;
    legend('Location','best');
    hold off;
    saveas(gcf, ['H3_PeakGraph_TBR_' participant_id '.png']);
    close(gcf);
    
    %% HB Peak Graph (per channel) with W baseline
    figure('Name', ['HB Peaks - ' participant_id], 'Position', [100, 100, 800, 600]);
    hold on;
    hb_peaks = zeros(num_channels, length(all_windows));
    % W segment
    hb_peaks_w = zeros(num_channels, length(peak_data(i).peaks.HighBeta.P200.channel(1).max_amplitude));
    for ch = 1:num_channels
        hb_peaks_w(ch, :) = peak_data(i).peaks.HighBeta.P200.channel(ch).max_amplitude;
    end
    hb_peaks(:, 1) = mean(hb_peaks_w(:, W_trials), 2, 'omitnan');
    % T2 segments
    for w = 2:length(all_windows)
        win = all_windows{w};
        hb_peaks_all = zeros(num_channels, length(peak_data(i).peaks.HighBeta.(win).channel(1).max_amplitude));
        for ch = 1:num_channels
            hb_peaks_all(ch, :) = peak_data(i).peaks.HighBeta.(win).channel(ch).max_amplitude;
        end
        hb_peaks(:, w) = mean(hb_peaks_all(:, T2_trials), 2, 'omitnan');
    end
    for ch = 1:num_channels
        plot(timeCenters, hb_peaks(ch, :), '-o', 'Color', channel_colors(ch, :), ...
            'LineWidth', 2, 'DisplayName', chanlocs(ch).labels);
    end
    % Shade ERP windows (skip baseline)
    for w = 2:length(all_windows)
        boundsIndex = w - 1;
        tStart = windowBounds(boundsIndex, 1);
        tEnd   = windowBounds(boundsIndex, 2);
        h = fill([tStart, tEnd, tEnd, tStart], [min(hb_peaks(:))-0.1, min(hb_peaks(:))-0.1, max(hb_peaks(:))+0.1, max(hb_peaks(:))+0.1], ...
            [0.8, 0.9, 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
    end
    xticks(timeCenters);
    xticklabels(all_windows);
    xlabel('Time (ms)');
    ylabel('Peak Z-Score (μV)');
    title(['HB Peaks Across ERP Windows - ' participant_id]);
    grid on;
    legend('Location','best');
    hold off;
    saveas(gcf, ['H3_PeakGraph_HB_' participant_id '.png']);
    close(gcf);
end

%% --- Z-Score Graphs for Individuals ---
for i = 1:num_participants
    participant_id = participants{i};
    
    %% FAA Z-Score Graph with W baseline
    figure('Name', ['FAA Z-Scores - ' participant_id], 'Position', [100, 100, 800, 600]);
    hold on;
    faa_zscores = zeros(1, length(all_windows));
    % W segment
    alpha_peaks_w = zeros(num_channels, length(peak_data(i).peaks.Alpha.P200.channel(1).max_amplitude));
    for ch = 1:num_channels
        alpha_peaks_w(ch, :) = peak_data(i).peaks.Alpha.P200.channel(ch).max_amplitude;
    end
    f3_z_w = mean(alpha_peaks_w(f3_idx, W_trials), 'omitnan');
    f4_z_w = mean(alpha_peaks_w(f4_idx, W_trials), 'omitnan');
    faa_zscores(1) = f4_z_w - f3_z_w;  % FAA Z-Score = F4 - F3
    % T2 segments
    for w = 2:length(all_windows)
        win = all_windows{w};
        alpha_peaks = zeros(num_channels, length(peak_data(i).peaks.Alpha.(win).channel(1).max_amplitude));
        for ch = 1:num_channels
            alpha_peaks(ch, :) = peak_data(i).peaks.Alpha.(win).channel(ch).max_amplitude;
        end
        f3_z = mean(alpha_peaks(f3_idx, T2_trials), 'omitnan');
        f4_z = mean(alpha_peaks(f4_idx, T2_trials), 'omitnan');
        faa_zscores(w) = f4_z - f3_z;
    end
    plot(timeCenters, faa_zscores, '-o', 'Color', [1, 0.5, 0], 'LineWidth', 2);
    % Shade ERP windows (skip baseline)
    for w = 2:length(all_windows)
        boundsIndex = w - 1;
        tStart = windowBounds(boundsIndex, 1);
        tEnd   = windowBounds(boundsIndex, 2);
        h = fill([tStart, tEnd, tEnd, tStart], [min(faa_zscores)-0.1, min(faa_zscores)-0.1, max(faa_zscores)+0.1, max(faa_zscores)+0.1], ...
            [0.8, 0.9, 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
    end
    xticks(timeCenters);
    xticklabels(all_windows);
    xlabel('Time (ms)');
    ylabel('FAA Z-Score');
    title(['FAA Z-Score Across ERP Windows - ' participant_id]);
    grid on;
    hold off;
    saveas(gcf, ['H3_ZScoreGraph_FAA_' participant_id '.png']);
    close(gcf);
    faa_zscores_all(i, :) = faa_zscores;  % Store for stats
    
    %% BAR Z-Score Graph (average across channels) with W baseline
    figure('Name', ['BAR Z-Scores - ' participant_id], 'Position', [100, 100, 800, 600]);
    hold on;
    bar_zscores = zeros(1, length(all_windows));
    % W segment
    alpha_peaks_w = zeros(num_channels, length(peak_data(i).peaks.Alpha.P200.channel(1).max_amplitude));
    beta_peaks_w  = zeros(num_channels, length(peak_data(i).peaks.Beta.P200.channel(1).max_amplitude));
    for ch = 1:num_channels
        alpha_peaks_w(ch, :) = peak_data(i).peaks.Alpha.P200.channel(ch).max_amplitude;
        beta_peaks_w(ch, :)  = peak_data(i).peaks.Beta.P200.channel(ch).max_amplitude;
    end
    bar_per_trial_w = beta_peaks_w(:, W_trials) ./ alpha_peaks_w(:, W_trials);
    bar_T2_w = mean(bar_per_trial_w, 2, 'omitnan');
    bar_zscores(1) = mean(bar_T2_w, 'omitnan');
    % T2 segments
    for w = 2:length(all_windows)
        win = all_windows{w};
        alpha_peaks = zeros(num_channels, length(peak_data(i).peaks.Alpha.(win).channel(1).max_amplitude));
        beta_peaks  = zeros(num_channels, length(peak_data(i).peaks.Beta.(win).channel(1).max_amplitude));
        for ch = 1:num_channels
            alpha_peaks(ch, :) = peak_data(i).peaks.Alpha.(win).channel(ch).max_amplitude;
            beta_peaks(ch, :)  = peak_data(i).peaks.Beta.(win).channel(ch).max_amplitude;
        end
        bar_per_trial = beta_peaks(:, T2_trials) ./ alpha_peaks(:, T2_trials);
        bar_T2 = mean(bar_per_trial, 2, 'omitnan');
        bar_zscores(w) = mean(bar_T2, 'omitnan');
    end
    plot(timeCenters, bar_zscores, '-o', 'Color', [1, 0.5, 0], 'LineWidth', 2);
    % Shade ERP windows (skip baseline)
    for w = 2:length(all_windows)
        boundsIndex = w - 1;
        tStart = windowBounds(boundsIndex, 1);
        tEnd   = windowBounds(boundsIndex, 2);
        h = fill([tStart, tEnd, tEnd, tStart], [min(bar_zscores)-0.1, min(bar_zscores)-0.1, max(bar_zscores)+0.1, max(bar_zscores)+0.1], ...
            [0.8, 0.9, 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
    end
    xticks(timeCenters);
    xticklabels(all_windows);
    xlabel('Time (ms)');
    ylabel('BAR Z-Score');
    title(['BAR Z-Score Across ERP Windows - ' participant_id]);
    grid on;
    hold off;
    saveas(gcf, ['H3_ZScoreGraph_BAR_' participant_id '.png']);
    close(gcf);
    bar_zscores_all(i, :) = bar_zscores;
    
    %% TBR Z-Score Graph (average across channels) with W baseline
    figure('Name', ['TBR Z-Scores - ' participant_id], 'Position', [100, 100, 800, 600]);
    hold on;
    tbr_zscores = zeros(1, length(all_windows));
    % W segment
    theta_peaks_w = zeros(num_channels, length(peak_data(i).peaks.Theta.P200.channel(1).max_amplitude));
    beta_peaks_w  = zeros(num_channels, length(peak_data(i).peaks.Beta.P200.channel(1).max_amplitude));
    for ch = 1:num_channels
        theta_peaks_w(ch, :) = peak_data(i).peaks.Theta.P200.channel(ch).max_amplitude;
        beta_peaks_w(ch, :)  = peak_data(i).peaks.Beta.P200.channel(ch).max_amplitude;
    end
    tbr_per_trial_w = theta_peaks_w(:, W_trials) ./ beta_peaks_w(:, W_trials);
    tbr_T2_w = mean(tbr_per_trial_w, 2, 'omitnan');
    tbr_zscores(1) = mean(tbr_T2_w, 'omitnan');
    % T2 segments
    for w = 2:length(all_windows)
        win = all_windows{w};
        theta_peaks = zeros(num_channels, length(peak_data(i).peaks.Theta.(win).channel(1).max_amplitude));
        beta_peaks  = zeros(num_channels, length(peak_data(i).peaks.Beta.(win).channel(1).max_amplitude));
        for ch = 1:num_channels
            theta_peaks(ch, :) = peak_data(i).peaks.Theta.(win).channel(ch).max_amplitude;
            beta_peaks(ch, :)  = peak_data(i).peaks.Beta.(win).channel(ch).max_amplitude;
        end
        tbr_per_trial = theta_peaks(:, T2_trials) ./ beta_peaks(:, T2_trials);
        tbr_T2 = mean(tbr_per_trial, 2, 'omitnan');
        tbr_zscores(w) = mean(tbr_T2, 'omitnan');
    end
    plot(timeCenters, tbr_zscores, '-o', 'Color', [1, 0.5, 0], 'LineWidth', 2);
    % Shade ERP windows (skip baseline)
    for w = 2:length(all_windows)
        boundsIndex = w - 1;
        tStart = windowBounds(boundsIndex, 1);
        tEnd   = windowBounds(boundsIndex, 2);
        h = fill([tStart, tEnd, tEnd, tStart], [min(tbr_zscores)-0.1, min(tbr_zscores)-0.1, max(tbr_zscores)+0.1, max(tbr_zscores)+0.1], ...
            [0.8, 0.9, 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
    end
    xticks(timeCenters);
    xticklabels(all_windows);
    xlabel('Time (ms)');
    ylabel('TBR Z-Score');
    title(['TBR Z-Score Across ERP Windows - ' participant_id]);
    grid on;
    hold off;
    saveas(gcf, ['H3_ZScoreGraph_TBR_' participant_id '.png']);
    close(gcf);
    tbr_zscores_all(i, :) = tbr_zscores;
    
    %% HB Z-Score Graph (average across channels) with W baseline
    figure('Name', ['HB Z-Scores - ' participant_id], 'Position', [100, 100, 800, 600]);
    hold on;
    hb_zscores = zeros(1, length(all_windows));
    % W segment
    hb_peaks_w = zeros(num_channels, length(peak_data(i).peaks.HighBeta.P200.channel(1).max_amplitude));
    for ch = 1:num_channels
        hb_peaks_w(ch, :) = peak_data(i).peaks.HighBeta.P200.channel(ch).max_amplitude;
    end
    hb_T2_w = mean(hb_peaks_w(:, W_trials), 2, 'omitnan');
    hb_zscores(1) = mean(hb_T2_w, 'omitnan');
    % T2 segments
    for w = 2:length(all_windows)
        win = all_windows{w};
        hb_peaks_all = zeros(num_channels, length(peak_data(i).peaks.HighBeta.(win).channel(1).max_amplitude));
        for ch = 1:num_channels
            hb_peaks_all(ch, :) = peak_data(i).peaks.HighBeta.(win).channel(ch).max_amplitude;
        end
        hb_T2 = mean(hb_peaks_all(:, T2_trials), 2, 'omitnan');
        hb_zscores(w) = mean(hb_T2, 'omitnan');
    end
    plot(timeCenters, hb_zscores, '-o', 'Color', [1, 0.5, 0], 'LineWidth', 2);
    % Shade ERP windows (skip baseline)
    for w = 2:length(all_windows)
        boundsIndex = w - 1;
        tStart = windowBounds(boundsIndex, 1);
        tEnd   = windowBounds(boundsIndex, 2);
        h = fill([tStart, tEnd, tEnd, tStart], [min(hb_zscores)-0.1, min(hb_zscores)-0.1, max(hb_zscores)+0.1, max(hb_zscores)+0.1], ...
            [0.8, 0.9, 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
    end
    xticks(timeCenters);
    xticklabels(all_windows);
    xlabel('Time (ms)');
    ylabel('HB Z-Score');
    title(['HB Z-Score Across ERP Windows - ' participant_id]);
    grid on;
    hold off;
    saveas(gcf, ['H3_ZScoreGraph_HB_' participant_id '.png']);
    close(gcf);
    hb_zscores_all(i, :) = hb_zscores;
end

%% --- Peak Graphs for Y* Group (PT2YY, PT3YY, PT4YY, PT7YN) ---
% FAA Peak Graph (Y* Group)
figure('Name', 'FAA Peaks - Y* Group', 'Position', [100, 100, 800, 600]);
hold on;
f3_peaks_y = zeros(length(y_group_indices), length(all_windows));
f4_peaks_y = zeros(length(y_group_indices), length(all_windows));
for p = 1:length(y_group_indices)
    i = y_group_indices(p);
    % W segment
    alpha_peaks_w = zeros(num_channels, length(peak_data(i).peaks.Alpha.P200.channel(1).max_amplitude));
    for ch = 1:num_channels
        alpha_peaks_w(ch, :) = peak_data(i).peaks.Alpha.P200.channel(ch).max_amplitude;
    end
    f3_peaks_y(p, 1) = mean(alpha_peaks_w(f3_idx, W_trials), 'omitnan');
    f4_peaks_y(p, 1) = mean(alpha_peaks_w(f4_idx, W_trials), 'omitnan');
    % T2 segments
    for w = 2:length(all_windows)
        win = all_windows{w};
        alpha_peaks = zeros(num_channels, length(peak_data(i).peaks.Alpha.(win).channel(1).max_amplitude));
        for ch = 1:num_channels
            alpha_peaks(ch, :) = peak_data(i).peaks.Alpha.(win).channel(ch).max_amplitude;
        end
        f3_peaks_y(p, w) = mean(alpha_peaks(f3_idx, T2_trials), 'omitnan');
        f4_peaks_y(p, w) = mean(alpha_peaks(f4_idx, T2_trials), 'omitnan');
    end
end
f3_peaks_y_mean = mean(f3_peaks_y, 1, 'omitnan');
f4_peaks_y_mean = mean(f4_peaks_y, 1, 'omitnan');
plot(timeCenters, f3_peaks_y_mean, '-o', 'Color', [1, 0.5, 0], 'LineWidth', 2, 'DisplayName', 'F3');
plot(timeCenters, f4_peaks_y_mean, '-o', 'Color', [0, 0.5, 1], 'LineWidth', 2, 'DisplayName', 'F4');
% Shade ERP windows (skip baseline)
for w = 2:length(all_windows)
    boundsIndex = w - 1;
    tStart = windowBounds(boundsIndex, 1);
    tEnd   = windowBounds(boundsIndex, 2);
    h = fill([tStart, tEnd, tEnd, tStart], [yMin_default, yMin_default, yMax_default, yMax_default], [0.8, 0.9, 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
end
xticks(timeCenters);
xticklabels(all_windows);
xlabel('Time (ms)');
ylabel('Peak Z-Score (μV)');
title('FAA Peaks Across ERP Components (F3 vs F4) - Y* Group (PT2YY, PT3YY, PT4YY, PT7YN)');
xlim([0, 800]);
ylim([yMin_default, yMax_default]);
grid on;
legend('Location','best');
hold off;
saveas(gcf, 'H3_PeakGraph_FAA_YStar.png');
close(gcf);

%% BAR Peak Graph (Y* Group)
figure('Name', 'BAR Peaks - Y* Group', 'Position', [100, 100, 800, 600]);
hold on;
bar_peaks_y = zeros(length(y_group_indices), num_channels, length(all_windows));
for p = 1:length(y_group_indices)
    i = y_group_indices(p);
    % W segment
    alpha_peaks_w = zeros(num_channels, length(peak_data(i).peaks.Alpha.P200.channel(1).max_amplitude));
    beta_peaks_w  = zeros(num_channels, length(peak_data(i).peaks.Beta.P200.channel(1).max_amplitude));
    for ch = 1:num_channels
        alpha_peaks_w(ch, :) = peak_data(i).peaks.Alpha.P200.channel(ch).max_amplitude;
        beta_peaks_w(ch, :)  = peak_data(i).peaks.Beta.P200.channel(ch).max_amplitude;
    end
    bar_per_trial_w = beta_peaks_w(:, W_trials) ./ alpha_peaks_w(:, W_trials);
    bar_peaks_y(p, :, 1) = mean(bar_per_trial_w, 2, 'omitnan');
    % T2 segments
    for w = 2:length(all_windows)
        win = all_windows{w};
        alpha_peaks = zeros(num_channels, length(peak_data(i).peaks.Alpha.(win).channel(1).max_amplitude));
        beta_peaks  = zeros(num_channels, length(peak_data(i).peaks.Beta.(win).channel(1).max_amplitude));
        for ch = 1:num_channels
            alpha_peaks(ch, :) = peak_data(i).peaks.Alpha.(win).channel(ch).max_amplitude;
            beta_peaks(ch, :)  = peak_data(i).peaks.Beta.(win).channel(ch).max_amplitude;
        end
        bar_per_trial = beta_peaks(:, T2_trials) ./ alpha_peaks(:, T2_trials);
        bar_peaks_y(p, :, w) = mean(bar_per_trial, 2, 'omitnan');
    end
end
bar_peaks_y_mean = mean(bar_peaks_y, 1, 'omitnan');
bar_peaks_y_mean = squeeze(bar_peaks_y_mean);
for ch = 1:num_channels
    plot(timeCenters, bar_peaks_y_mean(ch, :), '-o', 'Color', channel_colors(ch, :), ...
        'LineWidth', 2, 'DisplayName', chanlocs(ch).labels);
end
% Shade ERP windows (skip baseline)
for w = 2:length(all_windows)
    boundsIndex = w - 1;
    tStart = windowBounds(boundsIndex, 1);
    tEnd   = windowBounds(boundsIndex, 2);
    h = fill([tStart, tEnd, tEnd, tStart], [min(bar_peaks_y_mean(:))-0.1, min(bar_peaks_y_mean(:))-0.1, max(bar_peaks_y_mean(:))+0.1, max(bar_peaks_y_mean(:))+0.1], ...
        [0.8, 0.9, 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
end
xticks(timeCenters);
xticklabels(all_windows);
xlabel('Time (ms)');
ylabel('Peak BAR');
title('BAR Peaks Across ERP Windows - Y* Group (PT2YY, PT3YY, PT4YY, PT7YN)');
grid on;
legend('Location','best');
hold off;
saveas(gcf, 'H3_PeakGraph_BAR_YStar.png');
close(gcf);

%% TBR Peak Graph (Y* Group)
figure('Name', 'TBR Peaks - Y* Group', 'Position', [100, 100, 800, 600]);
hold on;
tbr_peaks_y = zeros(length(y_group_indices), num_channels, length(all_windows));
for p = 1:length(y_group_indices)
    i = y_group_indices(p);
    % W segment
    theta_peaks_w = zeros(num_channels, length(peak_data(i).peaks.Theta.P200.channel(1).max_amplitude));
    beta_peaks_w  = zeros(num_channels, length(peak_data(i).peaks.Beta.P200.channel(1).max_amplitude));
    for ch = 1:num_channels
        theta_peaks_w(ch, :) = peak_data(i).peaks.Theta.P200.channel(ch).max_amplitude;
        beta_peaks_w(ch, :)  = peak_data(i).peaks.Beta.P200.channel(ch).max_amplitude;
    end
    tbr_per_trial_w = theta_peaks_w(:, W_trials) ./ beta_peaks_w(:, W_trials);
    tbr_peaks_y(p, :, 1) = mean(tbr_per_trial_w, 2, 'omitnan');
    % T2 segments
    for w = 2:length(all_windows)
        win = all_windows{w};
        theta_peaks = zeros(num_channels, length(peak_data(i).peaks.Theta.(win).channel(1).max_amplitude));
        beta_peaks  = zeros(num_channels, length(peak_data(i).peaks.Beta.(win).channel(1).max_amplitude));
        for ch = 1:num_channels
            theta_peaks(ch, :) = peak_data(i).peaks.Theta.(win).channel(ch).max_amplitude;
            beta_peaks(ch, :)  = peak_data(i).peaks.Beta.(win).channel(ch).max_amplitude;
        end
        tbr_per_trial = theta_peaks(:, T2_trials) ./ beta_peaks(:, T2_trials);
        tbr_peaks_y(p, :, w) = mean(tbr_per_trial, 2, 'omitnan');
    end
end
tbr_peaks_y_mean = mean(tbr_peaks_y, 1, 'omitnan');
tbr_peaks_y_mean = squeeze(tbr_peaks_y_mean);
for ch = 1:num_channels
    plot(timeCenters, tbr_peaks_y_mean(ch, :), '-o', 'Color', channel_colors(ch, :), ...
        'LineWidth', 2, 'DisplayName', chanlocs(ch).labels);
end
% Shade ERP windows (skip baseline)
for w = 2:length(all_windows)
    boundsIndex = w - 1;
    tStart = windowBounds(boundsIndex, 1);
    tEnd   = windowBounds(boundsIndex, 2);
    h = fill([tStart, tEnd, tEnd, tStart], [min(tbr_peaks_y_mean(:))-0.1, min(tbr_peaks_y_mean(:))-0.1, max(tbr_peaks_y_mean(:))+0.1, max(tbr_peaks_y_mean(:))+0.1], ...
        [0.8, 0.9, 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
end
xticks(timeCenters);
xticklabels(all_windows);
xlabel('Time (ms)');
ylabel('Peak TBR');
title('TBR Peaks Across ERP Windows - Y* Group (PT2YY, PT3YY, PT4YY, PT7YN)');
grid on;
legend('Location','best');
hold off;
saveas(gcf, 'H3_PeakGraph_TBR_YStar.png');
close(gcf);

%% HB Peak Graph (Y* Group)
figure('Name', 'HB Peaks - Y* Group', 'Position', [100, 100, 800, 600]);
hold on;
hb_peaks_y = zeros(length(y_group_indices), num_channels, length(all_windows));
for p = 1:length(y_group_indices)
    i = y_group_indices(p);
    % W segment
    hb_peaks_w = zeros(num_channels, length(peak_data(i).peaks.HighBeta.P200.channel(1).max_amplitude));
    for ch = 1:num_channels
        hb_peaks_w(ch, :) = peak_data(i).peaks.HighBeta.P200.channel(ch).max_amplitude;
    end
    hb_peaks_y(p, :, 1) = mean(hb_peaks_w(:, W_trials), 2, 'omitnan');
    % T2 segments
    for w = 2:length(all_windows)
        win = all_windows{w};
        hb_peaks_all = zeros(num_channels, length(peak_data(i).peaks.HighBeta.(win).channel(1).max_amplitude));
        for ch = 1:num_channels
            hb_peaks_all(ch, :) = peak_data(i).peaks.HighBeta.(win).channel(ch).max_amplitude;
        end
        hb_peaks_y(p, :, w) = mean(hb_peaks_all(:, T2_trials), 2, 'omitnan');
    end
end
hb_peaks_y_mean = mean(hb_peaks_y, 1, 'omitnan');
hb_peaks_y_mean = squeeze(hb_peaks_y_mean);
for ch = 1:num_channels
    plot(timeCenters, hb_peaks_y_mean(ch, :), '-o', 'Color', channel_colors(ch, :), ...
        'LineWidth', 2, 'DisplayName', chanlocs(ch).labels);
end
% Shade ERP windows (skip baseline)
for w = 2:length(all_windows)
    boundsIndex = w - 1;
    tStart = windowBounds(boundsIndex, 1);
    tEnd   = windowBounds(boundsIndex, 2);
    h = fill([tStart, tEnd, tEnd, tStart], [min(hb_peaks_y_mean(:))-0.1, min(hb_peaks_y_mean(:))-0.1, max(hb_peaks_y_mean(:))+0.1, max(hb_peaks_y_mean(:))+0.1], ...
        [0.8, 0.9, 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
end
xticks(timeCenters);
xticklabels(all_windows);
xlabel('Time (ms)');
ylabel('Peak Z-Score (μV)');
title('HB Peaks Across ERP Windows - Y* Group (PT2YY, PT3YY, PT4YY, PT7YN)');
grid on;
legend('Location','best');
hold off;
saveas(gcf, 'H3_PeakGraph_HB_YStar.png');
close(gcf);

%% --- Z-Score Graphs for Y* Group ---
% FAA Z-Score Graph (Y* Group)
figure('Name', 'FAA Z-Scores - Y* Group', 'Position', [100, 100, 800, 600]);
hold on;
faa_zscores_y = zeros(length(y_group_indices), length(all_windows));
for p = 1:length(y_group_indices)
    i = y_group_indices(p);
    % W segment
    alpha_peaks_w = zeros(num_channels, length(peak_data(i).peaks.Alpha.P200.channel(1).max_amplitude));
    for ch = 1:num_channels
        alpha_peaks_w(ch, :) = peak_data(i).peaks.Alpha.P200.channel(ch).max_amplitude;
    end
    f3_z_w = mean(alpha_peaks_w(f3_idx, W_trials), 'omitnan');
    f4_z_w = mean(alpha_peaks_w(f4_idx, W_trials), 'omitnan');
    faa_zscores_y(p, 1) = f4_z_w - f3_z_w;
    % T2 segments
    for w = 2:length(all_windows)
        win = all_windows{w};
        alpha_peaks = zeros(num_channels, length(peak_data(i).peaks.Alpha.(win).channel(1).max_amplitude));
        for ch = 1:num_channels
            alpha_peaks(ch, :) = peak_data(i).peaks.Alpha.(win).channel(ch).max_amplitude;
        end
        f3_z = mean(alpha_peaks(f3_idx, T2_trials), 'omitnan');
        f4_z = mean(alpha_peaks(f4_idx, T2_trials), 'omitnan');
        faa_zscores_y(p, w) = f4_z - f3_z;
    end
end
faa_zscores_y_mean = mean(faa_zscores_y, 1, 'omitnan');
plot(timeCenters, faa_zscores_y_mean, '-o', 'Color', [1, 0.5, 0], 'LineWidth', 2);
% Shade ERP windows (skip baseline)
for w = 2:length(all_windows)
    boundsIndex = w - 1;
    tStart = windowBounds(boundsIndex, 1);
    tEnd   = windowBounds(boundsIndex, 2);
    h = fill([tStart, tEnd, tEnd, tStart], [min(faa_zscores_y_mean)-0.1, min(faa_zscores_y_mean)-0.1, max(faa_zscores_y_mean)+0.1, max(faa_zscores_y_mean)+0.1], ...
        [0.8, 0.9, 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
end
xticks(timeCenters);
xticklabels(all_windows);
xlabel('Time (ms)');
ylabel('FAA Z-Score');
title('FAA Z-Score Across ERP Windows - Y* Group (PT2YY, PT3YY, PT4YY, PT7YN)');
grid on;
hold off;
saveas(gcf, 'H3_ZScoreGraph_FAA_YStar.png');
close(gcf);

%% BAR Z-Score Graph (Y* Group)
figure('Name', 'BAR Z-Scores - Y* Group', 'Position', [100, 100, 800, 600]);
hold on;
bar_zscores_y = zeros(length(y_group_indices), length(all_windows));
for p = 1:length(y_group_indices)
    i = y_group_indices(p);
    % W segment
    alpha_peaks_w = zeros(num_channels, length(peak_data(i).peaks.Alpha.P200.channel(1).max_amplitude));
    beta_peaks_w  = zeros(num_channels, length(peak_data(i).peaks.Beta.P200.channel(1).max_amplitude));
    for ch = 1:num_channels
        alpha_peaks_w(ch, :) = peak_data(i).peaks.Alpha.P200.channel(ch).max_amplitude;
        beta_peaks_w(ch, :)  = peak_data(i).peaks.Beta.P200.channel(ch).max_amplitude;
    end
    bar_per_trial_w = beta_peaks_w(:, W_trials) ./ alpha_peaks_w(:, W_trials);
    bar_T2_w = mean(bar_per_trial_w, 2, 'omitnan');
    bar_zscores_y(p, 1) = mean(bar_T2_w, 'omitnan');
    % T2 segments
    for w = 2:length(all_windows)
        win = all_windows{w};
        alpha_peaks = zeros(num_channels, length(peak_data(i).peaks.Alpha.(win).channel(1).max_amplitude));
        beta_peaks  = zeros(num_channels, length(peak_data(i).peaks.Beta.(win).channel(1).max_amplitude));
        for ch = 1:num_channels
            alpha_peaks(ch, :) = peak_data(i).peaks.Alpha.(win).channel(ch).max_amplitude;
            beta_peaks(ch, :)  = peak_data(i).peaks.Beta.(win).channel(ch).max_amplitude;
        end
        bar_per_trial = beta_peaks(:, T2_trials) ./ alpha_peaks(:, T2_trials);
        bar_T2 = mean(bar_per_trial, 2, 'omitnan');
        bar_zscores_y(p, w) = mean(bar_T2, 'omitnan');
    end
end
bar_zscores_y_mean = mean(bar_zscores_y, 1, 'omitnan');
plot(timeCenters, bar_zscores_y_mean, '-o', 'Color', [1, 0.5, 0], 'LineWidth', 2);
% Shade ERP windows (skip baseline)
for w = 2:length(all_windows)
    boundsIndex = w - 1;
    tStart = windowBounds(boundsIndex, 1);
    tEnd   = windowBounds(boundsIndex, 2);
    h = fill([tStart, tEnd, tEnd, tStart], [min(bar_zscores_y_mean)-0.1, min(bar_zscores_y_mean)-0.1, max(bar_zscores_y_mean)+0.1, max(bar_zscores_y_mean)+0.1], ...
        [0.8, 0.9, 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
end
xticks(timeCenters);
xticklabels(all_windows);
xlabel('Time (ms)');
ylabel('BAR Z-Score');
title('BAR Z-Score Across ERP Windows - Y* Group (PT2YY, PT3YY, PT4YY, PT7YN)');
grid on;
hold off;
saveas(gcf, 'H3_ZScoreGraph_BAR_YStar.png');
close(gcf);

%% TBR Z-Score Graph (Y* Group)
figure('Name', 'TBR Z-Scores - Y* Group', 'Position', [100, 100, 800, 600]);
hold on;
tbr_zscores_y = zeros(length(y_group_indices), length(all_windows));
for p = 1:length(y_group_indices)
    i = y_group_indices(p);
    % W segment
    theta_peaks_w = zeros(num_channels, length(peak_data(i).peaks.Theta.P200.channel(1).max_amplitude));
    beta_peaks_w  = zeros(num_channels, length(peak_data(i).peaks.Beta.P200.channel(1).max_amplitude));
    for ch = 1:num_channels
        theta_peaks_w(ch, :) = peak_data(i).peaks.Theta.P200.channel(ch).max_amplitude;
        beta_peaks_w(ch, :)  = peak_data(i).peaks.Beta.P200.channel(ch).max_amplitude;
    end
    tbr_per_trial_w = theta_peaks_w(:, W_trials) ./ beta_peaks_w(:, W_trials);
    temp = mean(tbr_per_trial_w, 2, 'omitnan');
    tbr_zscores_y(p, 1) = mean(temp, 'omitnan');

    %tbr_zscores_y(p, 1) = mean(tbr_per_trial_w, 'omitnan');
    % T2 segments
    for w = 2:length(all_windows)
        win = all_windows{w};
        theta_peaks = zeros(num_channels, length(peak_data(i).peaks.Theta.(win).channel(1).max_amplitude));
        beta_peaks  = zeros(num_channels, length(peak_data(i).peaks.Beta.(win).channel(1).max_amplitude));
        for ch = 1:num_channels
            theta_peaks(ch, :) = peak_data(i).peaks.Theta.(win).channel(ch).max_amplitude;
            beta_peaks(ch, :)  = peak_data(i).peaks.Beta.(win).channel(ch).max_amplitude;
        end
        tbr_per_trial = theta_peaks(:, T2_trials) ./ beta_peaks(:, T2_trials);
        temp = mean(tbr_per_trial, 2, 'omitnan');
        tbr_zscores_y(p, w) = mean(temp, 'omitnan');
    end
end
tbr_zscores_y_mean = mean(tbr_zscores_y, 1, 'omitnan');
plot(timeCenters, tbr_zscores_y_mean, '-o', 'Color', [1, 0.5, 0], 'LineWidth', 2);
% Shade ERP windows (skip baseline)
for w = 2:length(all_windows)
    boundsIndex = w - 1;
    tStart = windowBounds(boundsIndex, 1);
    tEnd   = windowBounds(boundsIndex, 2);
    h = fill([tStart, tEnd, tEnd, tStart], [min(tbr_zscores_y_mean)-0.1, min(tbr_zscores_y_mean)-0.1, max(tbr_zscores_y_mean)+0.1, max(tbr_zscores_y_mean)+0.1], ...
        [0.8, 0.9, 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
end
xticks(timeCenters);
xticklabels(all_windows);
xlabel('Time (ms)');
ylabel('TBR Z-Score');
title('TBR Z-Score Across ERP Windows - Y* Group (PT2YY, PT3YY, PT4YY, PT7YN)');
grid on;
hold off;
saveas(gcf, 'H3_ZScoreGraph_TBR_YStar.png');
close(gcf);

%% HB Z-Score Graph (Y* Group)
figure('Name', 'HB Z-Scores - Y* Group', 'Position', [100, 100, 800, 600]);
hold on;
hb_zscores_y = zeros(length(y_group_indices), length(all_windows));
for p = 1:length(y_group_indices)
    i = y_group_indices(p);
    % W segment
    hb_peaks_w = zeros(num_channels, length(peak_data(i).peaks.HighBeta.P200.channel(1).max_amplitude));
    for ch = 1:num_channels
        hb_peaks_w(ch, :) = peak_data(i).peaks.HighBeta.P200.channel(ch).max_amplitude;
    end
    temp = mean(hb_peaks_w(:, W_trials), 'omitnan');
    hb_zscores_y(p, 1) =  mean(temp, 'omitnan');
    % T2 segments
    for w = 2:length(all_windows)
        win = all_windows{w};
        hb_peaks_all = zeros(num_channels, length(peak_data(i).peaks.HighBeta.(win).channel(1).max_amplitude));
        for ch = 1:num_channels
            hb_peaks_all(ch, :) = peak_data(i).peaks.HighBeta.(win).channel(ch).max_amplitude;
        end
        temp = mean(hb_peaks_all(:, T2_trials), 'omitnan')
        hb_zscores_y(p, w) = mean(temp, 'omitnan') ;
    end
end
hb_zscores_y_mean = mean(hb_zscores_y, 1, 'omitnan');
plot(timeCenters, hb_zscores_y_mean, '-o', 'Color', [1, 0.5, 0], 'LineWidth', 2);
% Shade ERP windows (skip baseline)
for w = 2:length(all_windows)
    boundsIndex = w - 1;
    tStart = windowBounds(boundsIndex, 1);
    tEnd   = windowBounds(boundsIndex, 2);
    h = fill([tStart, tEnd, tEnd, tStart], [min(hb_zscores_y_mean)-0.1, min(hb_zscores_y_mean)-0.1, max(hb_zscores_y_mean)+0.1, max(hb_zscores_y_mean)+0.1], ...
        [0.8, 0.9, 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
end
xticks(timeCenters);
xticklabels(all_windows);
xlabel('Time (ms)');
ylabel('HB Z-Score');
title('HB Z-Score Across ERP Windows - Y* Group (PT2YY, PT3YY, PT4YY, PT7YN)');
grid on;
hold off;
saveas(gcf, 'H3_ZScoreGraph_HB_YStar.png');
close(gcf);

%% --- Descriptive Statistics and Spearman’s Correlation ---
% Compute mean z-scores for T2 (average across ERP windows P200, N200, P300, LPP) and W
T2_windows = 2:length(all_windows);  % P200, N200, P300, LPP
faa_zscores_T2 = mean(faa_zscores_all(:, T2_windows), 2, 'omitnan');  % Mean across T2 windows
faa_zscores_w  = faa_zscores_all(:, 1);  % W segment
bar_zscores_T2 = mean(bar_zscores_all(:, T2_windows), 2, 'omitnan');
bar_zscores_w  = bar_zscores_all(:, 1);
tbr_zscores_T2 = mean(tbr_zscores_all(:, T2_windows), 2, 'omitnan');
tbr_zscores_w  = tbr_zscores_all(:, 1);
hb_zscores_T2  = mean(hb_zscores_all(:, T2_windows), 2, 'omitnan');
hb_zscores_w   = hb_zscores_all(:, 1);

% Descriptive Statistics
desc_stats = table();
desc_stats.Metric = {'FAA'; 'FAA'; 'BAR'; 'BAR'; 'TBR'; 'TBR'; 'HB'; 'HB'};
desc_stats.Condition = {'T2'; 'W'; 'T2'; 'W'; 'T2'; 'W'; 'T2'; 'W'};
desc_stats.Mean = [
    mean(faa_zscores_T2, 'omitnan');
    mean(faa_zscores_w, 'omitnan');
    mean(bar_zscores_T2, 'omitnan');
    mean(bar_zscores_w, 'omitnan');
    mean(tbr_zscores_T2, 'omitnan');
    mean(tbr_zscores_w, 'omitnan');
    mean(hb_zscores_T2, 'omitnan');
    mean(hb_zscores_w, 'omitnan')
];
desc_stats.SD = [
    std(faa_zscores_T2, 'omitnan');
    std(faa_zscores_w, 'omitnan');
    std(bar_zscores_T2, 'omitnan');
    std(bar_zscores_w, 'omitnan');
    std(tbr_zscores_T2, 'omitnan');
    std(tbr_zscores_w, 'omitnan');
    std(hb_zscores_T2, 'omitnan');
    std(hb_zscores_w, 'omitnan')
];
disp('Descriptive Statistics (All Participants):');
disp(desc_stats);
writetable(desc_stats, 'H3_Descriptive_Statistics.csv');

% Spearman’s Correlation (T2 vs W)
[rho_faa, p_faa] = corr(faa_zscores_T2, faa_zscores_w, 'Type', 'Spearman', 'Rows', 'complete');
[rho_bar, p_bar] = corr(bar_zscores_T2, bar_zscores_w, 'Type', 'Spearman', 'Rows', 'complete');
[rho_tbr, p_tbr] = corr(tbr_zscores_T2, tbr_zscores_w, 'Type', 'Spearman', 'Rows', 'complete');
[rho_hb, p_hb]   = corr(hb_zscores_T2, hb_zscores_w, 'Type', 'Spearman', 'Rows', 'complete');

spearman_stats = table();
spearman_stats.Metric = {'FAA'; 'BAR'; 'TBR'; 'HB'};
spearman_stats.Rho = [rho_faa; rho_bar; rho_tbr; rho_hb];
spearman_stats.PValue = [p_faa; p_bar; p_tbr; p_hb];
disp('Spearman’s Correlation (T2 vs W, All Participants):');
disp(spearman_stats);
writetable(spearman_stats, 'H3_Spearman_Correlation.csv');

%% --- Cohen’s d for Y* Group (and PT2NN for discussion) ---
% Y* Group (PT2YY, PT3YY, PT4YY, PT7YN)
faa_zscores_y_T2 = mean(faa_zscores_all(y_group_indices, T2_windows), 2, 'omitnan');
faa_zscores_y_w  = faa_zscores_all(y_group_indices, 1);
bar_zscores_y_T2 = mean(bar_zscores_all(y_group_indices, T2_windows), 2, 'omitnan');
bar_zscores_y_w  = bar_zscores_all(y_group_indices, 1);
tbr_zscores_y_T2 = mean(tbr_zscores_all(y_group_indices, T2_windows), 2, 'omitnan');
tbr_zscores_y_w  = tbr_zscores_all(y_group_indices, 1);
hb_zscores_y_T2  = mean(hb_zscores_all(y_group_indices, T2_windows), 2, 'omitnan');
hb_zscores_y_w   = hb_zscores_all(y_group_indices, 1);

% PT2NN (for discussion, N=1)
faa_zscores_pT2_T2 = mean(faa_zscores_all(1, T2_windows), 'omitnan');
faa_zscores_pT2_w  = faa_zscores_all(1, 1);
bar_zscores_pT2_T2 = mean(bar_zscores_all(1, T2_windows), 'omitnan');
bar_zscores_pT2_w  = bar_zscores_all(1, 1);
tbr_zscores_pT2_T2 = mean(tbr_zscores_all(1, T2_windows), 'omitnan');
tbr_zscores_pT2_w  = tbr_zscores_all(1, 1);
hb_zscores_pT2_T2  = mean(hb_zscores_all(1, T2_windows), 'omitnan');
hb_zscores_pT2_w   = hb_zscores_all(1, 1);

% Cohen’s d for Y* Group
cohens_d_y = zeros(4, 1);
% FAA
mean_diff_faa_y = mean(faa_zscores_y_T2) - mean(faa_zscores_y_w);
pooled_sd_faa_y = sqrt((std(faa_zscores_y_T2)^2 + std(faa_zscores_y_w)^2) / 2);
cohens_d_y(1) = mean_diff_faa_y / pooled_sd_faa_y;
% BAR
mean_diff_bar_y = mean(bar_zscores_y_T2) - mean(bar_zscores_y_w);
pooled_sd_bar_y = sqrt((std(bar_zscores_y_T2)^2 + std(bar_zscores_y_w)^2) / 2);
cohens_d_y(2) = mean_diff_bar_y / pooled_sd_bar_y;
% TBR
mean_diff_tbr_y = mean(tbr_zscores_y_T2) - mean(tbr_zscores_y_w);
pooled_sd_tbr_y = sqrt((std(tbr_zscores_y_T2)^2 + std(tbr_zscores_y_w)^2) / 2);
cohens_d_y(3) = mean_diff_tbr_y / pooled_sd_tbr_y;
% HB
mean_diff_hb_y = mean(hb_zscores_y_T2) - mean(hb_zscores_y_w);
pooled_sd_hb_y = sqrt((std(hb_zscores_y_T2)^2 + std(hb_zscores_y_w)^2) / 2);
cohens_d_y(4) = mean_diff_hb_y / pooled_sd_hb_y;

% Cohen’s d for PT2NN (N=1, for discussion only)
cohens_d_pT2 = zeros(4, 1);
% FAA
mean_diff_faa_pT2 = faa_zscores_pT2_T2 - faa_zscores_pT2_w;
pooled_sd_faa_pT2 = pooled_sd_faa_y;  % Using Y* group SD as proxy
cohens_d_pT2(1) = mean_diff_faa_pT2 / pooled_sd_faa_pT2;
% BAR
mean_diff_bar_pT2 = bar_zscores_pT2_T2 - bar_zscores_pT2_w;
pooled_sd_bar_pT2 = pooled_sd_bar_y;
cohens_d_pT2(2) = mean_diff_bar_pT2 / pooled_sd_bar_pT2;
% TBR
mean_diff_tbr_pT2 = tbr_zscores_pT2_T2 - tbr_zscores_pT2_w;
pooled_sd_tbr_pT2 = pooled_sd_tbr_y;
cohens_d_pT2(3) = mean_diff_tbr_pT2 / pooled_sd_tbr_pT2;
% HB
mean_diff_hb_pT2 = hb_zscores_pT2_T2 - hb_zscores_pT2_w;
pooled_sd_hb_pT2 = pooled_sd_hb_y;
cohens_d_pT2(4) = mean_diff_hb_pT2 / pooled_sd_hb_pT2;

% Cohen’s d Table
cohens_d_table = table();
cohens_d_table.Metric = {'FAA'; 'BAR'; 'TBR'; 'HB'; 'FAA_PT2NN'; 'BAR_PT2NN'; 'TBR_PT2NN'; 'HB_PT2NN'};
cohens_d_table.Group  = [cohens_d_y; cohens_d_pT2];
disp('Cohens (T2 vs W, All Participants):');
disp(cohens_d_table);
writetable(cohens_d_table, 'H3_Cohens_Correlation.csv');




%% H4 - Peak and Z-Score Graphs for All Participants and Y* Group with Baseline, Statistics, and Cohen's d
% Define parameters
T1T2_trials = segments.T1T2;      % Trial indices for T1T2 condition (241:300)
W_trials = segments.W;        % Trial indices for W condition (baseline)
chanlocs = all_eeg(1).EEG.chanlocs;  % Channel locations
num_channels = length(chanlocs);     % Number of EEG channels (assumed 8: AF3, AF4, F3, F4, FC5, FC6, P8, F8)
num_participants = length(participants);  % Number of participants (5)
all_windows = {'W', 'P200', 'N200', 'P300', 'LPP'};  % ERP windows including W as baseline
metrics = {'FAA', 'BAR', 'TBR', 'HB'};
y_group_indices = [2, 3, 4, 5];  % Indices for Y* group: PT1T2YY, PT3YY, PT4YY, PT7YN

% Find F3 and F4 channel indices (for FAA)
f3_idx = find(strcmp({chanlocs.labels}, 'F3'));
f4_idx = find(strcmp({chanlocs.labels}, 'F4'));
if isempty(f3_idx) || isempty(f4_idx)
    error('F3 or F4 channel not found in chanlocs.');
end

% Define colors for channels (for peak graphs with multiple channels)
channel_colors = lines(num_channels);  % MATLAB's 'lines' colormap for distinct colors

% Define time mapping for ERP windows (ms)
% Baseline (W): assume 0–100 ms with midpoint = 50 ms
% P200: 150–250 ms, midpoint = 200 ms
% N200: 180–240 ms, midpoint = 210 ms
% P300: 300–400 ms, midpoint = 350 ms
% LPP: 400–800 ms, midpoint = 600 ms
timeCenters = [50, 175, 250, 350, 510];  % one value per window (W, P200, N200, P300, LPP)
% For shading, we use the ERP windows (excluding baseline, index 1)
windowBounds = [150 200;   % P200 bounds
                210 290;   % N200 bounds
                310 390;   % P300 bounds
                420 600];  % LPP bounds
yMin_default = -1;
yMax_default = 2;

% Initialize arrays to store z-scores for statistical analysis
faa_zscores_all = zeros(num_participants, length(all_windows));  % [participants, windows]
bar_zscores_all = zeros(num_participants, length(all_windows));
tbr_zscores_all = zeros(num_participants, length(all_windows));
hb_zscores_all = zeros(num_participants, length(all_windows));

%% --- Peak Graphs for Individuals ---
for i = 1:num_participants
    participant_id = participants{i};
    
    %% FAA Peak Graph (F3 vs F4) with W baseline
    figure('Name', ['FAA Peaks - ' participant_id], 'Position', [100, 100, 800, 600]);
    hold on;
    f3_peaks = zeros(1, length(all_windows));
    f4_peaks = zeros(1, length(all_windows));
    % W segment (baseline)
    alpha_peaks_w = zeros(num_channels, length(peak_data(i).peaks.Alpha.P200.channel(1).max_amplitude));
    for ch = 1:num_channels
        alpha_peaks_w(ch, :) = peak_data(i).peaks.Alpha.P200.channel(ch).max_amplitude;
    end
    f3_peaks(1) = mean(alpha_peaks_w(f3_idx, W_trials), 'omitnan');
    f4_peaks(1) = mean(alpha_peaks_w(f4_idx, W_trials), 'omitnan');
    % T1T2 segments
    for w = 2:length(all_windows)
        win = all_windows{w};
        alpha_peaks = zeros(num_channels, length(peak_data(i).peaks.Alpha.(win).channel(1).max_amplitude));
        for ch = 1:num_channels
            alpha_peaks(ch, :) = peak_data(i).peaks.Alpha.(win).channel(ch).max_amplitude;
        end
        f3_peaks(w) = mean(alpha_peaks(f3_idx, T1T2_trials), 'omitnan');
        f4_peaks(w) = mean(alpha_peaks(f4_idx, T1T2_trials), 'omitnan');
    end
    plot(timeCenters, f3_peaks, '-o', 'Color', [1, 0.5, 0], 'LineWidth', 2, 'DisplayName', 'F3');
    plot(timeCenters, f4_peaks, '-o', 'Color', [0, 0.5, 1], 'LineWidth', 2, 'DisplayName', 'F4');
    % Shade ERP windows only (skip baseline: index 1)
    for w = 2:length(all_windows)
        boundsIndex = w - 1;
        tStart = windowBounds(boundsIndex, 1);
        tEnd   = windowBounds(boundsIndex, 2);
        h = fill([tStart, tEnd, tEnd, tStart], [yMin_default, yMin_default, yMax_default, yMax_default], [0.8, 0.9, 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
    end
    xticks(timeCenters);
    xticklabels(all_windows);
    xlabel('Time (ms)');
    ylabel('Peak Z-Score (μV)');
    title(['FAA Peaks Across ERP Components (F3 vs F4) - ' participant_id]);
    xlim([0, 800]);
    ylim([yMin_default, yMax_default]);
    grid on;
    legend('Location','best');
    hold off;
    saveas(gcf, ['H4_PeakGraph_FAA_' participant_id '.png']);
    close(gcf);
    
    %% BAR Peak Graph (per channel) with W baseline
    figure('Name', ['BAR Peaks - ' participant_id], 'Position', [100, 100, 800, 600]);
    hold on;
    bar_peaks = zeros(num_channels, length(all_windows));
    % W segment
    alpha_peaks_w = zeros(num_channels, length(peak_data(i).peaks.Alpha.P200.channel(1).max_amplitude));
    beta_peaks_w = zeros(num_channels, length(peak_data(i).peaks.Beta.P200.channel(1).max_amplitude));
    for ch = 1:num_channels
        alpha_peaks_w(ch, :) = peak_data(i).peaks.Alpha.P200.channel(ch).max_amplitude;
        beta_peaks_w(ch, :) = peak_data(i).peaks.Beta.P200.channel(ch).max_amplitude;
    end
    bar_per_trial_w = beta_peaks_w(:, W_trials) ./ alpha_peaks_w(:, W_trials);
    bar_peaks(:, 1) = mean(bar_per_trial_w, 2, 'omitnan');
    % T1T2 segments
    for w = 2:length(all_windows)
        win = all_windows{w};
        alpha_peaks = zeros(num_channels, length(peak_data(i).peaks.Alpha.(win).channel(1).max_amplitude));
        beta_peaks  = zeros(num_channels, length(peak_data(i).peaks.Beta.(win).channel(1).max_amplitude));
        for ch = 1:num_channels
            alpha_peaks(ch, :) = peak_data(i).peaks.Alpha.(win).channel(ch).max_amplitude;
            beta_peaks(ch, :)  = peak_data(i).peaks.Beta.(win).channel(ch).max_amplitude;
        end
        bar_per_trial = beta_peaks(:, T1T2_trials) ./ alpha_peaks(:, T1T2_trials);
        bar_peaks(:, w) = mean(bar_per_trial, 2, 'omitnan');
    end
    for ch = 1:num_channels
        plot(timeCenters, bar_peaks(ch, :), '-o', 'Color', channel_colors(ch, :), ...
            'LineWidth', 2, 'DisplayName', chanlocs(ch).labels);
    end
    % Shade ERP windows (skip baseline)
    for w = 2:length(all_windows)
        boundsIndex = w - 1;
        tStart = windowBounds(boundsIndex, 1);
        tEnd   = windowBounds(boundsIndex, 2);
        h = fill([tStart, tEnd, tEnd, tStart], [min(bar_peaks(:))-0.1, min(bar_peaks(:))-0.1, max(bar_peaks(:))+0.1, max(bar_peaks(:))+0.1], ...
            [0.8, 0.9, 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
    end
    xticks(timeCenters);
    xticklabels(all_windows);
    xlabel('Time (ms)');
    ylabel('Peak BAR');
    title(['BAR Peaks Across ERP Windows - ' participant_id]);
    grid on;
    legend('Location','best');
    hold off;
    saveas(gcf, ['H4_PeakGraph_BAR_' participant_id '.png']);
    close(gcf);
    
    %% TBR Peak Graph (per channel) with W baseline
    figure('Name', ['TBR Peaks - ' participant_id], 'Position', [100, 100, 800, 600]);
    hold on;
    tbr_peaks = zeros(num_channels, length(all_windows));
    % W segment
    theta_peaks_w = zeros(num_channels, length(peak_data(i).peaks.Theta.P200.channel(1).max_amplitude));
    beta_peaks_w  = zeros(num_channels, length(peak_data(i).peaks.Beta.P200.channel(1).max_amplitude));
    for ch = 1:num_channels
        theta_peaks_w(ch, :) = peak_data(i).peaks.Theta.P200.channel(ch).max_amplitude;
        beta_peaks_w(ch, :)  = peak_data(i).peaks.Beta.P200.channel(ch).max_amplitude;
    end
    tbr_per_trial_w = theta_peaks_w(:, W_trials) ./ beta_peaks_w(:, W_trials);
    tbr_peaks(:, 1) = mean(tbr_per_trial_w, 2, 'omitnan');
    % T1T2 segments
    for w = 2:length(all_windows)
        win = all_windows{w};
        theta_peaks = zeros(num_channels, length(peak_data(i).peaks.Theta.(win).channel(1).max_amplitude));
        beta_peaks  = zeros(num_channels, length(peak_data(i).peaks.Beta.(win).channel(1).max_amplitude));
        for ch = 1:num_channels
            theta_peaks(ch, :) = peak_data(i).peaks.Theta.(win).channel(ch).max_amplitude;
            beta_peaks(ch, :)  = peak_data(i).peaks.Beta.(win).channel(ch).max_amplitude;
        end
        tbr_per_trial = theta_peaks(:, T1T2_trials) ./ beta_peaks(:, T1T2_trials);
        tbr_peaks(:, w) = mean(tbr_per_trial, 2, 'omitnan');
    end
    for ch = 1:num_channels
        plot(timeCenters, tbr_peaks(ch, :), '-o', 'Color', channel_colors(ch, :), ...
            'LineWidth', 2, 'DisplayName', chanlocs(ch).labels);
    end
    % Shade ERP windows (skip baseline)
    for w = 2:length(all_windows)
        boundsIndex = w - 1;
        tStart = windowBounds(boundsIndex, 1);
        tEnd   = windowBounds(boundsIndex, 2);
        h = fill([tStart, tEnd, tEnd, tStart], [min(tbr_peaks(:))-0.1, min(tbr_peaks(:))-0.1, max(tbr_peaks(:))+0.1, max(tbr_peaks(:))+0.1], ...
            [0.8, 0.9, 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
    end
    xticks(timeCenters);
    xticklabels(all_windows);
    xlabel('Time (ms)');
    ylabel('Peak TBR');
    title(['TBR Peaks Across ERP Windows - ' participant_id]);
    grid on;
    legend('Location','best');
    hold off;
    saveas(gcf, ['H4_PeakGraph_TBR_' participant_id '.png']);
    close(gcf);
    
    %% HB Peak Graph (per channel) with W baseline
    figure('Name', ['HB Peaks - ' participant_id], 'Position', [100, 100, 800, 600]);
    hold on;
    hb_peaks = zeros(num_channels, length(all_windows));
    % W segment
    hb_peaks_w = zeros(num_channels, length(peak_data(i).peaks.HighBeta.P200.channel(1).max_amplitude));
    for ch = 1:num_channels
        hb_peaks_w(ch, :) = peak_data(i).peaks.HighBeta.P200.channel(ch).max_amplitude;
    end
    hb_peaks(:, 1) = mean(hb_peaks_w(:, W_trials), 2, 'omitnan');
    % T1T2 segments
    for w = 2:length(all_windows)
        win = all_windows{w};
        hb_peaks_all = zeros(num_channels, length(peak_data(i).peaks.HighBeta.(win).channel(1).max_amplitude));
        for ch = 1:num_channels
            hb_peaks_all(ch, :) = peak_data(i).peaks.HighBeta.(win).channel(ch).max_amplitude;
        end
        hb_peaks(:, w) = mean(hb_peaks_all(:, T1T2_trials), 2, 'omitnan');
    end
    for ch = 1:num_channels
        plot(timeCenters, hb_peaks(ch, :), '-o', 'Color', channel_colors(ch, :), ...
            'LineWidth', 2, 'DisplayName', chanlocs(ch).labels);
    end
    % Shade ERP windows (skip baseline)
    for w = 2:length(all_windows)
        boundsIndex = w - 1;
        tStart = windowBounds(boundsIndex, 1);
        tEnd   = windowBounds(boundsIndex, 2);
        h = fill([tStart, tEnd, tEnd, tStart], [min(hb_peaks(:))-0.1, min(hb_peaks(:))-0.1, max(hb_peaks(:))+0.1, max(hb_peaks(:))+0.1], ...
            [0.8, 0.9, 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
    end
    xticks(timeCenters);
    xticklabels(all_windows);
    xlabel('Time (ms)');
    ylabel('Peak Z-Score (μV)');
    title(['HB Peaks Across ERP Windows - ' participant_id]);
    grid on;
    legend('Location','best');
    hold off;
    saveas(gcf, ['H4_PeakGraph_HB_' participant_id '.png']);
    close(gcf);
end

%% --- Z-Score Graphs for Individuals ---
for i = 1:num_participants
    participant_id = participants{i};
    
    %% FAA Z-Score Graph with W baseline
    figure('Name', ['FAA Z-Scores - ' participant_id], 'Position', [100, 100, 800, 600]);
    hold on;
    faa_zscores = zeros(1, length(all_windows));
    % W segment
    alpha_peaks_w = zeros(num_channels, length(peak_data(i).peaks.Alpha.P200.channel(1).max_amplitude));
    for ch = 1:num_channels
        alpha_peaks_w(ch, :) = peak_data(i).peaks.Alpha.P200.channel(ch).max_amplitude;
    end
    f3_z_w = mean(alpha_peaks_w(f3_idx, W_trials), 'omitnan');
    f4_z_w = mean(alpha_peaks_w(f4_idx, W_trials), 'omitnan');
    faa_zscores(1) = f4_z_w - f3_z_w;  % FAA Z-Score = F4 - F3
    % T1T2 segments
    for w = 2:length(all_windows)
        win = all_windows{w};
        alpha_peaks = zeros(num_channels, length(peak_data(i).peaks.Alpha.(win).channel(1).max_amplitude));
        for ch = 1:num_channels
            alpha_peaks(ch, :) = peak_data(i).peaks.Alpha.(win).channel(ch).max_amplitude;
        end
        f3_z = mean(alpha_peaks(f3_idx, T1T2_trials), 'omitnan');
        f4_z = mean(alpha_peaks(f4_idx, T1T2_trials), 'omitnan');
        faa_zscores(w) = f4_z - f3_z;
    end
    plot(timeCenters, faa_zscores, '-o', 'Color', [1, 0.5, 0], 'LineWidth', 2);
    % Shade ERP windows (skip baseline)
    for w = 2:length(all_windows)
        boundsIndex = w - 1;
        tStart = windowBounds(boundsIndex, 1);
        tEnd   = windowBounds(boundsIndex, 2);
        h = fill([tStart, tEnd, tEnd, tStart], [min(faa_zscores)-0.1, min(faa_zscores)-0.1, max(faa_zscores)+0.1, max(faa_zscores)+0.1], ...
            [0.8, 0.9, 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
    end
    xticks(timeCenters);
    xticklabels(all_windows);
    xlabel('Time (ms)');
    ylabel('FAA Z-Score');
    title(['FAA Z-Score Across ERP Windows - ' participant_id]);
    grid on;
    hold off;
    saveas(gcf, ['H4_ZScoreGraph_FAA_' participant_id '.png']);
    close(gcf);
    faa_zscores_all(i, :) = faa_zscores;  % Store for stats
    
    %% BAR Z-Score Graph (average across channels) with W baseline
    figure('Name', ['BAR Z-Scores - ' participant_id], 'Position', [100, 100, 800, 600]);
    hold on;
    bar_zscores = zeros(1, length(all_windows));
    % W segment
    alpha_peaks_w = zeros(num_channels, length(peak_data(i).peaks.Alpha.P200.channel(1).max_amplitude));
    beta_peaks_w  = zeros(num_channels, length(peak_data(i).peaks.Beta.P200.channel(1).max_amplitude));
    for ch = 1:num_channels
        alpha_peaks_w(ch, :) = peak_data(i).peaks.Alpha.P200.channel(ch).max_amplitude;
        beta_peaks_w(ch, :)  = peak_data(i).peaks.Beta.P200.channel(ch).max_amplitude;
    end
    bar_per_trial_w = beta_peaks_w(:, W_trials) ./ alpha_peaks_w(:, W_trials);
    bar_T1T2_w = mean(bar_per_trial_w, 2, 'omitnan');
    bar_zscores(1) = mean(bar_T1T2_w, 'omitnan');
    % T1T2 segments
    for w = 2:length(all_windows)
        win = all_windows{w};
        alpha_peaks = zeros(num_channels, length(peak_data(i).peaks.Alpha.(win).channel(1).max_amplitude));
        beta_peaks  = zeros(num_channels, length(peak_data(i).peaks.Beta.(win).channel(1).max_amplitude));
        for ch = 1:num_channels
            alpha_peaks(ch, :) = peak_data(i).peaks.Alpha.(win).channel(ch).max_amplitude;
            beta_peaks(ch, :)  = peak_data(i).peaks.Beta.(win).channel(ch).max_amplitude;
        end
        bar_per_trial = beta_peaks(:, T1T2_trials) ./ alpha_peaks(:, T1T2_trials);
        bar_T1T2 = mean(bar_per_trial, 2, 'omitnan');
        bar_zscores(w) = mean(bar_T1T2, 'omitnan');
    end
    plot(timeCenters, bar_zscores, '-o', 'Color', [1, 0.5, 0], 'LineWidth', 2);
    % Shade ERP windows (skip baseline)
    for w = 2:length(all_windows)
        boundsIndex = w - 1;
        tStart = windowBounds(boundsIndex, 1);
        tEnd   = windowBounds(boundsIndex, 2);
        h = fill([tStart, tEnd, tEnd, tStart], [min(bar_zscores)-0.1, min(bar_zscores)-0.1, max(bar_zscores)+0.1, max(bar_zscores)+0.1], ...
            [0.8, 0.9, 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
    end
    xticks(timeCenters);
    xticklabels(all_windows);
    xlabel('Time (ms)');
    ylabel('BAR Z-Score');
    title(['BAR Z-Score Across ERP Windows - ' participant_id]);
    grid on;
    hold off;
    saveas(gcf, ['H4_ZScoreGraph_BAR_' participant_id '.png']);
    close(gcf);
    bar_zscores_all(i, :) = bar_zscores;
    
    %% TBR Z-Score Graph (average across channels) with W baseline
    figure('Name', ['TBR Z-Scores - ' participant_id], 'Position', [100, 100, 800, 600]);
    hold on;
    tbr_zscores = zeros(1, length(all_windows));
    % W segment
    theta_peaks_w = zeros(num_channels, length(peak_data(i).peaks.Theta.P200.channel(1).max_amplitude));
    beta_peaks_w  = zeros(num_channels, length(peak_data(i).peaks.Beta.P200.channel(1).max_amplitude));
    for ch = 1:num_channels
        theta_peaks_w(ch, :) = peak_data(i).peaks.Theta.P200.channel(ch).max_amplitude;
        beta_peaks_w(ch, :)  = peak_data(i).peaks.Beta.P200.channel(ch).max_amplitude;
    end
    tbr_per_trial_w = theta_peaks_w(:, W_trials) ./ beta_peaks_w(:, W_trials);
    tbr_T1T2_w = mean(tbr_per_trial_w, 2, 'omitnan');
    tbr_zscores(1) = mean(tbr_T1T2_w, 'omitnan');
    % T1T2 segments
    for w = 2:length(all_windows)
        win = all_windows{w};
        theta_peaks = zeros(num_channels, length(peak_data(i).peaks.Theta.(win).channel(1).max_amplitude));
        beta_peaks  = zeros(num_channels, length(peak_data(i).peaks.Beta.(win).channel(1).max_amplitude));
        for ch = 1:num_channels
            theta_peaks(ch, :) = peak_data(i).peaks.Theta.(win).channel(ch).max_amplitude;
            beta_peaks(ch, :)  = peak_data(i).peaks.Beta.(win).channel(ch).max_amplitude;
        end
        tbr_per_trial = theta_peaks(:, T1T2_trials) ./ beta_peaks(:, T1T2_trials);
        tbr_T1T2 = mean(tbr_per_trial, 2, 'omitnan');
        tbr_zscores(w) = mean(tbr_T1T2, 'omitnan');
    end
    plot(timeCenters, tbr_zscores, '-o', 'Color', [1, 0.5, 0], 'LineWidth', 2);
    % Shade ERP windows (skip baseline)
    for w = 2:length(all_windows)
        boundsIndex = w - 1;
        tStart = windowBounds(boundsIndex, 1);
        tEnd   = windowBounds(boundsIndex, 2);
        h = fill([tStart, tEnd, tEnd, tStart], [min(tbr_zscores)-0.1, min(tbr_zscores)-0.1, max(tbr_zscores)+0.1, max(tbr_zscores)+0.1], ...
            [0.8, 0.9, 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
    end
    xticks(timeCenters);
    xticklabels(all_windows);
    xlabel('Time (ms)');
    ylabel('TBR Z-Score');
    title(['TBR Z-Score Across ERP Windows - ' participant_id]);
    grid on;
    hold off;
    saveas(gcf, ['H4_ZScoreGraph_TBR_' participant_id '.png']);
    close(gcf);
    tbr_zscores_all(i, :) = tbr_zscores;
    
    %% HB Z-Score Graph (average across channels) with W baseline
    figure('Name', ['HB Z-Scores - ' participant_id], 'Position', [100, 100, 800, 600]);
    hold on;
    hb_zscores = zeros(1, length(all_windows));
    % W segment
    hb_peaks_w = zeros(num_channels, length(peak_data(i).peaks.HighBeta.P200.channel(1).max_amplitude));
    for ch = 1:num_channels
        hb_peaks_w(ch, :) = peak_data(i).peaks.HighBeta.P200.channel(ch).max_amplitude;
    end
    hb_T1T2_w = mean(hb_peaks_w(:, W_trials), 2, 'omitnan');
    hb_zscores(1) = mean(hb_T1T2_w, 'omitnan');
    % T1T2 segments
    for w = 2:length(all_windows)
        win = all_windows{w};
        hb_peaks_all = zeros(num_channels, length(peak_data(i).peaks.HighBeta.(win).channel(1).max_amplitude));
        for ch = 1:num_channels
            hb_peaks_all(ch, :) = peak_data(i).peaks.HighBeta.(win).channel(ch).max_amplitude;
        end
        hb_T1T2 = mean(hb_peaks_all(:, T1T2_trials), 2, 'omitnan');
        hb_zscores(w) = mean(hb_T1T2, 'omitnan');
    end
    plot(timeCenters, hb_zscores, '-o', 'Color', [1, 0.5, 0], 'LineWidth', 2);
    % Shade ERP windows (skip baseline)
    for w = 2:length(all_windows)
        boundsIndex = w - 1;
        tStart = windowBounds(boundsIndex, 1);
        tEnd   = windowBounds(boundsIndex, 2);
        h = fill([tStart, tEnd, tEnd, tStart], [min(hb_zscores)-0.1, min(hb_zscores)-0.1, max(hb_zscores)+0.1, max(hb_zscores)+0.1], ...
            [0.8, 0.9, 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
    end
    xticks(timeCenters);
    xticklabels(all_windows);
    xlabel('Time (ms)');
    ylabel('HB Z-Score');
    title(['HB Z-Score Across ERP Windows - ' participant_id]);
    grid on;
    hold off;
    saveas(gcf, ['H4_ZScoreGraph_HB_' participant_id '.png']);
    close(gcf);
    hb_zscores_all(i, :) = hb_zscores;
end

%% --- Peak Graphs for Y* Group (PT1T2YY, PT3YY, PT4YY, PT7YN) ---
% FAA Peak Graph (Y* Group)
figure('Name', 'FAA Peaks - Y* Group', 'Position', [100, 100, 800, 600]);
hold on;
f3_peaks_y = zeros(length(y_group_indices), length(all_windows));
f4_peaks_y = zeros(length(y_group_indices), length(all_windows));
for p = 1:length(y_group_indices)
    i = y_group_indices(p);
    % W segment
    alpha_peaks_w = zeros(num_channels, length(peak_data(i).peaks.Alpha.P200.channel(1).max_amplitude));
    for ch = 1:num_channels
        alpha_peaks_w(ch, :) = peak_data(i).peaks.Alpha.P200.channel(ch).max_amplitude;
    end
    f3_peaks_y(p, 1) = mean(alpha_peaks_w(f3_idx, W_trials), 'omitnan');
    f4_peaks_y(p, 1) = mean(alpha_peaks_w(f4_idx, W_trials), 'omitnan');
    % T1T2 segments
    for w = 2:length(all_windows)
        win = all_windows{w};
        alpha_peaks = zeros(num_channels, length(peak_data(i).peaks.Alpha.(win).channel(1).max_amplitude));
        for ch = 1:num_channels
            alpha_peaks(ch, :) = peak_data(i).peaks.Alpha.(win).channel(ch).max_amplitude;
        end
        f3_peaks_y(p, w) = mean(alpha_peaks(f3_idx, T1T2_trials), 'omitnan');
        f4_peaks_y(p, w) = mean(alpha_peaks(f4_idx, T1T2_trials), 'omitnan');
    end
end
f3_peaks_y_mean = mean(f3_peaks_y, 1, 'omitnan');
f4_peaks_y_mean = mean(f4_peaks_y, 1, 'omitnan');
plot(timeCenters, f3_peaks_y_mean, '-o', 'Color', [1, 0.5, 0], 'LineWidth', 2, 'DisplayName', 'F3');
plot(timeCenters, f4_peaks_y_mean, '-o', 'Color', [0, 0.5, 1], 'LineWidth', 2, 'DisplayName', 'F4');
% Shade ERP windows (skip baseline)
for w = 2:length(all_windows)
    boundsIndex = w - 1;
    tStart = windowBounds(boundsIndex, 1);
    tEnd   = windowBounds(boundsIndex, 2);
    h = fill([tStart, tEnd, tEnd, tStart], [yMin_default, yMin_default, yMax_default, yMax_default], [0.8, 0.9, 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
end
xticks(timeCenters);
xticklabels(all_windows);
xlabel('Time (ms)');
ylabel('Peak Z-Score (μV)');
title('FAA Peaks Across ERP Components (F3 vs F4) - Y* Group (PT1T2YY, PT3YY, PT4YY, PT7YN)');
xlim([0, 800]);
ylim([yMin_default, yMax_default]);
grid on;
legend('Location','best');
hold off;
saveas(gcf, 'H4_PeakGraph_FAA_YStar.png');
close(gcf);

%% BAR Peak Graph (Y* Group)
figure('Name', 'BAR Peaks - Y* Group', 'Position', [100, 100, 800, 600]);
hold on;
bar_peaks_y = zeros(length(y_group_indices), num_channels, length(all_windows));
for p = 1:length(y_group_indices)
    i = y_group_indices(p);
    % W segment
    alpha_peaks_w = zeros(num_channels, length(peak_data(i).peaks.Alpha.P200.channel(1).max_amplitude));
    beta_peaks_w  = zeros(num_channels, length(peak_data(i).peaks.Beta.P200.channel(1).max_amplitude));
    for ch = 1:num_channels
        alpha_peaks_w(ch, :) = peak_data(i).peaks.Alpha.P200.channel(ch).max_amplitude;
        beta_peaks_w(ch, :)  = peak_data(i).peaks.Beta.P200.channel(ch).max_amplitude;
    end
    bar_per_trial_w = beta_peaks_w(:, W_trials) ./ alpha_peaks_w(:, W_trials);
    bar_peaks_y(p, :, 1) = mean(bar_per_trial_w, 2, 'omitnan');
    % T1T2 segments
    for w = 2:length(all_windows)
        win = all_windows{w};
        alpha_peaks = zeros(num_channels, length(peak_data(i).peaks.Alpha.(win).channel(1).max_amplitude));
        beta_peaks  = zeros(num_channels, length(peak_data(i).peaks.Beta.(win).channel(1).max_amplitude));
        for ch = 1:num_channels
            alpha_peaks(ch, :) = peak_data(i).peaks.Alpha.(win).channel(ch).max_amplitude;
            beta_peaks(ch, :)  = peak_data(i).peaks.Beta.(win).channel(ch).max_amplitude;
        end
        bar_per_trial = beta_peaks(:, T1T2_trials) ./ alpha_peaks(:, T1T2_trials);
        bar_peaks_y(p, :, w) = mean(bar_per_trial, 2, 'omitnan');
    end
end
bar_peaks_y_mean = mean(bar_peaks_y, 1, 'omitnan');
bar_peaks_y_mean = squeeze(bar_peaks_y_mean);
for ch = 1:num_channels
    plot(timeCenters, bar_peaks_y_mean(ch, :), '-o', 'Color', channel_colors(ch, :), ...
        'LineWidth', 2, 'DisplayName', chanlocs(ch).labels);
end
% Shade ERP windows (skip baseline)
for w = 2:length(all_windows)
    boundsIndex = w - 1;
    tStart = windowBounds(boundsIndex, 1);
    tEnd   = windowBounds(boundsIndex, 2);
    h = fill([tStart, tEnd, tEnd, tStart], [min(bar_peaks_y_mean(:))-0.1, min(bar_peaks_y_mean(:))-0.1, max(bar_peaks_y_mean(:))+0.1, max(bar_peaks_y_mean(:))+0.1], ...
        [0.8, 0.9, 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
end
xticks(timeCenters);
xticklabels(all_windows);
xlabel('Time (ms)');
ylabel('Peak BAR');
title('BAR Peaks Across ERP Windows - Y* Group (PT1T2YY, PT3YY, PT4YY, PT7YN)');
grid on;
legend('Location','best');
hold off;
saveas(gcf, 'H4_PeakGraph_BAR_YStar.png');
close(gcf);

%% TBR Peak Graph (Y* Group)
figure('Name', 'TBR Peaks - Y* Group', 'Position', [100, 100, 800, 600]);
hold on;
tbr_peaks_y = zeros(length(y_group_indices), num_channels, length(all_windows));
for p = 1:length(y_group_indices)
    i = y_group_indices(p);
    % W segment
    theta_peaks_w = zeros(num_channels, length(peak_data(i).peaks.Theta.P200.channel(1).max_amplitude));
    beta_peaks_w  = zeros(num_channels, length(peak_data(i).peaks.Beta.P200.channel(1).max_amplitude));
    for ch = 1:num_channels
        theta_peaks_w(ch, :) = peak_data(i).peaks.Theta.P200.channel(ch).max_amplitude;
        beta_peaks_w(ch, :)  = peak_data(i).peaks.Beta.P200.channel(ch).max_amplitude;
    end
    tbr_per_trial_w = theta_peaks_w(:, W_trials) ./ beta_peaks_w(:, W_trials);
    tbr_peaks_y(p, :, 1) = mean(tbr_per_trial_w, 2, 'omitnan');
    % T1T2 segments
    for w = 2:length(all_windows)
        win = all_windows{w};
        theta_peaks = zeros(num_channels, length(peak_data(i).peaks.Theta.(win).channel(1).max_amplitude));
        beta_peaks  = zeros(num_channels, length(peak_data(i).peaks.Beta.(win).channel(1).max_amplitude));
        for ch = 1:num_channels
            theta_peaks(ch, :) = peak_data(i).peaks.Theta.(win).channel(ch).max_amplitude;
            beta_peaks(ch, :)  = peak_data(i).peaks.Beta.(win).channel(ch).max_amplitude;
        end
        tbr_per_trial = theta_peaks(:, T1T2_trials) ./ beta_peaks(:, T1T2_trials);
        tbr_peaks_y(p, :, w) = mean(tbr_per_trial, 2, 'omitnan');
    end
end
tbr_peaks_y_mean = mean(tbr_peaks_y, 1, 'omitnan');
tbr_peaks_y_mean = squeeze(tbr_peaks_y_mean);
for ch = 1:num_channels
    plot(timeCenters, tbr_peaks_y_mean(ch, :), '-o', 'Color', channel_colors(ch, :), ...
        'LineWidth', 2, 'DisplayName', chanlocs(ch).labels);
end
% Shade ERP windows (skip baseline)
for w = 2:length(all_windows)
    boundsIndex = w - 1;
    tStart = windowBounds(boundsIndex, 1);
    tEnd   = windowBounds(boundsIndex, 2);
    h = fill([tStart, tEnd, tEnd, tStart], [min(tbr_peaks_y_mean(:))-0.1, min(tbr_peaks_y_mean(:))-0.1, max(tbr_peaks_y_mean(:))+0.1, max(tbr_peaks_y_mean(:))+0.1], ...
        [0.8, 0.9, 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
end
xticks(timeCenters);
xticklabels(all_windows);
xlabel('Time (ms)');
ylabel('Peak TBR');
title('TBR Peaks Across ERP Windows - Y* Group (PT1T2YY, PT3YY, PT4YY, PT7YN)');
grid on;
legend('Location','best');
hold off;
saveas(gcf, 'H4_PeakGraph_TBR_YStar.png');
close(gcf);

%% HB Peak Graph (Y* Group)
figure('Name', 'HB Peaks - Y* Group', 'Position', [100, 100, 800, 600]);
hold on;
hb_peaks_y = zeros(length(y_group_indices), num_channels, length(all_windows));
for p = 1:length(y_group_indices)
    i = y_group_indices(p);
    % W segment
    hb_peaks_w = zeros(num_channels, length(peak_data(i).peaks.HighBeta.P200.channel(1).max_amplitude));
    for ch = 1:num_channels
        hb_peaks_w(ch, :) = peak_data(i).peaks.HighBeta.P200.channel(ch).max_amplitude;
    end
    hb_peaks_y(p, :, 1) = mean(hb_peaks_w(:, W_trials), 2, 'omitnan');
    % T1T2 segments
    for w = 2:length(all_windows)
        win = all_windows{w};
        hb_peaks_all = zeros(num_channels, length(peak_data(i).peaks.HighBeta.(win).channel(1).max_amplitude));
        for ch = 1:num_channels
            hb_peaks_all(ch, :) = peak_data(i).peaks.HighBeta.(win).channel(ch).max_amplitude;
        end
        hb_peaks_y(p, :, w) = mean(hb_peaks_all(:, T1T2_trials), 2, 'omitnan');
    end
end
hb_peaks_y_mean = mean(hb_peaks_y, 1, 'omitnan');
hb_peaks_y_mean = squeeze(hb_peaks_y_mean);
for ch = 1:num_channels
    plot(timeCenters, hb_peaks_y_mean(ch, :), '-o', 'Color', channel_colors(ch, :), ...
        'LineWidth', 2, 'DisplayName', chanlocs(ch).labels);
end
% Shade ERP windows (skip baseline)
for w = 2:length(all_windows)
    boundsIndex = w - 1;
    tStart = windowBounds(boundsIndex, 1);
    tEnd   = windowBounds(boundsIndex, 2);
    h = fill([tStart, tEnd, tEnd, tStart], [min(hb_peaks_y_mean(:))-0.1, min(hb_peaks_y_mean(:))-0.1, max(hb_peaks_y_mean(:))+0.1, max(hb_peaks_y_mean(:))+0.1], ...
        [0.8, 0.9, 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
end
xticks(timeCenters);
xticklabels(all_windows);
xlabel('Time (ms)');
ylabel('Peak Z-Score (μV)');
title('HB Peaks Across ERP Windows - Y* Group (PT1T2YY, PT3YY, PT4YY, PT7YN)');
grid on;
legend('Location','best');
hold off;
saveas(gcf, 'H4_PeakGraph_HB_YStar.png');
close(gcf);

%% --- Z-Score Graphs for Y* Group ---
% FAA Z-Score Graph (Y* Group)
figure('Name', 'FAA Z-Scores - Y* Group', 'Position', [100, 100, 800, 600]);
hold on;
faa_zscores_y = zeros(length(y_group_indices), length(all_windows));
for p = 1:length(y_group_indices)
    i = y_group_indices(p);
    % W segment
    alpha_peaks_w = zeros(num_channels, length(peak_data(i).peaks.Alpha.P200.channel(1).max_amplitude));
    for ch = 1:num_channels
        alpha_peaks_w(ch, :) = peak_data(i).peaks.Alpha.P200.channel(ch).max_amplitude;
    end
    f3_z_w = mean(alpha_peaks_w(f3_idx, W_trials), 'omitnan');
    f4_z_w = mean(alpha_peaks_w(f4_idx, W_trials), 'omitnan');
    faa_zscores_y(p, 1) = f4_z_w - f3_z_w;
    % T1T2 segments
    for w = 2:length(all_windows)
        win = all_windows{w};
        alpha_peaks = zeros(num_channels, length(peak_data(i).peaks.Alpha.(win).channel(1).max_amplitude));
        for ch = 1:num_channels
            alpha_peaks(ch, :) = peak_data(i).peaks.Alpha.(win).channel(ch).max_amplitude;
        end
        f3_z = mean(alpha_peaks(f3_idx, T1T2_trials), 'omitnan');
        f4_z = mean(alpha_peaks(f4_idx, T1T2_trials), 'omitnan');
        faa_zscores_y(p, w) = f4_z - f3_z;
    end
end
faa_zscores_y_mean = mean(faa_zscores_y, 1, 'omitnan');
plot(timeCenters, faa_zscores_y_mean, '-o', 'Color', [1, 0.5, 0], 'LineWidth', 2);
% Shade ERP windows (skip baseline)
for w = 2:length(all_windows)
    boundsIndex = w - 1;
    tStart = windowBounds(boundsIndex, 1);
    tEnd   = windowBounds(boundsIndex, 2);
    h = fill([tStart, tEnd, tEnd, tStart], [min(faa_zscores_y_mean)-0.1, min(faa_zscores_y_mean)-0.1, max(faa_zscores_y_mean)+0.1, max(faa_zscores_y_mean)+0.1], ...
        [0.8, 0.9, 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
end
xticks(timeCenters);
xticklabels(all_windows);
xlabel('Time (ms)');
ylabel('FAA Z-Score');
title('FAA Z-Score Across ERP Windows - Y* Group (PT1T2YY, PT3YY, PT4YY, PT7YN)');
grid on;
hold off;
saveas(gcf, 'H4_ZScoreGraph_FAA_YStar.png');
close(gcf);

%% BAR Z-Score Graph (Y* Group)
figure('Name', 'BAR Z-Scores - Y* Group', 'Position', [100, 100, 800, 600]);
hold on;
bar_zscores_y = zeros(length(y_group_indices), length(all_windows));
for p = 1:length(y_group_indices)
    i = y_group_indices(p);
    % W segment
    alpha_peaks_w = zeros(num_channels, length(peak_data(i).peaks.Alpha.P200.channel(1).max_amplitude));
    beta_peaks_w  = zeros(num_channels, length(peak_data(i).peaks.Beta.P200.channel(1).max_amplitude));
    for ch = 1:num_channels
        alpha_peaks_w(ch, :) = peak_data(i).peaks.Alpha.P200.channel(ch).max_amplitude;
        beta_peaks_w(ch, :)  = peak_data(i).peaks.Beta.P200.channel(ch).max_amplitude;
    end
    bar_per_trial_w = beta_peaks_w(:, W_trials) ./ alpha_peaks_w(:, W_trials);
    bar_T1T2_w = mean(bar_per_trial_w, 2, 'omitnan');
    bar_zscores_y(p, 1) = mean(bar_T1T2_w, 'omitnan');
    % T1T2 segments
    for w = 2:length(all_windows)
        win = all_windows{w};
        alpha_peaks = zeros(num_channels, length(peak_data(i).peaks.Alpha.(win).channel(1).max_amplitude));
        beta_peaks  = zeros(num_channels, length(peak_data(i).peaks.Beta.(win).channel(1).max_amplitude));
        for ch = 1:num_channels
            alpha_peaks(ch, :) = peak_data(i).peaks.Alpha.(win).channel(ch).max_amplitude;
            beta_peaks(ch, :)  = peak_data(i).peaks.Beta.(win).channel(ch).max_amplitude;
        end
        bar_per_trial = beta_peaks(:, T1T2_trials) ./ alpha_peaks(:, T1T2_trials);
        bar_T1T2 = mean(bar_per_trial, 2, 'omitnan');
        bar_zscores_y(p, w) = mean(bar_T1T2, 'omitnan');
    end
end
bar_zscores_y_mean = mean(bar_zscores_y, 1, 'omitnan');
plot(timeCenters, bar_zscores_y_mean, '-o', 'Color', [1, 0.5, 0], 'LineWidth', 2);
% Shade ERP windows (skip baseline)
for w = 2:length(all_windows)
    boundsIndex = w - 1;
    tStart = windowBounds(boundsIndex, 1);
    tEnd   = windowBounds(boundsIndex, 2);
    h = fill([tStart, tEnd, tEnd, tStart], [min(bar_zscores_y_mean)-0.1, min(bar_zscores_y_mean)-0.1, max(bar_zscores_y_mean)+0.1, max(bar_zscores_y_mean)+0.1], ...
        [0.8, 0.9, 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
end
xticks(timeCenters);
xticklabels(all_windows);
xlabel('Time (ms)');
ylabel('BAR Z-Score');
title('BAR Z-Score Across ERP Windows - Y* Group (PT1T2YY, PT3YY, PT4YY, PT7YN)');
grid on;
hold off;
saveas(gcf, 'H4_ZScoreGraph_BAR_YStar.png');
close(gcf);

%% TBR Z-Score Graph (Y* Group)
figure('Name', 'TBR Z-Scores - Y* Group', 'Position', [100, 100, 800, 600]);
hold on;
tbr_zscores_y = zeros(length(y_group_indices), length(all_windows));
for p = 1:length(y_group_indices)
    i = y_group_indices(p);
    % W segment
    theta_peaks_w = zeros(num_channels, length(peak_data(i).peaks.Theta.P200.channel(1).max_amplitude));
    beta_peaks_w  = zeros(num_channels, length(peak_data(i).peaks.Beta.P200.channel(1).max_amplitude));
    for ch = 1:num_channels
        theta_peaks_w(ch, :) = peak_data(i).peaks.Theta.P200.channel(ch).max_amplitude;
        beta_peaks_w(ch, :)  = peak_data(i).peaks.Beta.P200.channel(ch).max_amplitude;
    end
    tbr_per_trial_w = theta_peaks_w(:, W_trials) ./ beta_peaks_w(:, W_trials);
    temp = mean(tbr_per_trial_w, 2, 'omitnan');
    tbr_zscores_y(p, 1) = mean(temp, 'omitnan');

    %tbr_zscores_y(p, 1) = mean(tbr_per_trial_w, 'omitnan');
    % T1T2 segments
    for w = 2:length(all_windows)
        win = all_windows{w};
        theta_peaks = zeros(num_channels, length(peak_data(i).peaks.Theta.(win).channel(1).max_amplitude));
        beta_peaks  = zeros(num_channels, length(peak_data(i).peaks.Beta.(win).channel(1).max_amplitude));
        for ch = 1:num_channels
            theta_peaks(ch, :) = peak_data(i).peaks.Theta.(win).channel(ch).max_amplitude;
            beta_peaks(ch, :)  = peak_data(i).peaks.Beta.(win).channel(ch).max_amplitude;
        end
        tbr_per_trial = theta_peaks(:, T1T2_trials) ./ beta_peaks(:, T1T2_trials);
        temp = mean(tbr_per_trial, 2, 'omitnan');
        tbr_zscores_y(p, w) = mean(temp, 'omitnan');
    end
end
tbr_zscores_y_mean = mean(tbr_zscores_y, 1, 'omitnan');
plot(timeCenters, tbr_zscores_y_mean, '-o', 'Color', [1, 0.5, 0], 'LineWidth', 2);
% Shade ERP windows (skip baseline)
for w = 2:length(all_windows)
    boundsIndex = w - 1;
    tStart = windowBounds(boundsIndex, 1);
    tEnd   = windowBounds(boundsIndex, 2);
    h = fill([tStart, tEnd, tEnd, tStart], [min(tbr_zscores_y_mean)-0.1, min(tbr_zscores_y_mean)-0.1, max(tbr_zscores_y_mean)+0.1, max(tbr_zscores_y_mean)+0.1], ...
        [0.8, 0.9, 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
end
xticks(timeCenters);
xticklabels(all_windows);
xlabel('Time (ms)');
ylabel('TBR Z-Score');
title('TBR Z-Score Across ERP Windows - Y* Group (PT1T2YY, PT3YY, PT4YY, PT7YN)');
grid on;
hold off;
saveas(gcf, 'H4_ZScoreGraph_TBR_YStar.png');
close(gcf);

%% HB Z-Score Graph (Y* Group)
figure('Name', 'HB Z-Scores - Y* Group', 'Position', [100, 100, 800, 600]);
hold on;
hb_zscores_y = zeros(length(y_group_indices), length(all_windows));
for p = 1:length(y_group_indices)
    i = y_group_indices(p);
    % W segment
    hb_peaks_w = zeros(num_channels, length(peak_data(i).peaks.HighBeta.P200.channel(1).max_amplitude));
    for ch = 1:num_channels
        hb_peaks_w(ch, :) = peak_data(i).peaks.HighBeta.P200.channel(ch).max_amplitude;
    end
    temp = mean(hb_peaks_w(:, W_trials), 'omitnan');
    hb_zscores_y(p, 1) =  mean(temp, 'omitnan');
    % T1T2 segments
    for w = 2:length(all_windows)
        win = all_windows{w};
        hb_peaks_all = zeros(num_channels, length(peak_data(i).peaks.HighBeta.(win).channel(1).max_amplitude));
        for ch = 1:num_channels
            hb_peaks_all(ch, :) = peak_data(i).peaks.HighBeta.(win).channel(ch).max_amplitude;
        end
        temp = mean(hb_peaks_all(:, T1T2_trials), 'omitnan')
        hb_zscores_y(p, w) = mean(temp, 'omitnan') ;
    end
end
hb_zscores_y_mean = mean(hb_zscores_y, 1, 'omitnan');
plot(timeCenters, hb_zscores_y_mean, '-o', 'Color', [1, 0.5, 0], 'LineWidth', 2);
% Shade ERP windows (skip baseline)
for w = 2:length(all_windows)
    boundsIndex = w - 1;
    tStart = windowBounds(boundsIndex, 1);
    tEnd   = windowBounds(boundsIndex, 2);
    h = fill([tStart, tEnd, tEnd, tStart], [min(hb_zscores_y_mean)-0.1, min(hb_zscores_y_mean)-0.1, max(hb_zscores_y_mean)+0.1, max(hb_zscores_y_mean)+0.1], ...
        [0.8, 0.9, 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
end
xticks(timeCenters);
xticklabels(all_windows);
xlabel('Time (ms)');
ylabel('HB Z-Score');
title('HB Z-Score Across ERP Windows - Y* Group (PT1T2YY, PT3YY, PT4YY, PT7YN)');
grid on;
hold off;
saveas(gcf, 'H4_ZScoreGraph_HB_YStar.png');
close(gcf);

%% --- Descriptive Statistics and Spearman’s Correlation ---
% Compute mean z-scores for T1T2 (average across ERP windows P200, N200, P300, LPP) and W
T1T2_windows = 2:length(all_windows);  % P200, N200, P300, LPP
faa_zscores_T1T2 = mean(faa_zscores_all(:, T1T2_windows), 2, 'omitnan');  % Mean across T1T2 windows
faa_zscores_w  = faa_zscores_all(:, 1);  % W segment
bar_zscores_T1T2 = mean(bar_zscores_all(:, T1T2_windows), 2, 'omitnan');
bar_zscores_w  = bar_zscores_all(:, 1);
tbr_zscores_T1T2 = mean(tbr_zscores_all(:, T1T2_windows), 2, 'omitnan');
tbr_zscores_w  = tbr_zscores_all(:, 1);
hb_zscores_T1T2  = mean(hb_zscores_all(:, T1T2_windows), 2, 'omitnan');
hb_zscores_w   = hb_zscores_all(:, 1);

% Descriptive Statistics
desc_stats = table();
desc_stats.Metric = {'FAA'; 'FAA'; 'BAR'; 'BAR'; 'TBR'; 'TBR'; 'HB'; 'HB'};
desc_stats.Condition = {'T1T2'; 'W'; 'T1T2'; 'W'; 'T1T2'; 'W'; 'T1T2'; 'W'};
desc_stats.Mean = [
    mean(faa_zscores_T1T2, 'omitnan');
    mean(faa_zscores_w, 'omitnan');
    mean(bar_zscores_T1T2, 'omitnan');
    mean(bar_zscores_w, 'omitnan');
    mean(tbr_zscores_T1T2, 'omitnan');
    mean(tbr_zscores_w, 'omitnan');
    mean(hb_zscores_T1T2, 'omitnan');
    mean(hb_zscores_w, 'omitnan')
];
desc_stats.SD = [
    std(faa_zscores_T1T2, 'omitnan');
    std(faa_zscores_w, 'omitnan');
    std(bar_zscores_T1T2, 'omitnan');
    std(bar_zscores_w, 'omitnan');
    std(tbr_zscores_T1T2, 'omitnan');
    std(tbr_zscores_w, 'omitnan');
    std(hb_zscores_T1T2, 'omitnan');
    std(hb_zscores_w, 'omitnan')
];
disp('Descriptive Statistics (All Participants):');
disp(desc_stats);
writetable(desc_stats, 'H4_Descriptive_Statistics.csv');

% Spearman’s Correlation (T1T2 vs W)
[rho_faa, p_faa] = corr(faa_zscores_T1T2, faa_zscores_w, 'Type', 'Spearman', 'Rows', 'complete');
[rho_bar, p_bar] = corr(bar_zscores_T1T2, bar_zscores_w, 'Type', 'Spearman', 'Rows', 'complete');
[rho_tbr, p_tbr] = corr(tbr_zscores_T1T2, tbr_zscores_w, 'Type', 'Spearman', 'Rows', 'complete');
[rho_hb, p_hb]   = corr(hb_zscores_T1T2, hb_zscores_w, 'Type', 'Spearman', 'Rows', 'complete');

spearman_stats = table();
spearman_stats.Metric = {'FAA'; 'BAR'; 'TBR'; 'HB'};
spearman_stats.Rho = [rho_faa; rho_bar; rho_tbr; rho_hb];
spearman_stats.PValue = [p_faa; p_bar; p_tbr; p_hb];
disp('Spearman’s Correlation (T1T2 vs W, All Participants):');
disp(spearman_stats);
writetable(spearman_stats, 'H4_Spearman_Correlation.csv');

%% --- Cohen’s d for Y* Group (and PT1T2NN for discussion) ---
% Y* Group (PT1T2YY, PT3YY, PT4YY, PT7YN)
faa_zscores_y_T1T2 = mean(faa_zscores_all(y_group_indices, T1T2_windows), 2, 'omitnan');
faa_zscores_y_w  = faa_zscores_all(y_group_indices, 1);
bar_zscores_y_T1T2 = mean(bar_zscores_all(y_group_indices, T1T2_windows), 2, 'omitnan');
bar_zscores_y_w  = bar_zscores_all(y_group_indices, 1);
tbr_zscores_y_T1T2 = mean(tbr_zscores_all(y_group_indices, T1T2_windows), 2, 'omitnan');
tbr_zscores_y_w  = tbr_zscores_all(y_group_indices, 1);
hb_zscores_y_T1T2  = mean(hb_zscores_all(y_group_indices, T1T2_windows), 2, 'omitnan');
hb_zscores_y_w   = hb_zscores_all(y_group_indices, 1);

% PT1T2NN (for discussion, N=1)
faa_zscores_pT1T2_T1T2 = mean(faa_zscores_all(1, T1T2_windows), 'omitnan');
faa_zscores_pT1T2_w  = faa_zscores_all(1, 1);
bar_zscores_pT1T2_T1T2 = mean(bar_zscores_all(1, T1T2_windows), 'omitnan');
bar_zscores_pT1T2_w  = bar_zscores_all(1, 1);
tbr_zscores_pT1T2_T1T2 = mean(tbr_zscores_all(1, T1T2_windows), 'omitnan');
tbr_zscores_pT1T2_w  = tbr_zscores_all(1, 1);
hb_zscores_pT1T2_T1T2  = mean(hb_zscores_all(1, T1T2_windows), 'omitnan');
hb_zscores_pT1T2_w   = hb_zscores_all(1, 1);

% Cohen’s d for Y* Group
cohens_d_y = zeros(4, 1);
% FAA
mean_diff_faa_y = mean(faa_zscores_y_T1T2) - mean(faa_zscores_y_w);
pooled_sd_faa_y = sqrt((std(faa_zscores_y_T1T2)^2 + std(faa_zscores_y_w)^2) / 2);
cohens_d_y(1) = mean_diff_faa_y / pooled_sd_faa_y;
% BAR
mean_diff_bar_y = mean(bar_zscores_y_T1T2) - mean(bar_zscores_y_w);
pooled_sd_bar_y = sqrt((std(bar_zscores_y_T1T2)^2 + std(bar_zscores_y_w)^2) / 2);
cohens_d_y(2) = mean_diff_bar_y / pooled_sd_bar_y;
% TBR
mean_diff_tbr_y = mean(tbr_zscores_y_T1T2) - mean(tbr_zscores_y_w);
pooled_sd_tbr_y = sqrt((std(tbr_zscores_y_T1T2)^2 + std(tbr_zscores_y_w)^2) / 2);
cohens_d_y(3) = mean_diff_tbr_y / pooled_sd_tbr_y;
% HB
mean_diff_hb_y = mean(hb_zscores_y_T1T2) - mean(hb_zscores_y_w);
pooled_sd_hb_y = sqrt((std(hb_zscores_y_T1T2)^2 + std(hb_zscores_y_w)^2) / 2);
cohens_d_y(4) = mean_diff_hb_y / pooled_sd_hb_y;

% Cohen’s d for PT1T2NN (N=1, for discussion only)
cohens_d_pT1T2 = zeros(4, 1);
% FAA
mean_diff_faa_pT1T2 = faa_zscores_pT1T2_T1T2 - faa_zscores_pT1T2_w;
pooled_sd_faa_pT1T2 = pooled_sd_faa_y;  % Using Y* group SD as proxy
cohens_d_pT1T2(1) = mean_diff_faa_pT1T2 / pooled_sd_faa_pT1T2;
% BAR
mean_diff_bar_pT1T2 = bar_zscores_pT1T2_T1T2 - bar_zscores_pT1T2_w;
pooled_sd_bar_pT1T2 = pooled_sd_bar_y;
cohens_d_pT1T2(2) = mean_diff_bar_pT1T2 / pooled_sd_bar_pT1T2;
% TBR
mean_diff_tbr_pT1T2 = tbr_zscores_pT1T2_T1T2 - tbr_zscores_pT1T2_w;
pooled_sd_tbr_pT1T2 = pooled_sd_tbr_y;
cohens_d_pT1T2(3) = mean_diff_tbr_pT1T2 / pooled_sd_tbr_pT1T2;
% HB
mean_diff_hb_pT1T2 = hb_zscores_pT1T2_T1T2 - hb_zscores_pT1T2_w;
pooled_sd_hb_pT1T2 = pooled_sd_hb_y;
cohens_d_pT1T2(4) = mean_diff_hb_pT1T2 / pooled_sd_hb_pT1T2;

% Cohen’s d Table
cohens_d_table = table();
cohens_d_table.Metric = {'FAA'; 'BAR'; 'TBR'; 'HB'; 'FAA_PT1T2NN'; 'BAR_PT1T2NN'; 'TBR_PT1T2NN'; 'HB_PT1T2NN'};
cohens_d_table.Group  = [cohens_d_y; cohens_d_pT1T2];
disp('Cohens (T1T2 vs W, All Participants):');
disp(cohens_d_table);
writetable(cohens_d_table, 'H4_Cohens_Correlation.csv');





%% PCA code
% Define participants and conditions
participants = {'PT1NN', 'PT2YY', 'PT3YY', 'PT4YY', 'PT7YN'};
conditions = {'W', 'T1', 'T2', 'T1T2'};
num_participants = length(participants);
num_conditions = length(conditions);

% Initialize data matrix: 20 rows (participant-condition pairs), 4 columns (metrics)
data_matrix = zeros(num_participants * num_conditions, 4);

% Initialize labels for participants and conditions
participant_labels = repelem(participants, num_conditions);
condition_labels = repmat(conditions, 1, num_participants);

% Extract z-scored values from P300 window
for i = 1:num_participants
    for j = 1:num_conditions
        seg_name = conditions{j};
        idx = (i-1) * num_conditions + j;
        
        % Extract values from structures
        faa = faa_results(i).windows.P300.(seg_name).FAA;
        hb = hb_results(i).windows.P300.(seg_name).z_hb;
        bar = segment_bar_results(i).windows.P300.(seg_name).BAR;
        tbr = tbr_results(i).windows.N200.(seg_name).TBR;
        
        % Store in data matrix
        data_matrix(idx, :) = [faa, hb, bar, tbr];
    end
end

% Perform PCA
[~, score, ~, ~, explained] = pca(data_matrix);

% Extract PC1 and PC2
pc1 = score(:, 1);
pc2 = score(:, 2);

% Create PCA array table
pca_table = table(participant_labels', condition_labels', ...
                  data_matrix(:,1), data_matrix(:,2), data_matrix(:,3), data_matrix(:,4), ...
                  pc1, pc2, ...
                  'VariableNames', {'Participant', 'Condition', 'FAA', 'HighBeta', 'BAR', 'TBR', 'PC1', 'PC2'});

% Display the table in the MATLAB command window
disp('PCA Array Table:');
disp(pca_table);

% Save the table to an Excel file
writetable(pca_table, 'PCA_Array_Table.xlsx');

% Create scatter plot
figure;
hold on;

% Define colors as a matrix of RGB triplets
colors = [1, 0.75, 0.8;  % Pink
          0.5, 0, 0.5;   % Purple
          1, 0.5, 0;     % Orange
          1, 0, 0];      % Red

% Plot points for each condition
for j = 1:num_conditions
    cond = conditions{j};
    idx = strcmp(condition_labels, cond);
    scatter(pc1(idx), pc2(idx), 50, 'x', 'MarkerEdgeColor', colors(j, :), 'LineWidth', 1.5);
end

% Label axes with variance explained
xlabel(['PC1 (', num2str(explained(1), '%.1f'), '% Variance)'], 'FontSize', 12);
ylabel(['PC2 (', num2str(explained(2), '%.1f'), '% Variance)'], 'FontSize', 12);
title('PCA of EEG Metrics (FAA, High Beta, BAR, TBR)', 'FontSize', 14);

% Set axis limits
xlim([min(pc1) - 0.5, max(pc1) + 0.5]);
ylim([min(pc2) - 0.5, max(pc2) + 0.5]);

% Add grid and legend
grid on;
legend(conditions, 'Location', 'northeast', 'FontSize', 10);

hold off;


%% H2 - Consolidated Topoplots for High Beta (HB) and FAA using F3, F4, AF3, AF4
% Define parameters
window = 'P300';              % Time window of interest (kept as P300 per original code)
T1_trials = segments.T1;      % Trial indices for T1 condition
chanlocs = all_eeg(1).EEG.chanlocs;  % Channel locations
num_participants = length(participants);  % Number of participants
num_all_channels = length(chanlocs);  % Total number of channels in chanlocs

% Define the channels of interest: F3, F4, AF3, AF4
channels_of_interest = {'F3', 'F4', 'AF3', 'AF4'};
channel_indices = zeros(1, length(channels_of_interest));
for ch = 1:length(channels_of_interest)
    idx = find(strcmp({chanlocs.labels}, channels_of_interest{ch}));
    if isempty(idx)
        error(['Channel ' channels_of_interest{ch} ' not found in chanlocs.']);
    end
    channel_indices(ch) = idx;
end
num_channels = length(channel_indices);  % Only 4 channels: F3, F4, AF3, AF4

% Create a single figure with subplots for all participants and metrics
figure('Name', 'H2 - Topoplots for High Beta and FAA (P300, T1)', ...
       'Position', [100, 100, 800, 1200]);  % Adjusted size for 5x2 grid

% Loop over all participants to generate topoplots for HB and FAA
for i = 1:num_participants
    participant_id = participants{i};
    
    % --- High Beta (HB) Topoplot ---
    % Initialize matrix for high beta peaks [num_channels, num_trials]
    hb_peaks = zeros(num_channels, length(peak_data(i).peaks.HighBeta.(window).channel(1).max_amplitude));
    for ch = 1:num_channels
        hb_peaks(ch, :) = peak_data(i).peaks.HighBeta.(window).channel(channel_indices(ch)).max_amplitude;
    end
    % Compute mean across T1 trials for each channel
    hb_T1 = mean(hb_peaks(:, T1_trials), 2, 'omitnan');
    % Create a vector matching the length of chanlocs, with zeros elsewhere
    hb_vector = zeros(num_all_channels, 1);
    for ch = 1:num_channels
        hb_vector(channel_indices(ch)) = hb_T1(ch);
    end
    % Plot HB topoplot in the first column
    subplot(num_participants, 2, (i-1)*2 + 1);
    topoplot(hb_vector, chanlocs, 'maplimits', 'absmax', 'electrodes', 'on');
    colormap(jet);
    colorbar;
    caxis([0 1])
    title(['HB - ' participant_id], 'FontSize', 10);
    
    % --- FAA Topoplot ---
    % Initialize matrix for alpha peaks [num_channels, num_trials]
    alpha_peaks = zeros(num_channels, length(peak_data(i).peaks.Alpha.(window).channel(1).max_amplitude));
    for ch = 1:num_channels
        alpha_peaks(ch, :) = peak_data(i).peaks.Alpha.(window).channel(channel_indices(ch)).max_amplitude;
    end
    % Compute mean across T1 trials for each channel
    alpha_T1 = mean(alpha_peaks(:, T1_trials), 2, 'omitnan');
    % Compute FAA: Difference between right (F4, AF4) and left (F3, AF3) channels
    f3_idx = find(strcmp(channels_of_interest, 'F3'));
    f4_idx = find(strcmp(channels_of_interest, 'F4'));
    %af3_idx = find(strcmp(channels_of_interest, 'AF3'));
    %af4_idx = find(strcmp(channels_of_interest, 'AF4'));
    %faa_T1 = mean(alpha_T1(f4_idx), 'omitnan') - ...
    %        mean(alpha_T1(f3_idx), 'omitnan');  % Right - Left
    % Create a vector matching the length of chanlocs, with FAA values at relevant channels
    faa_vector = zeros(num_all_channels, 1);
    faa_vector(channel_indices(f3_idx)) = -mean(alpha_T1(f3_idx), 'omitnan');  % Left hemisphere negative
    faa_vector(channel_indices(f4_idx)) = mean(alpha_T1(f4_idx), 'omitnan');   % Right hemisphere positive
    %faa_vector(channel_indices(af3_idx)) = -faa_T1/2;
    %faa_vector(channel_indices(af4_idx)) = faa_T1/2;
    % Plot FAA topoplot in the second column
    subplot(num_participants, 2, (i-1)*2 + 2);
    topoplot(faa_vector, chanlocs, 'maplimits', 'absmax', 'electrodes', 'on');
    colormap(jet);
    colorbar;
    %caxis([-1 1])
    title(['FAA - ' participant_id], 'FontSize', 10);
    
end

% Add a super title for the entire figure
sgtitle('Topoplots for High Beta and FAA (P300, T1)', 'FontSize', 10);

% Set figure background color to white
set(gcf, 'Color', 'w');

% Add column labels for metrics (HB, FAA)
axes('Position', [0, 0.95, 1, 0.05], 'Visible', 'off');
text(0.25, 0.5, 'High Beta (HB)', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
text(0.75, 0.5, 'FAA', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

% Add row labels for participants
for i = 1:num_participants
    axes('Position', [0.02, 0.8 - (i-1)*0.18, 0.05, 0.05], 'Visible', 'off');
    text(0.5, 0.5, participants{i}, 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
end

% Save the figure as an image
saveas(gcf, 'H2_T1_Consolidated_Topoplots_HB_FAA.png');


%% H2 - Consolidated Topoplots for High Beta (HB) and FAA using F3, F4, AF3, AF4
% Define parameters
window = 'P300';              % Time window of interest (kept as P300 per original code)
W_trials = segments.W;      % Trial indices for W condition
chanlocs = all_eeg(1).EEG.chanlocs;  % Channel locations
num_participants = length(participants);  % Number of participants
num_all_channels = length(chanlocs);  % Total number of channels in chanlocs

% Define the channels of interest: F3, F4, AF3, AF4
channels_of_interest = {'F3', 'F4', 'AF3', 'AF4'};
channel_indices = zeros(1, length(channels_of_interest));
for ch = 1:length(channels_of_interest)
    idx = find(strcmp({chanlocs.labels}, channels_of_interest{ch}));
    if isempty(idx)
        error(['Channel ' channels_of_interest{ch} ' not found in chanlocs.']);
    end
    channel_indices(ch) = idx;
end
num_channels = length(channel_indices);  % Only 4 channels: F3, F4, AF3, AF4

% Create a single figure with subplots for all participants and metrics
figure('Name', 'H2 - Topoplots for High Beta and FAA (P300, W)', ...
       'Position', [100, 100, 800, 1200]);  % Adjusted size for 5x2 grid

% Loop over all participants to generate topoplots for HB and FAA
for i = 1:num_participants
    participant_id = participants{i};
    
    % --- High Beta (HB) Topoplot ---
    % Initialize matrix for high beta peaks [num_channels, num_trials]
    hb_peaks = zeros(num_channels, length(peak_data(i).peaks.HighBeta.(window).channel(1).max_amplitude));
    for ch = 1:num_channels
        hb_peaks(ch, :) = peak_data(i).peaks.HighBeta.(window).channel(channel_indices(ch)).max_amplitude;
    end
    % Compute mean across W trials for each channel
    %hb_W = mean(hb_peaks(:, W_trials), 2, 'omitnan');
    hb_W = hb_results(i).windows.(window).W.z_hb/4;
    % Create a vector matching the length of chanlocs, with zeros elsewhere
    hb_vector = zeros(num_all_channels, 1);
    for ch = 1:num_channels
        hb_vector(channel_indices(ch)) = hb_W;
    end
    % Plot HB topoplot in the first column
    subplot(num_participants, 2, (i-1)*2 + 1);
    topoplot(hb_vector, chanlocs, 'maplimits', 'absmax', 'electrodes', 'on');
    colormap(jet);
    colorbar;
    caxis([-0.001 0.001])
    title(['HB - ' participant_id], 'FontSize', 10);
    
    % --- FAA Topoplot ---
    % Initialize matrix for alpha peaks [num_channels, num_trials]
    alpha_peaks = zeros(num_channels, length(peak_data(i).peaks.Alpha.(window).channel(1).max_amplitude));
    for ch = 1:num_channels
        alpha_peaks(ch, :) = peak_data(i).peaks.Alpha.(window).channel(channel_indices(ch)).max_amplitude;
    end
    % Compute mean across W trials for each channel
    alpha_W = mean(alpha_peaks(:, W_trials), 2, 'omitnan');
    % Compute FAA: Difference between right (F4, AF4) and left (F3, AF3) channels
    f3_idx = find(strcmp(channels_of_interest, 'F3'));
    f4_idx = find(strcmp(channels_of_interest, 'F4'));
    %af3_idx = find(strcmp(channels_of_interest, 'AF3'));
    %af4_idx = find(strcmp(channels_of_interest, 'AF4'));
    %faa_W = mean(alpha_W(f4_idx), 'omitnan') - ...
    %        mean(alpha_W(f3_idx), 'omitnan');  % Right - Left
    % Create a vector matching the length of chanlocs, with FAA values at relevant channels
    faa_vector = zeros(num_all_channels, 1);
    faa_vector(channel_indices(f3_idx)) = faa_results(i).windows.(window).W.z_f3;%-mean(alpha_W(f3_idx), 'omitnan');  % Left hemisphere negative
    faa_vector(channel_indices(f4_idx)) = faa_results(i).windows.(window).W.z_f4;%mean(alpha_W(f4_idx), 'omitnan');   % Right hemisphere positive
    %faa_vector(channel_indices(af3_idx)) = -faa_W/2;
    %faa_vector(channel_indices(af4_idx)) = faa_W/2;
    % Plot FAA topoplot in the second column
    subplot(num_participants, 2, (i-1)*2 + 2);
    topoplot(faa_vector, chanlocs, 'maplimits', 'absmax', 'electrodes', 'on');
    colormap(jet);
    colorbar;
    caxis([-0.001 0.001])
    title(['FAA - ' participant_id], 'FontSize', 10);
    
end

% Add a super title for the entire figure
sgtitle('Topoplots for High Beta and FAA (P300, W)', 'FontSize', 10);

% Set figure background color to white
set(gcf, 'Color', 'w');

% Add column labels for metrics (HB, FAA)
axes('Position', [0, 0.95, 1, 0.05], 'Visible', 'off');
text(0.25, 0.5, 'High Beta (HB)', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
text(0.75, 0.5, 'FAA', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

% Add row labels for participants
for i = 1:num_participants
    axes('Position', [0.02, 0.8 - (i-1)*0.18, 0.05, 0.05], 'Visible', 'off');
    text(0.5, 0.5, participants{i}, 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
end

% Save the figure as an image
saveas(gcf, 'H2_W_Consolidated_Topoplots_HB_FAA.png');



%% H3 - Consolidated Topoplots for High Beta (HB) and FAA using F3, F4, AF3, AF4
% Define parameters
window = 'P300';              % Time window of interest (kept as P300 per original code)
T2_trials = segments.T2;      % Trial indices for T2 condition
chanlocs = all_eeg(1).EEG.chanlocs;  % Channel locations
num_participants = length(participants);  % Number of participants
num_all_channels = length(chanlocs);  % Total number of channels in chanlocs

% Define the channels of interest: F3, F4, AF3, AF4
channels_of_interest = {'F3', 'F4', 'AF3', 'AF4'};
channel_indices = zeros(1, length(channels_of_interest));
for ch = 1:length(channels_of_interest)
    idx = find(strcmp({chanlocs.labels}, channels_of_interest{ch}));
    if isempty(idx)
        error(['Channel ' channels_of_interest{ch} ' not found in chanlocs.']);
    end
    channel_indices(ch) = idx;
end
num_channels = length(channel_indices);  % Only 4 channels: F3, F4, AF3, AF4

% Create a single figure with subplots for all participants and metrics
figure('Name', 'H2 - Topoplots for High Beta and FAA (P300, T2)', ...
       'Position', [100, 100, 800, 1200]);  % Adjusted size for 5x2 grid

% Loop over all participants to generate topoplots for HB and FAA
for i = 1:num_participants
    participant_id = participants{i};
    
    % --- High Beta (HB) Topoplot ---
    % Initialize matrix for high beta peaks [num_channels, num_trials]
    hb_peaks = zeros(num_channels, length(peak_data(i).peaks.HighBeta.(window).channel(1).max_amplitude));
    for ch = 1:num_channels
        hb_peaks(ch, :) = peak_data(i).peaks.HighBeta.(window).channel(channel_indices(ch)).max_amplitude;
    end
    % Compute mean across T2 trials for each channel
    hb_T2 = mean(hb_peaks(:, T2_trials), 2, 'omitnan');
    % Create a vector matching the length of chanlocs, with zeros elsewhere
    hb_vector = zeros(num_all_channels, 1);
    for ch = 1:num_channels
        hb_vector(channel_indices(ch)) = hb_T2(ch);
    end
    % Plot HB topoplot in the first column
    subplot(num_participants, 2, (i-1)*2 + 1);
    topoplot(hb_vector, chanlocs, 'maplimits', 'absmax', 'electrodes', 'on');
    colormap(jet);
    colorbar;
    caxis([0 1])
    title(['HB - ' participant_id], 'FontSize', 10);
    
    % --- FAA Topoplot ---
    % Initialize matrix for alpha peaks [num_channels, num_trials]
    alpha_peaks = zeros(num_channels, length(peak_data(i).peaks.Alpha.(window).channel(1).max_amplitude));
    for ch = 1:num_channels
        alpha_peaks(ch, :) = peak_data(i).peaks.Alpha.(window).channel(channel_indices(ch)).max_amplitude;
    end
    % Compute mean across T2 trials for each channel
    alpha_T2 = mean(alpha_peaks(:, T2_trials), 2, 'omitnan');
    % Compute FAA: Difference between right (F4, AF4) and left (F3, AF3) channels
    f3_idx = find(strcmp(channels_of_interest, 'F3'));
    f4_idx = find(strcmp(channels_of_interest, 'F4'));
    %af3_idx = find(strcmp(channels_of_interest, 'AF3'));
    %af4_idx = find(strcmp(channels_of_interest, 'AF4'));
    %faa_T2 = mean(alpha_T2(f4_idx), 'omitnan') - ...
    %        mean(alpha_T2(f3_idx), 'omitnan');  % Right - Left
    % Create a vector matching the length of chanlocs, with FAA values at relevant channels
    faa_vector = zeros(num_all_channels, 1);
    faa_vector(channel_indices(f3_idx)) = -mean(alpha_T2(f3_idx), 'omitnan');  % Left hemisphere negative
    faa_vector(channel_indices(f4_idx)) = mean(alpha_T2(f4_idx), 'omitnan');   % Right hemisphere positive
    %faa_vector(channel_indices(af3_idx)) = -faa_T2/2;
    %faa_vector(channel_indices(af4_idx)) = faa_T2/2;
    % Plot FAA topoplot in the second column
    subplot(num_participants, 2, (i-1)*2 + 2);
    topoplot(faa_vector, chanlocs, 'maplimits', 'absmax', 'electrodes', 'on');
    colormap(jet);
    colorbar;
    %caxis([-1 1])
    title(['FAA - ' participant_id], 'FontSize', 10);
    
end

% Add a super title for the entire figure
sgtitle('Topoplots for High Beta and FAA (P300, T2)', 'FontSize', 10);

% Set figure background color to white
set(gcf, 'Color', 'w');

% Add column labels for metrics (HB, FAA)
axes('Position', [0, 0.95, 1, 0.05], 'Visible', 'off');
text(0.25, 0.5, 'High Beta (HB)', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
text(0.75, 0.5, 'FAA', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

% Add row labels for participants
for i = 1:num_participants
    axes('Position', [0.02, 0.8 - (i-1)*0.18, 0.05, 0.05], 'Visible', 'off');
    text(0.5, 0.5, participants{i}, 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
end

% Save the figure as an image
saveas(gcf, 'H3_T2_Consolidated_Topoplots_HB_FAA.png');




%% H4 - Consolidated Topoplots for High Beta (HB) and FAA using F3, F4, AF3, AF4
% Define parameters
window = 'P300';              % Time window of interest (kept as P300 per original code)
T1T2_trials = segments.T1T2;      % Trial indices for T1T2 condition
chanlocs = all_eeg(1).EEG.chanlocs;  % Channel locations
num_participants = length(participants);  % Number of participants
num_all_channels = length(chanlocs);  % Total number of channels in chanlocs

% Define the channels of interest: F3, F4, AF3, AF4
channels_of_interest = {'F3', 'F4', 'AF3', 'AF4'};
channel_indices = zeros(1, length(channels_of_interest));
for ch = 1:length(channels_of_interest)
    idx = find(strcmp({chanlocs.labels}, channels_of_interest{ch}));
    if isempty(idx)
        error(['Channel ' channels_of_interest{ch} ' not found in chanlocs.']);
    end
    channel_indices(ch) = idx;
end
num_channels = length(channel_indices);  % Only 4 channels: F3, F4, AF3, AF4

% Create a single figure with subplots for all participants and metrics
figure('Name', 'H2 - Topoplots for High Beta and FAA (P300, T1T2)', ...
       'Position', [100, 100, 800, 1200]);  % Adjusted size for 5x2 grid

% Loop over all participants to generate topoplots for HB and FAA
for i = 1:num_participants
    participant_id = participants{i};
    
    % --- High Beta (HB) Topoplot ---
    % Initialize matrix for high beta peaks [num_channels, num_trials]
    hb_peaks = zeros(num_channels, length(peak_data(i).peaks.HighBeta.(window).channel(1).max_amplitude));
    for ch = 1:num_channels
        hb_peaks(ch, :) = peak_data(i).peaks.HighBeta.(window).channel(channel_indices(ch)).max_amplitude;
    end
    % Compute mean across T1T2 trials for each channel
    hb_T1T2 = mean(hb_peaks(:, T1T2_trials), 2, 'omitnan');
    % Create a vector matching the length of chanlocs, with zeros elsewhere
    hb_vector = zeros(num_all_channels, 1);
    for ch = 1:num_channels
        hb_vector(channel_indices(ch)) = hb_T1T2(ch);
    end
    % Plot HB topoplot in the first column
    subplot(num_participants, 2, (i-1)*2 + 1);
    topoplot(hb_vector, chanlocs, 'maplimits', 'absmax', 'electrodes', 'on');
    colormap(jet);
    colorbar;
    %caxis([0 1])
    title(['HB - ' participant_id], 'FontSize', 10);
    
    % --- FAA Topoplot ---
    % Initialize matrix for alpha peaks [num_channels, num_trials]
    alpha_peaks = zeros(num_channels, length(peak_data(i).peaks.Alpha.(window).channel(1).max_amplitude));
    for ch = 1:num_channels
        alpha_peaks(ch, :) = peak_data(i).peaks.Alpha.(window).channel(channel_indices(ch)).max_amplitude;
    end
    % Compute mean across T1T2 trials for each channel
    alpha_T1T2 = mean(alpha_peaks(:, T1T2_trials), 2, 'omitnan');
    % Compute FAA: Difference between right (F4, AF4) and left (F3, AF3) channels
    f3_idx = find(strcmp(channels_of_interest, 'F3'));
    f4_idx = find(strcmp(channels_of_interest, 'F4'));
    %af3_idx = find(strcmp(channels_of_interest, 'AF3'));
    %af4_idx = find(strcmp(channels_of_interest, 'AF4'));
    %faa_T1T2 = mean(alpha_T1T2(f4_idx), 'omitnan') - ...
    %        mean(alpha_T1T2(f3_idx), 'omitnan');  % Right - Left
    % Create a vector matching the length of chanlocs, with FAA values at relevant channels
    faa_vector = zeros(num_all_channels, 1);
    faa_vector(channel_indices(f3_idx)) = -mean(alpha_T1T2(f3_idx), 'omitnan');  % Left hemisphere negative
    faa_vector(channel_indices(f4_idx)) = mean(alpha_T1T2(f4_idx), 'omitnan');   % Right hemisphere positive
    %faa_vector(channel_indices(af3_idx)) = -faa_T1T2/2;
    %faa_vector(channel_indices(af4_idx)) = faa_T1T2/2;
    % Plot FAA topoplot in the second column
    subplot(num_participants, 2, (i-1)*2 + 2);
    topoplot(faa_vector, chanlocs, 'maplimits', 'absmax', 'electrodes', 'on');
    colormap(jet);
    colorbar;
    %caxis([-1 1])
    title(['FAA - ' participant_id], 'FontSize', 10);
    
end

% Add a super title for the entire figure
sgtitle('Topoplots for High Beta and FAA (P300, T1T2)', 'FontSize', 10);

% Set figure background color to white
set(gcf, 'Color', 'w');

% Add column labels for metrics (HB, FAA)
axes('Position', [0, 0.95, 1, 0.05], 'Visible', 'off');
text(0.25, 0.5, 'High Beta (HB)', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
text(0.75, 0.5, 'FAA', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

% Add row labels for participants
for i = 1:num_participants
    axes('Position', [0.02, 0.8 - (i-1)*0.18, 0.05, 0.05], 'Visible', 'off');
    text(0.5, 0.5, participants{i}, 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
end

% Save the figure as an image
saveas(gcf, 'H4_T1T2_Consolidated_Topoplots_HB_FAA.png');
