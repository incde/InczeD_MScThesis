close all
clear all

addpath('/home/incde/thesis/npy-matlab')
addpath('/home/incde/thesis/L8/L8_400mA_DC_2mON_4mOFF_230614_095342')
addpath('/home/incde/thesis/L8/L8_400mA_DC_2mON_4mOFF_230614_095342/results_2')

save_folder = '/home/incde/thesis/FiringRate/L8_stim';
fs = 20000;
num_channels = 16;
num_channels_analog = 2;

%% Processing and visualising the trigger

% Based on RHD Application Note with slight modifications
fileinfo = dir('/home/incde/thesis/L8/L8_400mA_DC_2mON_4mOFF_230614_095342/analogin.dat');
num_samples_analog = fileinfo.bytes ./ (num_channels_analog * 2);
fid  = fopen('/home/incde/thesis/L8/L8_400mA_DC_2mON_4mOFF_230614_095342/analogin.dat', 'r');
w = fread(fid, [num_channels_analog, num_samples_analog], 'uint16');
fclose(fid);
analog_data = w(2,:) * 0.000050354;
time_values = (1:length(analog_data)) / fs;

% Plot the analog signal
figure;
plot(time_values, analog_data);
ylim([-0.05 3.5])
xlabel('Time (s)');
ylabel('Voltage');
title('Analog Input Signal');

% Set the trigger threshold
threshold = 3;
% Find trigger events using findpeaks
trigger_events = find(analog_data > threshold);
trigger_times = trigger_events / fs;
 
% Plot vetical lines at trigger event locations
hold on;
for i = 1:length(trigger_events)
    line([trigger_times(i), trigger_times(i)], ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5);
end
% 
legend('Analog Signal', 'Trigger Events');

% Saving the occurence of the trigger in s
differences = diff(trigger_events);
threshold2 = 1000;
positions = find(differences > threshold2) + 1;
positions = [1, positions];
single_trigger_times = trigger_times(positions);

% Marking the trigger times on the x axis
xticks(round(single_trigger_times, 2));
hold off;

% Save the figures to the specified folder
saveas(gcf, fullfile(save_folder, 'AnalogInputSignal.png'));
% Close the figure to avoid overlapping
close(gcf);

%% Reading in the needed files and selecting only MUA

cluster_info = readtable('/home/incde/thesis/L8/L8_400mA_DC_2mON_4mOFF_230614_095342/results_2/cluster_info.tsv', 'FileType', 'text', 'Delimiter', '\t');
spike_times = readNPY('/home/incde/thesis/L8/L8_400mA_DC_2mON_4mOFF_230614_095342/results_2/spike_times.npy');
spike_clusters = readNPY('/home/incde/thesis/L8/L8_400mA_DC_2mON_4mOFF_230614_095342/results_2/spike_clusters.npy');
amplitudes = readNPY('/home/incde/thesis/L8/L8_400mA_DC_2mON_4mOFF_230614_095342/results_2/amplitudes.npy');
all_clusters = table(spike_clusters, spike_times, amplitudes,'VariableNames', {'cluster_id', 'spike_times', 'amplitude'});

mua_cluster_ids = cluster_info.cluster_id(cluster_info.group == "mua");
noise_cluster_ids = cluster_info.cluster_id(cluster_info.group == "noise");
keep_rows = ~ismember(all_clusters.cluster_id, noise_cluster_ids);
MUA = all_clusters(keep_rows, :);

%% Creating the struct for the clusters to handle them individually bbut store them in one variable

unique_ids = unique(MUA.cluster_id);
cluster_data = struct();

for i = 1:length(unique_ids)
    cluster_id = unique_ids(i);
    current_cluster = MUA(MUA.cluster_id == cluster_id, :);
    cluster_data.(['cluster_' num2str(cluster_id)]) = current_cluster;
end

% Deleting the first column of each field
fields = fieldnames(cluster_data);
for i = 1:numel(fields)
    current_data = cluster_data.(fields{i});
    cluster_data.(fields{i}) = current_data(:, 2:end);
end

%% Calculating the duration of the recording

fileinfo2 = dir('/home/incde/thesis/L8/L8_400mA_DC_2mON_4mOFF_230614_095342/amplifier.dat');
num_samples = fileinfo2.bytes / (num_channels * 2);
fid  = fopen('/home/incde/thesis/L8/L8_400mA_DC_2mON_4mOFF_230614_095342/amplifier.dat', 'r');
v = fread(fid, [num_channels, num_samples], 'int16');
fclose(fid);

v = v(1,:) * 0.195;
duration = length(v) ./ fs;
duration_in_mins = duration/60;

%% Calculating the changes in firing rate over time in 10s bin windows

bin_size = 10;
num_bins = floor(duration / bin_size);

fields = fieldnames(cluster_data);
% Creating empty cells to store results
all_firing_rates = cell(1, numel(fields));
all_time_bins = cell(1, numel(fields));
all_firing_rates_granular = cell(1, numel(fields));

trigger_bins = round(single_trigger_times ./ bin_size);

% Iterate through each cluster field and visualising the changes in firing rate with respect to the ON and OFF periods

for i = 1:numel(fields)
    field_name = fields{i};
    
    % Extract spike times from the current cluster field
    spike_times = cluster_data.(field_name){:,1} ./ fs;
    
    % Initialise arrays to store firing rates and time bins
    firing_rates = zeros(1, num_bins);
    time_bins = zeros(1, num_bins);
    
    % Calculate firing rate for each bin
    for j = 1:num_bins
        bin_start = (j - 1) * bin_size;
        bin_end = j * bin_size;
        
        % Extract spikes within the current bin
        spikes_in_bin = spike_times((spike_times >= bin_start) &(spike_times < bin_end));
        
        % Calculate firing rate for the bin
        firing_rates(j) = length(spikes_in_bin) / bin_size;
        time_bins(j) = bin_start;
    end
    
    % Store firing rates and time bins for each field
    all_firing_rates{i} = firing_rates;
    all_time_bins{i} = time_bins;
    
    % Replicate ech firing rate value 10 times
    all_firing_rates_granular{i} = repelem(firing_rates, bin_size);
    
    % Plot the firing rates over time using a bar plot for each filed
    figure;
    bar(all_time_bins{i}, all_firing_rates{i}, 'barwidth', 0.8);
    xlabel('Time(s)');
    ylabel('Firing Rate Changes (Hz)');
    title(['Changes in Firing Rate Over Time for Cluster ' field_name(end:end) ]);
    
    % Plot ON and OFF periods at the bottom
    hold on;
    line([0, duration], [0, 0], 'Color', 'k', 'LineWidth', 7);
    
    for k = 1:length(trigger_bins)
        ON_start = (trigger_bins(k) - 1) * bin_size;
        ON_end = ON_start + 12 * bin_size;
        line([ON_start, ON_end], [0, 0], 'Color', 'r', 'LineWidth', 7);
        
    end
    
    hold off;
    
    % Save the figures to the specified folder
    saveas(gcf, fullfile(save_folder, ['CorrectedDuration_' field_name '.png']));
    % Close the figure to avoid overlapping
    close(gcf);
    
end

%% FURTHER CALCULATIONS

%% Average firing rates in ON +-12 bin periods

% Duration of each ON period in bins
stimulation_duration_bins = 12;
pre_post_bins = 12;

% Initialize a cell array to store the average firing rates for each cluster
all_average_ON_firing_rates = cell(1, numel(all_firing_rates));

% Iterate through each cluster
for cluster_idx = 1:numel(all_firing_rates)
    % Extract firing rates for the current cluster
    current_firing_rates = all_firing_rates{cluster_idx};
    
    % Initialize a vector to store the average firing rates for each stimulation period
    average_firing_rates = zeros(1, numel(trigger_bins));
    
    % Iterate through each ON +-12 bins period
    for i = 1:numel(trigger_bins)
        % Define the indices for the current stimulation period for the current cluster
        period_indices = (trigger_bins(i) - pre_post_bins):(trigger_bins(i) + stimulation_duration_bins + pre_post_bins - 1);
        
        % Extract firing rates in the current stimulation period for the current cluster
        stimulation_firing_rates = current_firing_rates(period_indices);
        
        % Calculate the average firing rate for the current stimulation period
        average_firing_rates(i) = mean(stimulation_firing_rates);
    end
    
    % Store the average firing rates for the current cluster
    all_average_ON_firing_rates{cluster_idx} = average_firing_rates;
end


%% Baseline firing rates (2 min OFF before ON)

BL_period_bins = 12;

all_BL_firing_rates = cell(1, numel(all_firing_rates));

for cluster_idx = 1:numel(all_firing_rates)
    
    current_BL_rate = all_firing_rates{cluster_idx};
    average_BL_firing_rates = zeros(1, numel(trigger_bins));
    
    for i = 1:numel(trigger_bins)
        BL_indices = (trigger_bins(i) - BL_period_bins):(trigger_bins(i) - 1);
        
        BL_firing_rates = current_BL_rate(BL_indices);
        average_BL_firing_rates(i) = mean(BL_firing_rates);
    end
    
    all_BL_firing_rates{cluster_idx} = average_BL_firing_rates;
end


%% Calculate relative change in firing rates for each cluster

all_relative_change_firing_rates = cell(1, numel(all_BL_firing_rates));

for cluster_idx = 1:numel(all_BL_firing_rates)
    % Extract firing rates for the current cluster
    BL_firing_rates = all_BL_firing_rates{cluster_idx};
    relative_firing_rates = all_average_ON_firing_rates{cluster_idx};
    
    % Calculate relative change in firing rates
    all_relative_change_firing_rates{cluster_idx} = ((relative_firing_rates - BL_firing_rates) ./ BL_firing_rates) .* 100;
end 

%% Average of relative changes in firing rates for all clusters

average_relative_change = cell(size(all_relative_change_firing_rates));
std_relative_change = cell(size(all_relative_change_firing_rates));

for cluster_idx = 1:numel(all_relative_change_firing_rates)
   
    average_relative_change{cluster_idx} = mean(all_relative_change_firing_rates{cluster_idx});
    std_relative_change{cluster_idx} = std(all_relative_change_firing_rates{cluster_idx});
end

%% The baseline increases over the stimulation

average_BL_firing_rates = [16.96 17.49 17.11 17.15 25.38]; 

figure;
for cluster_idx = 1:numel(all_BL_firing_rates)
   plot(all_BL_firing_rates{cluster_idx}, '-o', 'LineWidth', 2.5);
   hold on;
    
end

plot(average_BL_firing_rates, '-ok', 'LineWidth', 2.5);

xlabel('Number of OFF period');
ylabel('Baseline Firing Rate during OFF period');
title('Changes in Baseline Firing Rate for Each Cluster');
legend('Cluster 2', 'Cluster 5', 'Cluster 6', 'Cluster 7', 'Cluster 10', 'Cluster 12', 'Cluster 14', 'Cluster 15', 'Cluster 16', 'Average', 'Location', 'northwest');
xticks(1:numel(all_BL_firing_rates{1}));
grid on;
hold off;



