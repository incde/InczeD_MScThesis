close all
clear all

addpath('/home/incde/thesis/altatott/ARat_2022_11_29_221129_124249')

save_folder = '/home/incde/thesis/LFP/anaesthetized';
fs = 20000;
num_channels = 16;
num_channels_digital = 1;

% Setting the frequency range of the main brain waves
alpha_range = [8 12];
beta_range = [12 35];
gamma_range = [35 100];
delta_range = [0.5 4];
theta_range = [4 8];

%% Processing and visualising the trigger

fileinfo2 = dir('/home/incde/thesis/altatott/ARat_2022_11_29_221129_124249/digitalin.dat');
num_samples2 = fileinfo2.bytes / 2;
fid  = fopen('/home/incde/thesis/altatott/ARat_2022_11_29_221129_124249/digitalin.dat', 'r');
w = fread(fid, num_samples2, 'uint16');
fclose(fid);

time_values = (1:length(w)) / fs;

%% Plot the digital signal
figure;
plot(time_values, w);
ylim([-0.05 1.2])
xlabel('Time (s)');
ylabel('Voltage');
title('Digital Input Signal');

% Set the trigger threshold
threshold = 0.5;
% Find trigger events using findpeaks
trigger_events = find(w > threshold);
trigger_times = trigger_events / fs;

differences = diff(trigger_times);
threshold2 = 100;
positions = find(differences > threshold2) + 1;
positions = [1, positions(:)'];
single_trigger_times = trigger_times(positions);

xticks(round(single_trigger_times, 2));
xtickangle(45);  % Rotate the x-axis labels by 45 degrees
hold off;

% Save the figures to the specified folder
saveas(gcf, fullfile(save_folder, 'DigitalInputSignal.png'));
% Close the figure to avoid overlapping
close(gcf);

%% Calculating the duration of the recording

fileinfo = dir('/home/incde/thesis/altatott/ARat_2022_11_29_221129_124249/amplifier.dat');
num_samples = fileinfo.bytes / (num_channels * 2);
fid  = fopen('/home/incde/thesis/altatott/ARat_2022_11_29_221129_124249/amplifier.dat', 'r');
raw_recording = fread(fid, [num_channels, num_samples], 'int16');
fclose(fid);

v2 = raw_recording(1,:) * 0.195;
duration = length(v2) ./ fs;
duration_in_mins = duration/60;

%% Extracting the first OFF-ON-OFF period of the recording and saving the channels in sepaarte varibales and performing the filtering

index_of_second_trigger = find(time_values == single_trigger_times(2));
period_1 = time_values(1 : index_of_second_trigger-1);
indicies_period_1 = round(period_1 * fs);
num_channels = size(raw_recording, 1);

for channel = 1:size(raw_recording,1)
    
    % Extracting the current channel
    channel_data = raw_recording(channel, 1:indicies_period_1(end));
    % Save the variables separately
    variable_name = sprintf('channel_%d', channel);
    assignin('base', variable_name, channel_data);
    % Save the raw recordings data separately
    assignin('base', variable_name, channel_data);
    save(fullfile('/home/incde/thesis/LFP/anaesthetized', sprintf('channel_%d_.mat', channel)), variable_name);
    
    % Alpha bandpass filtering
    alpha_filtered_data = bandpass(channel_data, alpha_range, fs);
    % Save the results separately
    variable_name_alpha = sprintf('channel_%d_alpha', channel);
    assignin('base', variable_name_alpha, alpha_filtered_data);
    % Save the alpha-filtered data separately
    assignin('base', variable_name_alpha, alpha_filtered_data);
    save(fullfile('/home/incde/thesis/LFP/anaesthetized/alpha', sprintf('channel_%d_alpha.mat', channel)), variable_name_alpha);
    
    % Beta bandpass filtering
    beta_filtered_data = bandpass(channel_data, beta_range, fs);
    % Save the results separately
    variable_name_beta = sprintf('channel_%d_beta', channel);
    assignin('base', variable_name_beta, beta_filtered_data);
    % Save the beta-filtered data separately
    assignin('base', variable_name_beta, beta_filtered_data);
    save(fullfile('/home/incde/thesis/LFP/anaesthetized/beta', sprintf('channel_%d_beta.mat', channel)), variable_name_beta);
    
    % Gamma bandpass filtering
    gamma_filtered_data = bandpass(channel_data, gamma_range, fs);
    % Save the results separately
    variable_name_gamma = sprintf('channel_%d_gamma', channel);
    assignin('base', variable_name_gamma, gamma_filtered_data);
    % Save the gamma-filtered data separately
    assignin('base', variable_name_gamma, gamma_filtered_data);
    save(fullfile('/home/incde/thesis/LFP/anaesthetized/gamma', sprintf('channel_%d_gamma.mat', channel)), variable_name_gamma);
    
    % Theta bandpass filtering
    theta_filtered_data = bandpass(channel_data, theta_range, fs);
    % Save the results separately
    variable_name_theta = sprintf('channel_%d_theta', channel);
    assignin('base', variable_name_theta, theta_filtered_data);
    % Save the theta-filtered data separately
    assignin('base', variable_name_theta, theta_filtered_data);
    save(fullfile('/home/incde/thesis/LFP/anaesthetized/theta', sprintf('channel_%d_theta.mat', channel)), variable_name_theta);

   
end
%% Plotting channel_(1:16): raw and alpha (manually change alpha to beta, gamma or theta)

raw_folder = '/home/incde/thesis/LFP/anaesthetized';
alpha_folder = '/home/incde/thesis/LFP/anaesthetized/alpha';
save_folder_alpha = '/home/incde/thesis/LFP/anaesthetized/beta';
time_vector = single_trigger_times(1) + (0:fs-1) / fs; %to visualise only 1 second of the recording
trigger_time = single_trigger_times(1);

% Loop through each channel
for channel = 1:32
    
    % Load raw recording
    raw_filename = fullfile(raw_folder, sprintf('channel_%d.mat', channel));
    raw_data_struct = load(raw_filename);
    % Get field names for raw_data_struct
    field_names_raw = fieldnames(raw_data_struct);
    % Access the field containing raw data
    raw_data_field_name = field_names_raw{1};  % Assuming there's only one field
    raw_data = raw_data_struct.(raw_data_field_name);

    % Load alpha filtered signal
    alpha_filename = fullfile(alpha_folder, sprintf('channel_%d_alpha.mat', channel));
    alpha_filtered_data_struct = load(alpha_filename);
    % Get field names for alpha_filtered_data_struct
    field_names_alpha = fieldnames(alpha_filtered_data_struct);
    % Access the field containing alpha filtered data
    alpha_filtered_data_field_name = field_names_alpha{1};  % Assuming there's only one field
    alpha_filtered_data = alpha_filtered_data_struct.(alpha_filtered_data_field_name);
    
    start_index = round(single_trigger_times(1) * fs);

    % Plot raw and alpha filtered signals
    figure;
    plot(time_vector, raw_data(start_index:start_index+fs-1), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Raw Signal');
    hold on;
    plot(time_vector, alpha_filtered_data(start_index:start_index+fs-1), 'g-', 'LineWidth', 1.5, 'DisplayName', 'Alpha Filtered');
    xlabel('Time (seconds)');
    xlim([time_vector(1), time_vector(end)]);
    ylabel('Amplitude');
    title(['Channel ', num2str(channel), ' - Raw vs Alpha Filtered']);
    legend('show', 'Location', 'northwest');
    
% To visualise the whole recording, it is helpful to see the periods, however for 1 second it is not relevant    
%     % Plot thick black line until the first element of single_trigger_times
%     plot([time_vector(1), time_vector(end)], [min(raw_data), min(raw_data)], 'k-', 'LineWidth', 5, 'DisplayName', 'OFF Period');
% 
%     % Plot thick red line for the next 120 seconds
%     plot([trigger_time, trigger_time + 120], [min(raw_data), min(raw_data)], 'r-', 'LineWidth', 5, 'DisplayName', 'ON period');
%     
%     hold off;

    % Save the plot to the output folder
    filename = sprintf('channel_%d_plot.png', channel);
    saveas(gcf, fullfile(save_folder_alpha, filename));
    close(gcf);  % Close the figure to avoid displaying multiple figures
end

