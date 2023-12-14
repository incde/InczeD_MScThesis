close all
clear all

addpath('/home/incde/thesis/L8/L8_400mA_DC_2mON_4mOFF_230614_095342')

save_folder = '/home/incde/thesis/LFP/L8';
fs = 20000;

% Setting the frequency range of the main brain waves
alpha_range = [8 12];
beta_range = [12 35];
gamma_range = [35 100];
delta_range = [0.5 4];
theta_range = [4 8];

num_channels = 16;
num_channels_analog = 2;

%% Processing and visualising the trigger

fileinfo2 = dir('/home/incde/thesis/L8/L8_400mA_DC_2mON_4mOFF_230614_095342/analogin.dat');
num_samples_analog = fileinfo2.bytes ./ (num_channels_analog * 2);
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
    line([trigger_times(i), trigger_times(i)], ylim, 'Color', 'r', 'LineStyle', '--');
end
% 
legend('Analog Signal', 'Trigger Events');

differences = diff(trigger_events);
threshold2 = 1000;
positions = find(differences > threshold2) + 1;
positions = [1, positions];
single_trigger_times = trigger_times(positions);

xticks(round(single_trigger_times, 2));
hold off;

% Save the figures to the specified folder
saveas(gcf, fullfile(save_folder, 'AnalogInputSignal.png'));
% Close the figure to avoid overlapping
close(gcf);

%% Calculating the duration of the recording

fileinfo = dir('/home/incde/thesis/L8/L8_400mA_DC_2mON_4mOFF_230614_095342/amplifier.dat');
num_samples = fileinfo.bytes / (num_channels * 2);
fid  = fopen('/home/incde/thesis/L8/L8_400mA_DC_2mON_4mOFF_230614_095342/amplifier.dat', 'r');
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
    save(fullfile('/home/incde/thesis/LFP/L8', sprintf('channel_%d_.mat', channel)), variable_name);
    
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
    save(fullfile('/home/incde/thesis/LFP/L8', sprintf('channel_%d_theta.mat', channel)), variable_name_theta);

   
end


%% Plotting channel_(1:16): raw and alpha (manually change alpha to beta, gamma or theta)

raw_folder = '/home/incde/thesis/LFP/L8';
alpha_folder = '/home/incde/thesis/LFP/L8';
save_folder_alpha = '/home/incde/thesis/LFP/L8/alpha';
time_vector = single_trigger_times(1) + (0:fs-1) / fs; %to visualise only 1 second of the recording
trigger_time = single_trigger_times(1);

% Loop through each channel
for channel = 1:16
    
    % Load raw recording
    raw_filename = fullfile(raw_folder, sprintf('channel_%d.mat', channel));
    raw_data_struct = load(raw_filename);
    % Get field names for raw_data_struct
    field_names_raw = fieldnames(raw_data_struct);
    % Access the field containing raw data
    raw_data_field_name = field_names_raw{1};  % Assuming there's only one field
    raw_data = raw_data_struct.(raw_data_field_name);

    % Load alpha filtered signal
    alpha_filename = fullfile(alpha_folder, sprintf('channel_%d_gamma.mat', channel));
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
    plot(time_vector, alpha_filtered_data(start_index:start_index+fs-1), 'g-', 'LineWidth', 1.5, 'DisplayName', 'Gamma Filtered');
    xlabel('Time (seconds)');
    xlim([time_vector(1), time_vector(end)]);
    ylabel('Amplitude');
    title(['Channel ', num2str(channel), ' - Raw vs Gamma Filtered']);
    legend('show');
    
% To visualise the whole recording, it is helpful to see the periods, however for 1 second it is not relevant    
%     % Plot thick black line until the first element of single_trigger_times
%     plot([time_vector(1), time_vector(end)], [min(raw_data), min(raw_data)], 'k-', 'LineWidth', 5, 'DisplayName', 'OFF Period');
% 
%     % Plot thick red line for the next 120 seconds
%     plot([trigger_time, trigger_time + 120], [min(raw_data), min(raw_data)], 'r-', 'LineWidth', 5, 'DisplayName', 'ON period');
%     
%     hold off;
% 
%     % Save the plot to the output folder
%     filename = sprintf('channel_%d_plot.png', channel);
%     saveas(gcf, fullfile(save_folder_alpha, filename));
%     close(gcf);  % Close the figure to avoid displaying multiple figures
end

%% Loading the previously saved files for CWT

alpha_filename_7 = 'channel_7_alpha.mat';
beta_filename_7 = 'channel_7_beta.mat';
gamma_filename_7 = 'channel_7_gamma.mat';
theta_filename_7 = 'channel_7_theta.mat';

alpha_filename_11 = 'channel_11_alpha.mat';
beta_filename_11 = 'channel_11_beta.mat';
gamma_filename_11 = 'channel_11_gamma.mat';
theta_filename_11 = 'channel_11_theta.mat';

alpha_filename_9 = 'channel_9_alpha.mat';
beta_filename_9 = 'channel_9_beta.mat';
gamma_filename_9 = 'channel_9_gamma.mat';
theta_filename_9 = 'channel_9_theta.mat';

% Load the mat files
alpha_data_7 = load(alpha_filename_7);
beta_data_7 = load(beta_filename_7);
gamma_data_7 = load(gamma_filename_7);
theta_data_7 = load(theta_filename_7);

alpha_data_11 = load(alpha_filename_11);
beta_data_11 = load(beta_filename_11);
gamma_data_11 = load(gamma_filename_11);
theta_data_11 = load(theta_filename_11);

alpha_data_9 = load(alpha_filename_9);
beta_data_9 = load(beta_filename_9);
gamma_data_9 = load(gamma_filename_9);
theta_data_9 = load(theta_filename_9);

% Extract the signals from the structures
fieldnames_alpha_7 = fieldnames(alpha_data_7);
fn_a_7 = fieldnames_alpha_7{1};
alpha_channel_7 = alpha_data_7.(fn_a_7);

fieldnames_beta_7 = fieldnames(beta_data_7);
fn_b_7 = fieldnames_beta_7{1};
beta_channel_7 = beta_data_7.(fn_b_7);

fieldnames_gamma_7 = fieldnames(gamma_data_7);
fn_g_7 = fieldnames_gamma_7{1};
gamma_channel_7 = gamma_data_7.(fn_g_7);

fieldnames_theta_7 = fieldnames(theta_data_7);
fn_t_7 = fieldnames_theta_7{1};
theta_channel_7 = theta_data_7.(fn_t_7);


fieldnames_alpha_11 = fieldnames(alpha_data_11);
fn_a_11 = fieldnames_alpha_11{1};
alpha_channel_11 = alpha_data_11.(fn_a_11);

fieldnames_beta_11 = fieldnames(beta_data_11);
fn_b_11 = fieldnames_beta_11{1};
beta_channel_11 = beta_data_11.(fn_b_11);

fieldnames_gamma_11 = fieldnames(gamma_data_11);
fn_g_11 = fieldnames_gamma_11{1};
gamma_channel_11 = gamma_data_11.(fn_g_11);

fieldnames_theta_11 = fieldnames(theta_data_11);
fn_t_11 = fieldnames_theta_11{1};
theta_channel_11 = theta_data_11.(fn_t_11);

fieldnames_alpha_9 = fieldnames(alpha_data_9);
fn_a_9 = fieldnames_alpha_9{1};
alpha_channel_9 = alpha_data_9.(fn_a_9);

fieldnames_beta_9 = fieldnames(beta_data_9);
fn_b_9 = fieldnames_beta_9{1};
beta_channel_9 = beta_data_9.(fn_b_9);

fieldnames_gamma_9 = fieldnames(gamma_data_9);
fn_g_9 = fieldnames_gamma_9{1};
gamma_channel_9 = gamma_data_9.(fn_g_9);

fieldnames_theta_9 = fieldnames(theta_data_9);
fn_t_9 = fieldnames_theta_9{1};
theta_channel_9 = theta_data_9.(fn_t_9);
%%
% Define parameters for the CWT
%fs = 20000; % Adjust the sampling frequency based on your data
alpha_scales = 1:64;
beta_scales = 1:128;
gamma_scales = 1:256;
theta_scales = 1:32;
wavelet = 'morl'; % Adjust the wavelet based on your analysis goals
%%
% Perform CWT for each channel and each frequency band
cwt_channel_7_alpha = cwt(alpha_channel_7, alpha_scales, wavelet, 'SamplingFrequency', fs);
cwt_channel_7_beta = cwt(beta_channel_7, beta_scales, wavelet, 'SamplingFrequency', fs);
cwt_channel_7_gamma = cwt(gamma_channel_7, gamma_scales, wavelet, 'SamplingFrequency', fs);
cwt_channel_7_theta = cwt(theta_channel_7, delta_scales, wavelet, 'SamplingFrequency', fs);
%%
cwt_channel_11_alpha = cwt(alpha_channel_11, scales, wavelet, 'SamplingFrequency', fs);
cwt_channel_11_beta = cwt(beta_channel_11, scales, wavelet, 'SamplingFrequency', fs);
cwt_channel_11_gamma = cwt(gamma_channel_11, scales, wavelet, 'SamplingFrequency', fs);
cwt_channel_11_theta = cwt(theta_channel_11, scales, wavelet, 'SamplingFrequency', fs);

cwt_channel_9_alpha = cwt(alpha_channel_9, scales, wavelet, 'SamplingFrequency', fs);
cwt_channel_9_beta = cwt(beta_channel_9, scales, wavelet, 'SamplingFrequency', fs);
cwt_channel_9_gamma = cwt(gamma_channel_9, scales, wavelet, 'SamplingFrequency', fs);
cwt_channel_9_theta = cwt(theta_channel_9, scales, wavelet, 'SamplingFrequency', fs);

% Visualize the results
figure;

subplot(3, 4, 1); imagesc(abs(cwt_channel_7_alpha)); title('Channel 7 - Alpha');
subplot(3, 4, 2); imagesc(abs(cwt_channel_7_beta)); title('Channel 7 - Beta');
subplot(3, 4, 3); imagesc(abs(cwt_channel_7_gamma)); title('Channel 7 - Gamma');
subplot(3, 4, 4); imagesc(abs(cwt_channel_7_theta)); title('Channel 7 - Theta');

subplot(3, 4, 5); imagesc(abs(cwt_channel_11_alpha)); title('Channel 11 - Alpha');
subplot(3, 4, 6); imagesc(abs(cwt_channel_11_beta)); title('Channel 11 - Beta');
subplot(3, 4, 7); imagesc(abs(cwt_channel_11_gamma)); title('Channel 11 - Gamma');
subplot(3, 4, 8); imagesc(abs(cwt_channel_11_theta)); title('Channel 11 - Theta');

subplot(3, 4, 9); imagesc(abs(cwt_channel_9_alpha)); title('Channel 9 - Alpha');
subplot(3, 4, 10); imagesc(abs(cwt_channel_9_beta)); title('Channel 9 - Beta');
subplot(3, 4, 11); imagesc(abs(cwt_channel_9_gamma)); title('Channel 9 - Gamma');
subplot(3, 4, 12); imagesc(abs(cwt_channel_9_theta)); title('Channel 9 - Theta');





