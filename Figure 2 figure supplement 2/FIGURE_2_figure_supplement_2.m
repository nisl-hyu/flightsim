% Figure 2–figure supplement 2

% 2024/08/05 Angel Canelo & Anmo Kim

clear all; clc; close all;
% error_vec = [];
% inter_sacc = 0;

t0=0:0.001:5;
theta_bar = 20*pi/180 + 30*pi/180*t0;
dot_theta_bar = diff(theta_bar([1 1:end]))/(t0(2)-t0(1));

P = @(theta) 10e-12 * sin(theta);  % pos to torque e-11
V = @(theta) 1e-13 * (20 - 2 * theta^2);  % vel to torque e-12

y1(1) = 0; y2(1) = 0; y3(1) = 0; error_int(1) = 0;
h = 0.001;
I = 6e-14;  % moment of inertia
beta = 1e-12;  % drag coefficient
tau_bar = 0.05;


best_alpha = 0;
best_beta_th = 0;
min_error = inf;

% Grid search for alpha and beta_th
% alpha_values = 0.01:0.01:0.1;
% beta_th_values = 500:100:1400;
alpha_values = 0.01:0.001:0.08;
beta_th_values = 4:1:45;
hfig=figure();clf;set(gcf,'Color','w'); ii = 1;
for ii = 1:length(alpha_values)
    alpha = alpha_values(ii);
    for j = 1:length(beta_th_values)
        beta_th = beta_th_values(j);
        error_int(1) = 0; % Reset error_int for each parameter combination
        torque_start_index = 1;
        saccade_amplitudes = [];
        saccade_times = [];
        prevsac = 0;
        for i = 1:(length(t0) - 1)
            if error_int(i) * 180 / pi >= beta_th % + 0.25 * rand
                torque(i) = 1*(P(theta_bar(i) - y1(i)) + (dot_theta_bar(i) - y2(i)) * V(theta_bar(i) - y2(i)));
                if torque(i-1) == 0 % Detect the start of torque
                    saccade_amplitude = abs(y1(i) - y1(torque_start_index))/pi*180; % Calculate saccade amplitude
                    saccade_amplitudes = [saccade_amplitudes saccade_amplitude]; % Store saccade amplitude
                    [mval,idx] = max(y1(torque_start_index:i));
                    % saccade_time = t0(torque_start_index+idx) - t0(torque_start_index);
                    saccade_times = [saccade_times t0(i)-prevsac];
                    % saccade_times = [saccade_times saccade_time];
                    % saccade_vel = max(y2(torque_start_index:i))/pi*180; % Calculate saccade amplitude
                    % saccade_vels = [saccade_vels saccade_vel]; % Store saccade amplitude
                    torque_start_index = i;
                    prevsac = t0(i);
                end
                %prevsac = t0(i);
            else
                torque(i) = 0;
            end
            y1(i+1) = y1(i) + h * y2(i);
            y2(i+1) = y2(i) + h * (-beta * y2(i) + torque(i)) / I;
            %y3(i+1) = y3(i) + h * (-y3(i) + torque(i)) / tau_bar;
            error_int(i+1) = error_int(i) + h * (-error_int(i) + (theta_bar(i) - y1(i))) / alpha;
        end
        % error_opt(ii, j) = sqrt((0.45*mean(error_int * 180 / pi) - (mean(1e12*torque)))^2); % Store error angle
        error_opt(ii, j) = sqrt((mean(saccade_amplitudes)-9)^2 + (mean(saccade_times)-0.300)^2); % Store error angle
        error_iss(ii,j) = sqrt((mean(saccade_times)-0.300)^2);
        error_sa(ii,j) = sqrt((mean(saccade_amplitudes)-9)^2);
        if error_opt(ii, j) < min_error
            min_error = error_opt(ii, j);
            best_alpha = alpha;
            best_beta_th = beta_th;
            
            figure(hfig); clf(hfig)
            haxes(1)=subplot(3,1,1);
            plot(t0,y1/pi*180,'LineWidth',1.2, 'Color', 'blue');
            hold on
            plot(t0,theta_bar/pi*180,'black', 'LineWidth',1.2);
            %plot(t0,(theta_bar-y1)/pi*180,'red', 'LineWidth',1.2);
            plot(t0,error_int/pi*180,'red', 'LineWidth',1.2);
            legend('fly-bar','bar','error');
            ylabel('Heading (deg)');
            xlabel('Time (s)');
            haxes(2)=subplot(3,1,2);
            plot(t0,y2/pi*180,'LineWidth',1.2, 'Color', 'blue');
            ylabel('Angular velocity (deg/s)');
            xlabel('Time (s)');
            haxes(3)=subplot(3,1,3);
            plot(t0(1:end-1),torque,'LineWidth',1.2, 'Color', 'blue');
            ylabel('Torque');
            xlabel('Time (s)');
        end
    end
end

disp(['Best alpha: ', num2str(best_alpha)]);
disp(['Best beta: ', num2str(best_beta_th)]);
disp(['Objective error: ', num2str(min(error_opt,[],"all")),' degrees']);
% %%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%
% Plot heatmap
figure();
h = heatmap(beta_th_values, alpha_values, error_opt, 'Colormap',jet, 'CellLabelColor', 'None');
h.YDisplayData = flipud(h.YDisplayData);
xlabel('C (degree)');
ylabel('Tau (s)');
title('Heatmap for objective function');
colorbar;

% Plot heatmap
figure();
h = heatmap(beta_th_values, alpha_values, error_iss, 'Colormap',jet, 'CellLabelColor', 'None');
h.YDisplayData = flipud(h.YDisplayData);
xlabel('C (degree)');
ylabel('Tau (s)');
title('Heatmap for ISS');
colorbar;

% Plot heatmap
figure();
h = heatmap(beta_th_values, alpha_values, error_sa, 'Colormap',jet, 'CellLabelColor', 'None');
h.YDisplayData = flipud(h.YDisplayData);
xlabel('C (degree)');
ylabel('Tau (s)');
title('Heatmap for SA');
colorbar;
%clim([0, 0.5]);
%print('heatmap_figure.eps', '-depsc');

%% Simulation for different bar velocities
clear all;
t0=0:0.001:5;
theta_bar = 10*pi/180+20*pi/180*t0;
dot_theta_bar = diff(theta_bar([1 1:end]))/(t0(2)-t0(1));
speed = [20 30 40];
offs = [speed(1)-10 speed(2)-10 speed(3)-10];

P = @(theta) 10e-12 * sin(theta);  % pos to torque e-11
V = @(theta) 1e-13 * (20 - 2 * theta^2);  % vel to torque e-12

y1(1) = 0; y2(1) = 0; y3(1) = 0; error_int(1) = 0;
h = 0.001;
I = 6e-14;  % moment of inertia
beta = 1e-12;  % drag coefficient
tau_bar = 0.05;

hfig2=figure();clf;set(gcf,'Color','w'); ii = 1;

alpha = 0.044; %alpha_values(ii);
beta_th = 18; %beta_th_values(j);

for j=1:3   % loop for bar speeds
    theta_bar = offs(j)*pi/180+speed(j)*pi/180*t0;
    dot_theta_bar = diff(theta_bar([1 1:end]))/(t0(2)-t0(1));
    saccade_amplitudes = [];
    saccade_times = [];
    saccade_vels = [];
    y1(1) = 0; y2(1) = 0; y3(1) = 0; error_int(1) = 0;
    torque_start_index = 1;
    torque_start_detected = false;
    prevsac = 0;
    for i = 1:(length(t0) - 1)
        if abs(error_int(i) * 180 / pi) >= beta_th %+ 0.25 * rand
            torque(i) = 1*(P(theta_bar(i) - y1(i)) + (dot_theta_bar(i) - y2(i)) * V(theta_bar(i) - y2(i)));
            if torque(i-1) == 0 % Detect the start of torque
                saccade_amplitude = abs(y1(i) - y1(torque_start_index))/pi*180; % Calculate saccade amplitude
                saccade_amplitudes = [saccade_amplitudes saccade_amplitude]; % Store saccade amplitude
                [mval,idx] = max(y1(torque_start_index:i));
                saccade_time = t0(torque_start_index+idx) - t0(torque_start_index);
                % saccade_times = [saccade_times saccade_time];
                saccade_times = [saccade_times t0(i)-prevsac];
                saccade_vel = max(y2(torque_start_index:i))/pi*180; % Calculate saccade amplitude
                saccade_vels = [saccade_vels saccade_vel]; % Store saccade amplitude
                torque_start_index = i;
                prevsac = t0(i);
            end
            % nsacc = nsacc + 1;
        else
            torque(i) = 0;
        end
            y1(i+1) = y1(i) + h * y2(i);
            y2(i+1) = y2(i) + h * (-beta * y2(i) + torque(i)) / I;
            %y3(i+1) = y3(i) + h * (-y3(i) + torque(i)) / tau_bar;
            error_int(i+1) = error_int(i) + h * (-error_int(i) + (theta_bar(i) - y1(i))) / alpha;
    end
cell_amplitude{j} = saccade_amplitudes;
cell_time{j} = saccade_times;
cell_vel{j} = saccade_vels;
end
figure(hfig2); clf(hfig2)
haxes(1)=subplot(3,1,1);
plot(t0,y1/pi*180,'LineWidth',1.2, 'Color', 'blue');
hold on
plot(t0,theta_bar/pi*180,'black', 'LineWidth',1.2);
%plot(t0,(theta_bar-y1)/pi*180,'red', 'LineWidth',1.2);
plot(t0,error_int/pi*180,'red', 'LineWidth',1.2);
legend('fly-bar','bar','error');
ylabel('Heading (deg)');
xlabel('Time (s)');
haxes(2)=subplot(3,1,2);
plot(t0,y2/pi*180,'LineWidth',1.2, 'Color', 'blue');
ylabel('Angular velocity (deg/s)');
xlabel('Time (s)');
haxes(2)=subplot(3,1,3);
plot(t0(1:end-1),torque,'LineWidth',1.2, 'Color', 'blue');
ylabel('Torque');
xlabel('Time (s)');

whk = 3;
figure(); clf()
subplot(3, 1, 1);
boxplot(cell_amplitude{1}, 'Positions', 1, 'Whisker',whk);
hold on;
boxplot(cell_amplitude{2}, 'Positions', 2, 'Whisker',whk);
hold on;
boxplot(cell_amplitude{3}, 'Positions', 3, 'Whisker',whk); 
% errorbar(1:length(mean_amplitude), mean_amplitude, std_ampl, 'ro', 'LineWidth', 1.5);
xlabel('Bar speed');
ylabel('Amplitude (degree)'); ylim([6 13])
title('Saccade Amplitude');
set(gca, 'Box', 'off', 'TickDir', 'out');

% figure(); clf()
subplot(3, 1, 2);
boxplot(cell_time{1}, 'Positions', 1, 'Whisker',whk);
hold on;
boxplot(cell_time{2}, 'Positions', 2, 'Whisker',whk);
hold on;
boxplot(cell_time{3}, 'Positions', 3, 'Whisker',whk); 
% errorbar(1:length(mean_amplitude), mean_amplitude, std_ampl, 'ro', 'LineWidth', 1.5);
xlabel('Bar speed');
ylabel('Duration (s)'); ylim([0.25 0.45])
title('Saccade Duration');
set(gca, 'Box', 'off', 'TickDir', 'out');

% figure(); clf()
subplot(3, 1, 3);
boxplot(cell_vel{1}, 'Positions', 1, 'Whisker',whk);
hold on;
boxplot(cell_vel{2}, 'Positions', 2, 'Whisker',whk);
hold on;
boxplot(cell_vel{3}, 'Positions', 3, 'Whisker',whk); 
% errorbar(1:length(mean_amplitude), mean_amplitude, std_ampl, 'ro', 'LineWidth', 1.5);
xlabel('Bar speed');
ylabel('Peak velocity (deg/s)'); ylim([60 110])
title('Saccade Angular Velocity');
set(gca, 'Box', 'off', 'TickDir', 'out');