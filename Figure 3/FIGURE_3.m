% FIGURE 3

% 2023/07/24 Angel Canelo & Anmo Kim

clear all; clc; close all;

PAT_NAMES = {'uniform_gs4','grating_for_magno_short',...
    'bar_rw_for_magno_short','bar_lw_for_magno_short',...
    'grtng_rw_for_magno_short','grtng_lw_for_magno_short',...
    'spot_rw_for_magno_short','spot_lw_for_magno_short',...
    'loom_right_for_magno_short','loom_left_for_magno_short',...
    'loom_center_for_magno_short','3vBars_CW','3vBars_CCW'};

PAT_NAMES2 = {'uniform_gs4','grating_for_magno_short',...
    'bar_RWLW_combined','grtng__RWLW_combined',...
    'spot_RWLW_combined','loom_RWLW_combined',...
    'loom_center_for_magno_short','3vBars_CWCCW_combined'};

downsampleRatio = 10;
PATID = 3:13;
STIMDUR = ones(size(PATID))*7.5;    %stimulus duration for each pattern [sec]
STIMDUR((end-1):end)=12;    %12sec for the last pattern (3vBars)

%step1: load fmf file results (timestamp + body_angle) and EDR file
%datapath = '/Volumes/magnet_exp';

fig = figure(); clf;
mat_file_dir = dir('movie202304*_2.mat');
edr_file_dir = dir('2304*.EDR');

peak_ampli2{1} = [];
latency_onset2{1} = [];
num_trials = [];
for kk=1:length(mat_file_dir)
    load([mat_file_dir(kk).name], 'body_angle');
    body_angle_t = body_angle;
    load(['time_' mat_file_dir(kk).name(1:(end-6)) '.fmf.mat'], 'timestamp');
    timestamp = timestamp - timestamp(1);
    [data, h] = import_edr(edr_file_dir(kk).name);
    disp(['Opening ' mat_file_dir(kk).name(1:(end-6))])
    
    
    %% step2: combine the data files according to the timestamp
    t = data(:,1); 
    t = linspace(t(1), t(end), length(t));
    dt = t(2) - t(1);
    xpos = data(:,2)/1000;
    ypos = data(:,3)/1000;
    patid = round(2*1.03*(data(:,4)/1000-0.6));
    camtrig = data(:,5)/1000;
    
    %extract the camtrig
    trig_idx_EDR = find(camtrig>4 & camtrig([1 1:(end-1)])<=4); % get trigger beyond threshold
    % figure();plot(camtrig)
    % figure();plot(diff(trig_idx_EDR)*dt)
    
    % Get the first trigger of each EDR and fmf. Experiment start mark
    sync1_idx_fmf = find(diff(timestamp)>1.5,1,'first')+1;
    %%%% DIFFERENT CASES TO SYNC EDR %%%%%%%
    triggg = find(diff(trig_idx_EDR*dt)>1.5);
    howmany = length(triggg);   % how many triggers
    if howmany==1
        sync1_idx_EDR = find(diff(trig_idx_EDR*dt)>1.5,1,'first')+1;
        len = min(length(trig_idx_EDR) - sync1_idx_EDR, length(timestamp) - sync1_idx_fmf);
    elseif howmany==2 && triggg(2)<=20000
        sync1_idx_EDR = find(diff(trig_idx_EDR*dt)>1.5,1,'last')+1;
        len = min(length(trig_idx_EDR) - sync1_idx_EDR, length(timestamp) - sync1_idx_fmf);
    elseif howmany==2 && triggg(2)>=60000
        sync1_idx_EDR = find(diff(trig_idx_EDR*dt)>1.5,1,'first')+1;
        sync2_idx_EDR = find(diff(trig_idx_EDR*dt)>1.5,1,'last')+1;
        len = min(sync2_idx_EDR - sync1_idx_EDR, length(timestamp) - sync1_idx_fmf);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
    % t is the global time vector of EDR. dt*trig_idx_EDR is the time vector associated to the trigger
    body_angle_edr = interp1(dt*trig_idx_EDR(sync1_idx_EDR+(0:(len-1))),body_angle_t(sync1_idx_fmf+(0:(len-1))), t);  % Interpolate body angle to global EDR time vector from trigger time vector
    
    %downsample the traces
    t = t(1:downsampleRatio:end);
    dt = dt*downsampleRatio;
    xpos = xpos(1:downsampleRatio:end);
    ypos = ypos(1:downsampleRatio:end);
    patid = patid(1:downsampleRatio:end);
    body_angle_edr = body_angle_edr(1:downsampleRatio:end);
    
    
    %% step3: extract the stimulus-triggered averages
    %detect onsets & collect stim-triggered averages
    stim_trig_tr = cell(length(PATID), 3);  % Cell array to store traces for each patid
    for i = 1:length(PATID)
        segment_idx = round(-1/dt):round(STIMDUR(i)/dt+1/dt);   % block duration
        segment_idx_vec{i} = round(-1/dt):round(STIMDUR(i)/dt+1/dt);
    
        onset_idx = find(patid([1 1:(end-1)])<PATID(i) & patid==PATID(i) & patid([2:end end])==PATID(i) );
        offset_idx = find(patid([1 1 1:(end-2)])==PATID(i) & patid([1 1:(end-1)])==PATID(i) & patid<PATID(i))-1;
    
        %figure(1);clf;plot(t,patid);hold on;plot(t(onset_idx),patid(onset_idx),'ro');
        %plot(t(offset_idx),patid(offset_idx),'md');
        
        % Discarding blocks upon condition
        isValid = true(size(onset_idx));
        for j = 1:length(onset_idx)
            k = find(offset_idx>onset_idx(j),1,'first');
            if(isempty(k))
                isValid(j) = false;
            end
    
            if(abs((offset_idx(k)-onset_idx(j))*dt - STIMDUR(i))>0.2)
                isValid(j) = false;
            end
            if isnan(body_angle_edr(onset_idx(j))) % discard blocks with big sudden changes
                isValid(j) = false;
            end
        end
        onset_idx = onset_idx(isValid);
        
        baseline_idx = find(segment_idx*dt>5.2 & segment_idx*dt<5.5);   % starting from 5.2 sec
        after_onset = find(segment_idx*dt>=5.5);
        after_onset_mark = find(segment_idx*dt>=5.5, 1, 'first');
        after_onset_mark2 = find(segment_idx*dt>=5.6, 1, 'first');
        onset_segment = find(segment_idx*dt>=5.5 & segment_idx*dt<=5.7);
        response_segment = find(segment_idx*dt>=5.6 & segment_idx*dt<=5.8);
        boundary_segment = find(segment_idx*dt>=5.6 & segment_idx*dt<=8);
        for j=1:length(onset_idx)
            stim_trig_tr{i,1}(j,:) = xpos(onset_idx(j)+segment_idx);
            stim_trig_tr{i,2}(j,:) = patid(onset_idx(j)+segment_idx);
            body_angle_mean_before_stim = mean(body_angle_edr(onset_idx(j)+segment_idx(baseline_idx)));
            
            good_bad = body_angle_edr(onset_idx(j)+segment_idx)-360*round(body_angle_mean_before_stim/360);
            bad_trace{i}(j,:) = zeros(length(good_bad),1);
            
            %%%% Xpos onset latency %%%%%
            xpos_onset_segment = xpos(onset_idx(j)+segment_idx);
            xpos_onset_segment = xpos_onset_segment(onset_segment);
            threshold = mean(xpos_onset_segment); %min(xpos_onset_segment) + (max(xpos_onset_segment)-min(xpos_onset_segment))/2;
            % Detect the onsets Stim
            onsets = find(xpos_onset_segment >= threshold,1,'first');
            %plot(xpos_onset_segment)
            t_onset = dt*segment_idx(after_onset_mark+onsets); 
            t_onset_arr(kk,i) = dt*segment_idx(after_onset_mark); %dt*segment_idx(after_onset_mark+onsets);    
            %%%%%%%%%
            peak_ampli{kk}(i,j) = 0;
            latency_onset{kk}(i,j) = 0;
            %%%% CHECKING LOCK-IN AND SUDDEN CHANGE CONDITION %%%%%%%%
            if mod(i, 2) ~= 0 %&& mean(good_bad(baseline_idx)) >= -115 && mean(good_bad(baseline_idx)) <= -80 %&& length(find(abs(diff(good_bad))>1.5))==0 && min(good_bad(segment_idx(boundary_segment)))>-150 && max(good_bad(segment_idx(boundary_segment)))<150
                if i==11
                    stim_trig_tr{i-1,3} = [stim_trig_tr{i-1,3}; good_bad];
                elseif i==7
                    stim_trig_tr{i,3}(j,:) = good_bad - mean(good_bad(baseline_idx));
                    body_angle_onset_segment = good_bad(response_segment)-mean(good_bad(baseline_idx));
                    peak_ampli{kk}(i,j) = mean(body_angle_onset_segment);

                    val50 = mean(body_angle_onset_segment);
                    ind_b = find(body_angle_onset_segment<=val50,1,'first');
                    t50stim = dt*segment_idx(after_onset_mark2+ind_b);
                    if val50>0
                        latency_onset{kk}(i,j) = abs(t50stim - t_onset);
                    end
                else
                    stim_trig_tr{i,3}(j,:) = good_bad - mean(good_bad(baseline_idx));
                    body_angle_onset_segment = good_bad(response_segment)-mean(good_bad(baseline_idx));
                    peak_ampli{kk}(i,j) = mean(body_angle_onset_segment);

                    val50 = mean(body_angle_onset_segment);
                    ind_b = find(body_angle_onset_segment>=val50,1,'first');
                    t50stim = dt*segment_idx(after_onset_mark2+ind_b);
                    if val50>0
                        latency_onset{kk}(i,j) = abs(t50stim - t_onset);
                    end

                end
            elseif mod(i, 2) == 0 %&& mean(good_bad(baseline_idx)) >= 80 && mean(good_bad(baseline_idx)) <= 115 %&& length(find(abs(diff(good_bad))>1.5))==0 && min(good_bad(segment_idx(boundary_segment)))>-150 && max(good_bad(segment_idx(boundary_segment)))<150
                if i==8
                    stim_trig_tr{i-1,3} = [stim_trig_tr{i-1,3}; good_bad - mean(good_bad(baseline_idx))];
                    body_angle_onset_segment = good_bad(response_segment)-mean(good_bad(baseline_idx));
                    peak_ampli{kk}(i-1,end+1) = mean(body_angle_onset_segment);

                    val50 = mean(body_angle_onset_segment);
                    ind_b = find(body_angle_onset_segment<=val50,1,'first');
                    t50stim = dt*segment_idx(after_onset_mark2+ind_b);
                    if val50>0
                        latency_onset{kk}(i-1,end+1) = abs(t50stim - t_onset);
                    end

                elseif i==10
                    stim_trig_tr{i,3}(j,:) = -good_bad;
                else
                    stim_trig_tr{i-1,3} = [stim_trig_tr{i-1,3}; -good_bad + mean(good_bad(segment_idx(baseline_idx)))];
                    body_angle_onset_segment = - good_bad(response_segment)+mean(good_bad(baseline_idx));
                    peak_ampli{kk}(i-1,end+1) = mean(body_angle_onset_segment);

                    val50 = mean(body_angle_onset_segment);
                    ind_b = find(body_angle_onset_segment>=val50,1,'first');
                    t50stim = dt*segment_idx(after_onset_mark2+ind_b);
                    if val50>0
                        latency_onset{kk}(i-1,end+1) = abs(t50stim - t_onset);
                    end

                end
            elseif i==11 %&& mean(good_bad) >= -150 && mean(good_bad) <= -20 && length(find(abs(diff(good_bad))>1))==0
                 stim_trig_tr{i,3}(j,:) = good_bad;
            elseif i==10 %&& mean(good_bad) >= 20 && mean(good_bad) <= 150 && length(find(abs(diff(good_bad))>1))==0
                 stim_trig_tr{i,3}(j,:) = -good_bad;
            else
                bad_trace{i}(j,:) = good_bad;
                continue
            end
        end
        % Mean of body angle
    end
    ii = 1;
    iii = 1;
    body_angle_org{kk,ii} = [];
    body_angle_mean{kk,ii} = [];
    for i = 1:length(PATID)-1
        segment_idx = round(-1/dt):round(STIMDUR(i)/dt+1/dt);
        response_segment = find(segment_idx*dt>=5.5 & segment_idx*dt<=5.8);
        after_onset_mark2 = find(segment_idx*dt>=5.5, 1, 'first');
        baseline_idx = find(segment_idx*dt>5.5 & segment_idx*dt<5.6);
        response_condit = find(segment_idx*dt>=5.7 & segment_idx*dt<=5.8);
        
        stim_trig_tr{i,3}(stim_trig_tr{i,3}==0) = nan;
        check_mean = mean(stim_trig_tr{i,3},'omitnan');
        if isempty(stim_trig_tr{i,3})
            num_trials(kk,i) = 0;
        else
            num_trials(kk,i) = sum(~isnan(stim_trig_tr{i,3}(:,1)));
      
            if mod(i, 2) ~= 0 && sum(~isnan(stim_trig_tr{i,3}(:,1)))>=4 && max(check_mean(baseline_idx))<=16 && min(check_mean(baseline_idx))>=-10
                body_angle_org{kk,ii} = stim_trig_tr{i,3};
                body_angle_mean{kk,ii} = mean(stim_trig_tr{i,3},'omitnan');
                flies_valid(kk,i) = sum(~isnan(stim_trig_tr{i,3}(:,1)));
                if mean(body_angle_mean{kk,ii}(response_condit)) < 0
                    peak_ampli_mean{kk}(ii) = nan;
                    val_50 = 0.5*min(body_angle_mean{kk,ii}(response_segment));
                    ind_50 = find(body_angle_mean{kk,ii}(response_segment)<=val_50,1,'first');
                    t50onset = dt*segment_idx(after_onset_mark2+ind_50);
                    latency_onset_mean{kk}(ii) = nan;%abs(t50onset-t_onset_arr(kk,i));
                else
                    peak_ampli_mean{kk}(ii) = max(body_angle_mean{kk,ii}(response_segment));
                    val_50 = 0.5*max(body_angle_mean{kk,ii}(response_segment));
                    ind_50 = find(body_angle_mean{kk,ii}(response_segment)>=val_50,1,'first');
                    t50onset = dt*segment_idx(after_onset_mark2+ind_50);
                    latency_onset_mean{kk}(ii) = abs(t50onset-t_onset_arr(kk,i));    % 2.5 ms delay due to camera open and patid voltage
                    if latency_onset_mean{kk}(ii)<0.05
                        peak_ampli_mean{kk}(ii) = nan;
                        latency_onset_mean{kk}(ii) = nan;
                    end
                end
                ii = ii + 1;
            end
        end
    end
    %% step4: plot the results    
    downsamplePlot = 100;
    if edr_file_dir(kk).name == '230413_008.EDR'
    figure(fig);clf;set(gcf,'color','w','DefaultAxesColorOrder',gray(10));
    haxes=[];
    for i=1:length(PATID)-8
        if i==6 || i==7
            segment_idx = round(-1/dt):round(STIMDUR(i+4)/dt+1/dt);
        else
            segment_idx = round(-1/dt):round(STIMDUR(i)/dt+1/dt);
        end
        onset_segment = find(segment_idx*dt>=5.0 & segment_idx*dt<=6.5);
        for j = 1:1
            haxes(j,i) = subplot(1, length(PATID)-8, i+(j-1)*(length(PATID)-8));            
            if(j==1)
                if size(body_angle_mean,2)<=i || isempty(body_angle_org{kk,i})
                    continue
                else
                    plot(dt*segment_idx(onset_segment), body_angle_org{kk,i}(:,onset_segment));hold on;
                    plot(dt*segment_idx(onset_segment), body_angle_mean{kk,i}(:,onset_segment),'r','LineWidth',1.5)
                end
                axis tight;
                plot([5.5 5.5],[-50 50],'black');
                plot([5 6.5],[0 0],'black');
                title(PAT_NAMES2{PATID(i)},'FontSize',7,'FontWeight','normal','Interpreter','none');
            end
        end
    end
    set(haxes(1,:),'YLim',[-50 50]);
    set(haxes(:,2:end),'YTickLabel',[]);
    set(haxes(1,:),'XTickLabel',[]);
    set(haxes,'box','off','tickdir','out');
    axes('position',[0 0.95 1 0.001]);set(gca,'XTick',[],'YTick',[],'XColor','w','YColor','w');
    title(['Single fly ' edr_file_dir(kk).name],'Interpreter','none');
    end
end
count_val = 0;
for jjj=1:size(flies_valid,1)
    if sum(flies_valid(jjj,1))>=4
        count_val = count_val+1;
    end
end
%%%%%%%%%%%%%% POPULATION %%%%%%%%%%%%%%%%
% Traces per pattern ID
    downsamplePlot = 100;
    figure();clf;set(gcf,'DefaultAxesColorOrder',gray(10));
    for i=1:length(PATID)-8
        if i==6 || i==7
            segment_idx = round(-1/dt):round(STIMDUR(i+4)/dt+1/dt);
        else
            segment_idx = round(-1/dt):round(STIMDUR(i)/dt+1/dt);
        end
        onset_segment = find(segment_idx*dt>=5.0 & segment_idx*dt<=6.5);        
        for j = 1:1
            haxes2(j,i) = subplot(1, length(PATID)-8, i+(j-1)*(length(PATID)-8));       
                nn = 1;
                for kk=1:length(mat_file_dir)
                    if isempty(body_angle_mean{kk,i})
                        continue
                    else
                        plot(dt*segment_idx(onset_segment), body_angle_mean{kk,i}(onset_segment));hold on;
                        body_angle_to_plot(nn,:) = body_angle_mean{kk,i}(onset_segment);
                        nn = nn+1;
                    end
                end
            if(j==1)
                plot(dt*segment_idx(onset_segment), mean(body_angle_to_plot,1,'omitnan'),'r','LineWidth',1.5);
                axis tight;
                ylim =get(gca,'YLim');
                plot([5.5 5.5],[-50 50],'b');
                plot([5 6.5],[0 0],'b');
                to_compare(i,:) = mean(body_angle_to_plot,1,'omitnan');
                title(PAT_NAMES2{PATID(i)},'FontSize',7,'FontWeight','normal','Interpreter','none');
            end
        end
    end
    set(haxes2(1,:),'YLim',[-50 50]);
    set(haxes2(:,2:end),'YTickLabel',[]);
    set(haxes2(1,:),'XTickLabel',[]);
    set(haxes2,'box','off','tickdir','out');
    axes('position',[0 0.95 1 0.001]);set(gca,'XTick',[],'YTick',[],'XColor','w','YColor','w');
    title(['Population (n=' num2str(count_val) ' flies)'],'Interpreter','none');
    %print('-painters','-depsc', 'Magnet_population_v2.eps')

    figure();clf;
    plot(dt*segment_idx(onset_segment), to_compare(1,:),'blue','LineWidth',1.5); hold on
    plot(dt*segment_idx(onset_segment), to_compare(2,:),'red','LineWidth',1.5);
    plot(dt*segment_idx(onset_segment), to_compare(3,:),'green','LineWidth',1.5);
    title('Mean comparison'); xlabel('Time (s)'); ylabel('Body angle (degree)')

% Plotting error bars for amplitude and latency
x = 1:5; %1:11;
y_pre = NaN(numel(peak_ampli_mean), 5); % Preallocate matrix
for jjj = 1:numel(peak_ampli_mean) % Loop over cell elements
    element = peak_ampli_mean{jjj}; % Get current cell element
    element_size = size(element, 2);
    y_pre(jjj, 1:element_size) = element; % Insert cell element into matrix
end
y_pre(y_pre==0) = nan;
y = mean(y_pre,1, 'omitnan');   % Mean of population for each patid, excluding 3 bars patid
yError =  std(y_pre,1, 'omitnan'); % Error values

% Calculate confidence interval
alpha = 0.05; % Significance level (5% confidence)
n = length(x);
t = tinv(1 - alpha/2, n - 1);
ci = t * yError / sqrt(n); % Confidence interval

% Plot the data with error bars
figure();clf;
boxplot(y_pre(:,1:3),'Notch','on','Whisker',1,'Labels',{'Bar','Grating','Spot'}); hold on

jitter = 0.1;

xlabel('Stimulus type');
ylabel('Body angle (degrees)');
title('Peak amplitude (ci95)');
grid on;

% Onset Latency
y_pre2 = NaN(numel(latency_onset_mean), 5); % Preallocate matrix
for jjj = 1:numel(latency_onset_mean) % Loop over cell elements
    element = latency_onset_mean{jjj}; % Get current cell element
    element_size = size(element, 2);
    y_pre2(jjj, 1:element_size) = element; % Insert cell element into matrix
end
y_pre2(y_pre2==0) = nan;
y2 = mean(y_pre2,1, 'omitnan');
yError2 = std(y_pre2,1, 'omitnan'); % Error values
ci2 = t * yError2 / sqrt(n); % Confidence interval

% Plot the data with error bars
figure();clf;
boxplot(y_pre2(:,1:3),'Notch','on','Whisker',1,'Labels',{'Bar','Grating','Spot'}); hold on

xlabel('Stimulus type');
ylabel('Time (s)');
title('Onset 50% latency (ci95)');
%set(barplt2, 'XTickLabel', PAT_NAMES2{3:6})
grid on;