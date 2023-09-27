% FIGURE 1
% Deriving POS VEL
% 2023/07/21    Angel Canelo & Anmo Kim

clear all; clc; close all;

save_pos_vel = ['POSVEL_data\position_velocity_components_spot.mat POSVEL_data\position_velocity_components_bar.mat POSVEL_data\position_velocity_components_grating.mat'];
save_pos_vel = split(save_pos_vel);

V2deg = 135/5;

lbound = [14 20 6]; 
hbound = [15 21 8]; % Bounds to pick between bar, spot, grating in pattern id
bbound = [13 19 5]; % Base bound
lboundcc = [-14 -20 -6]; 
hboundcc = [-15 -21 -8]; % Bounds to pick between bar, spot, grating in pattern id
bboundcc = [-13 -19 -5];
hz = [1 2 3];
vt = [4 5 6];
plotplot = [5 6 1];     % choose fly with similar response than population
plotplot2 = [6 6 1];
fg = figure; clf; fg3 = figure; clf;

folder = '..\data\exp\**\';
filen = dir([folder, '\*.EDR']);
downs = 500;% 500 or 1000 to save as eps;    % downsampling scale for data storage
kk=0;
for jj=1:length(filen)
    if ~contains(filen(jj).folder, 'discards')
        kk = kk+1;
        myfile{kk} = filen(jj).name;
        myfolder{kk} = filen(jj).folder;
    end
end

flyn_plt = [6 6 5];
tit = ["spot" "bar" "grating"];
for k=1:3   % (k=pattern, j=fly number)
    clear wba_tr; clear wba_m; clear lwba_tr; clear rwba_tr; clear patid_tr;
    clear wbacc_tr; clear wbacc_m; clear pos; clear vel;
wba_cell = []; wbacc_cell = []; time_tr = [];
jjj = 1; m = 1; jjjcc = 1;
for j=1:kk
    [data,h] = import_edr([myfolder{j} '\' myfile{j}]);
    dt=h.DT;
    t=(0:(size(data,1)-1))*dt;
    xpos{j}=data(:,2);
    lwba{j}= data(:,4);
    rwba{j}= data(:,5);
    wba{j} = data(:,4) - data(:,5);     % Clockwise wba
    patid{j}=round((data(:,6)+0.2)*3);    % normalizing pattern id
    patidcc{j}=round((data(:,6)-0.2)*3);

    onset_idx=find(patid{j}([1 1:(end-1)])<=bbound(k) & (patid{j}>=lbound(k) & patid{j}<=hbound(k)));
    onsetcc_idx=find(patidcc{j}([1 1:(end-1)])>=bboundcc(k) & (patidcc{j}<=lboundcc(k) & patidcc{j}>=hboundcc(k)));

    dur=10.1;   % from posvel file
    trial_idx=round(-1/dt):(dur+1)/dt;   %round(-1/dt):(dur+1)/dt;   % block duration (open loop)
    %%%%% Clockwise %%%%%
    jj=1;valid=0;
    for i=1:length(onset_idx)
        if lwba{j}(trial_idx+onset_idx(i))~=0 & rwba{j}(trial_idx+onset_idx(i))~=0 % Checking blocks with non zero elements
            if patid{j}(trial_idx(end)+2+onset_idx(i))<lbound(k) % Checking open loop block
                jj = jj+1;
                if jj == 4      % Check whether there is more than 3 valid blocks
                    valid = 1;
                    break
                end
            end
        end
    end
    jj=1;
    for i=1:length(onset_idx)
        if lwba{j}(trial_idx+onset_idx(i))~=0 & rwba{j}(trial_idx+onset_idx(i))~=0 % Checking blocks with non zero elements
            if patid{j}(trial_idx(end)+2+onset_idx(i))<lbound(k) % Checking open loop block
                if valid==1
                    lwba_tr(jjj,jj,:)=lwba{j}(trial_idx+onset_idx(i))';    % Saving wingbeat for each file (j) and each onset interval (i)
                    rwba_tr(jjj,jj,:)=rwba{j}(trial_idx+onset_idx(i))';
                    wba_tr(jjj,jj,:)=bpf_ak(wba{j}(trial_idx+onset_idx(i))',dt, 0, 20);
                    patid_tr(jjj,jj,:)=patid{j}(trial_idx+onset_idx(i));
                    xpos_tr(jjj,jj,:)=xpos{j}(trial_idx+onset_idx(i))';
                    jj = jj+1;
                    if jj==7
                        plt(k,m) = jjj;     % Choose a fly to plot with at least 6 blocks
                        m = m+1;
                    end
                end
            end
        end
    end
    if valid == 1 % To not fill with empty elements in the middle of the array (For some files it might not to enter here)
        jjj = jjj+1;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Counter-Clockwise %%%%%
        jjcc=1;validcc=0;
    for i=1:length(onsetcc_idx)
        if lwba{j}(trial_idx+onsetcc_idx(i))~=0 & rwba{j}(trial_idx+onsetcc_idx(i))~=0 % Checking blocks with non zero elements
            if patid{j}(trial_idx(end)+2+onsetcc_idx(i))>lboundcc(k) % Checking open loop block
                jjcc = jjcc+1;
                if jjcc == 4      % Check whether there is more than 3 valid blocks
                    validcc = 1;
                    break
                end
            end
        end
    end
    jjcc=1;
    for i=1:length(onsetcc_idx)
        if lwba{j}(trial_idx+onsetcc_idx(i))~=0 & rwba{j}(trial_idx+onsetcc_idx(i))~=0 % Checking blocks with non zero elements
            if patid{j}(trial_idx(end)+2+onsetcc_idx(i))>lboundcc(k) % Checking open loop block
                if validcc==1
                    wbaccc_tr(jjjcc,jjcc,:)= bpf_ak(wba{j}(trial_idx+onsetcc_idx(i))',dt, 0, 20);
                    jjcc = jjcc+1;
                end
            end
        end
    end
    if validcc == 1 % To not fill with empty elements in the middle of the array (For some files it might not to enter here)
        jjjcc = jjjcc+1;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
end
until = min(size(wba_tr,1),size(wbaccc_tr,1));
for j=1:until
    % Computing mean for each fly
    if length(find(lwba_tr(j,:,1)~=0))==1   % Avoiding to compute mean in the case of single element
        wba_m(j,:) = wba_tr(j,find(wba_tr(j,:,1)~=0),:);    % We use "find" to pick non zero blocks only
    else
        wba_m(j,:) = mean(squeeze(wba_tr(j,find(wba_tr(j,:,1)~=0),:)),1);
    end
end
for j=1:until  % (k=pattern, j=fly number)
    % Computing mean for each fly
    if length(find(wbaccc_tr(j,:,1)~=0))==1   % Avoiding to compute mean in the case of single element
        wbaccc_m(j,:) = wbaccc_tr(j,find(wbaccc_tr(j,:,1)~=0),:);
    else
        wbaccc_m(j,:) = mean(squeeze(wbaccc_tr(j,find(wbaccc_tr(j,:,1)~=0),:)),1);
    end
%     if k==1
        for_spot = [wba_m(j,1:downs:end); -wbaccc_m(j,1:downs:end)];
        for_spot = mean(for_spot,1);
        pos(k,j,:) = (for_spot+flip(-for_spot))/2;
        vel(k,j,:) = (for_spot-flip(-for_spot))/2;
%     else
%         pos(k,j,:) = (wba_m(j,1:downs:end)+flip(wbaccc_m(j,1:downs:end)))/2;
%         vel(k,j,:) = (wba_m(j,1:downs:end)-flip(wbaccc_m(j,1:downs:end)))/2;
%     end
end
wba_tm = mean(wba_m,1);     % Total mean out of 20 flies
wbaccc_tm = mean(wbaccc_m,1);

pos_tm = mean(squeeze(pos(k,:,:)),1);
vel_tm = mean(squeeze(vel(k,:,:)),1);

% Plotting
flyn = plt(k,plotplot(k)); flyn2 = plt(k,plotplot2(k));

%%%%% CONFIDENCE INTERVAL %%%%%
wba_std = std(V2deg*squeeze(wba_tr(flyn,find(wba_tr(flyn,:,1)~=0),1:downs:end)),0,1);
nel = size(squeeze(wba_tr(flyn,find(wba_tr(flyn,:,1)~=0),1:downs:end)),1);
x = trial_idx(1:downs:end)*dt;
conf = V2deg*wba_m(flyn,1:downs:end) - 1.96*(wba_std/sqrt(nel));
conf2 = V2deg*wba_m(flyn,1:downs:end) + 1.96*(wba_std/sqrt(nel));

wba_std_t = std(V2deg*wba_m(:,1:downs:end),0,1);
nel2 = size(wba_m(:,1:downs:end),1);
conft = V2deg*wba_tm(1:downs:end) - 1.96*(wba_std_t/sqrt(nel2));
conft2 = V2deg*wba_tm(1:downs:end) + 1.96*(wba_std_t/sqrt(nel2));

wbaccc_std = std(V2deg*squeeze(wbaccc_tr(flyn2,find(wbaccc_tr(flyn2,:,1)~=0),1:downs:end)),0,1);
nelcc = size(squeeze(wbaccc_tr(flyn2,find(wbaccc_tr(flyn2,:,1)~=0),1:downs:end)),1);
confcc = V2deg*wbaccc_m(flyn2,1:downs:end) - 1.96*(wbaccc_std/sqrt(nelcc));
confcc2 = V2deg*wbaccc_m(flyn2,1:downs:end) + 1.96*(wbaccc_std/sqrt(nelcc));

wbaccc_std_t = std(V2deg*wbaccc_m(:,1:downs:end),0,1);
nelcc2 = size(wbaccc_m(:,1:downs:end),1);
conftcc = V2deg*wbaccc_tm(1:downs:end) - 1.96*(wbaccc_std_t/sqrt(nelcc2));
conftcc2 = V2deg*wbaccc_tm(1:downs:end) + 1.96*(wbaccc_std_t/sqrt(nelcc2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clockwise wba
figure(fg);
haxes3=subplot(2,3,hz(k)); hold on; ylim([-50 50])
plot(trial_idx(1:downs:end)*dt,V2deg*wba_m(flyn,1:downs:end),'r','Linewidth',2)
patch([x fliplr(x)], [conf  fliplr(conf2)], 'red','FaceColor','red','FaceAlpha',.3, 'LineStyle', 'none')
plot(trial_idx(1:downs:end)*dt,V2deg*wbaccc_m(flyn2,1:downs:end),'blue','Linewidth',2)
patch([x fliplr(x)], [confcc  fliplr(confcc2)], 'blue','FaceColor','blue','FaceAlpha',.3, 'LineStyle', 'none')
title(sprintf('Single fly >=3 trials (%s)',tit(k)));
ylabel('L-R wba (degree)');
set(gca,'TickDir','out','Box','off');
haxes4=subplot(2,3,vt(k)); hold on; ylim([-50 50])
plot(trial_idx(1:downs:end)*dt,V2deg*wba_tm(1:downs:end),'r','Linewidth',2)
patch([x fliplr(x)], [conft  fliplr(conft2)], 'red','FaceColor','red','FaceAlpha',.3, 'LineStyle', 'none')
plot(trial_idx(1:downs:end)*dt,V2deg*wbaccc_tm(1:downs:end),'blue','Linewidth',2)
patch([x fliplr(x)], [conftcc  fliplr(conftcc2)], 'blue','FaceColor','blue','FaceAlpha',.3, 'LineStyle', 'none')
title(sprintf('Population of flies (%s)',tit(k))); xlabel('Time (s)')
ylabel('L-R wba (degree)'); sgtitle('Wing response (L-R WBA)')
set(gca,'TickDir','out','Box','off');

%%%%% CONFIDENCE INTERVAL %%%%%
psi = linspace(-180,180,length(trial_idx(1:downs:end)));
pos_std = std(V2deg*squeeze(pos(k,:,:)),0,1);
x = psi;
nell = size(squeeze(pos(k,:,:)),1);
confp = V2deg*pos_tm - 1.96*(pos_std/sqrt(nell));
confp2 = V2deg*pos_tm + 1.96*(pos_std/sqrt(nell));

vel_std = std(V2deg*squeeze(vel(k,:,:)),0,1);
nell2 = size(squeeze(vel(k,:,:)),1);
confv = V2deg*vel_tm - 1.96*(vel_std/sqrt(nell2));
confv2 = V2deg*vel_tm + 1.96*(vel_std/sqrt(nell2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Position and Velocity functions
figure(fg3);
haxesp(1)=subplot(2,3,hz(k)); hold on; ylim([-15 15])
plot(psi,V2deg*pos_tm,'magenta','Linewidth',2)
patch([x fliplr(x)], [confp  fliplr(confp2)], 'red','FaceColor','magenta','FaceAlpha',.3, 'LineStyle', 'none')
title(sprintf('Position component (%s)',tit(k)));
ylabel('P (degree)');
set(gca,'TickDir','out','Box','off','XTick',-180:90:180,...
    'XTickLabel',{'-\pi','-\pi/2','0','\pi/2','\pi'});
haxesp(2)=subplot(2,3,vt(k)); hold on; ylim([-10 30])
plot(psi,V2deg*vel_tm,'magenta','Linewidth',2)
patch([x fliplr(x)], [confv  fliplr(confv2)], 'red','FaceColor','magenta','FaceAlpha',.3, 'LineStyle', 'none')
title(sprintf('Velocity component (%s)',tit(k))); xlabel('Pattern position (degree)')
ylabel('V (degree)'); sgtitle('Position and Velocity components')
set(gca, 'TickDir','out','Box','off','XTick',-180:90:180,...
    'XTickLabel',{'-\pi','-\pi/2','0','\pi/2','\pi'});
end