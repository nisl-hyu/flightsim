% FIGURE 1
% One trial of a single fly
% 2023/07/21    Angel Canelo & Anmo Kim

clear all; clc; close all;

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
fg = figure; clf;

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
                    lwba_tr(jjj,jj,:)=bpf_ak(lwba{j}(trial_idx+onset_idx(i))',dt, 0, 20);
                    rwba_tr(jjj,jj,:)=bpf_ak(rwba{j}(trial_idx+onset_idx(i))',dt, 0, 20);
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
    pos(k,j,:) = (wba_m(j,1:downs:end)+flip(wbaccc_m(j,1:downs:end)))/2;
    vel(k,j,:) = (wba_m(j,1:downs:end)-flip(wbaccc_m(j,1:downs:end)))/2;
end
wba_tm = mean(wba_m,1);     % Total mean out of 20 flies
wbaccc_tm = mean(wbaccc_m,1);

pos_tm = mean(squeeze(pos(k,:,:)),1);
vel_tm = mean(squeeze(vel(k,:,:)),1);

% Plotting
flyn = plt(k,plotplot(k)); flyn2 = plt(k,plotplot2(k));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clockwise wba
figure(fg);
haxes3=subplot(1,3,hz(k)); hold on; %ylim([-150 100])
plot(trial_idx(1:downs:end)*dt,V2deg*squeeze(lwba_tr(flyn,1,1:downs:end)),'black','Linewidth',1)
plot(trial_idx(1:downs:end)*dt,V2deg*squeeze(rwba_tr(flyn,1,1:downs:end)),'black','Linewidth',1)
plot(trial_idx(1:downs:end)*dt,V2deg*squeeze(wba_tr(flyn,1,1:downs:end)),'red','Linewidth',1)
title(sprintf('Single fly >=3 trials (%s)',tit(k)));
ylabel('L-R wba (degree)');
set(gca,'TickDir','out','Box','off');

end