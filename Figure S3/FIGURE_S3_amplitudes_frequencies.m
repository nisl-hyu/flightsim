% Figure S3
% Heading responses for different stimulus amplitudes and frequencies
% 2023/09/27 Angel Canelo & Anmo Kim

clear all; clc; close all;

t0=0:0.01:3.5;
hfig=figure();clf;set(gcf,'Color','w'); ii = 1;
fib = figure();clf;set(gcf,'Color','w'); barcolor = [0 0 1; 0 1 0; 1 0 0]; % Blue, Green, Red
ampl = [60 90 120 150];
freq = [0.5 1 1.5 2];
for jj = 1:4
for kk = 1:2
    if kk==1        
        theta_bar = ampl(jj)*pi/180./(1+exp((-t0+0.5)*30));  % for sigmoid-like stimulus
    elseif kk==2
        theta_bar = sin(2*pi*freq(jj)*t0);       % for sinusoidal stimulus
    end
dot_theta_bar = diff(theta_bar([1 1:end]))/(t0(2)-t0(1));
    
[t, y]=ode45(@(t,y) sim_bar1(t,y,theta_bar,dot_theta_bar,t0),t0, [0 0 0]');
[t1, y1]=ode45(@(t,y) sim_spot1(t,y,theta_bar,dot_theta_bar,t0),t0, [-pi 0 0]');
[t2, y2]=ode45(@(t,y) sim_grat1(t,y,theta_bar,dot_theta_bar,t0),t0, [0 0 0]');

%%%%%%%%%%%% LATENCIES %%%%%%%%%%%%%%%%%%
    if kk ==1
        %theta_bar_now = interp1(t0,theta_bar,t);
        % Calculate the baseline level of the signal
        baseline = mean(theta_bar);
        % Set the threshold
        %threshold = 0 + 0.5*std(theta_bar);
        threshold = 0.5*max(theta_bar);
        % Detect the onsets
        onsets = find(theta_bar >= threshold,1,'first');
        t_onset = t0(onsets);
        % bar
        val50 = 0.5*max(y(:,1));    % 50% of peak value
        ind_b = find(y(:,1)>=val50,1,'first');
        t50bar = t(ind_b);
        latency_bar(jj) = abs(t50bar - t_onset);
        % spot
        val50 = 0.5*max(y1(:,1)+pi);    % 50% of peak value
        ind_b = find(y1(:,1)+pi>=val50,1,'first');
        t50spot = t1(ind_b);
        latency_spot(jj) = abs(t50spot - t_onset);% - 0.6;
        % grating
        val50 = 0.5*max(y2(:,1));    % 50% of peak value
        ind_b = find(y2(:,1)>=val50,1,'first');
        t50grat = t2(ind_b);
        latency_grat(jj) = abs(t50grat - t_onset);
    end
%%%%%%%%%%%% PHASE DIFFERENCE %%%%%%%%%%%
    if kk==2
        Fs = 1/0.01; % sampling frequency (f/time step)
        % bar
        [bar_xcross,lags] = xcorr(theta_bar(1:end), y(1:end,1));
        %stem(lags,phase_bar)
        [~,idx] = max(abs(bar_xcross));
        %latency_bar = 180/pi*(-lags(idx) / (2*pi));
        phase_bar(jj) = -lags(idx) / Fs;
        % spot
        [spot_xcross,lags] = xcorr(theta_bar(1:end), y1(1:end,1)+pi);
        [~,idx] = max(abs(spot_xcross));
        %latency_spot = 180/pi*(-lags(idx) / (2*pi));
        phase_spot(jj) = -lags(idx) / Fs;
        % grating
        [grat_xcross,lags] = xcorr(theta_bar(1:end), y2(1:end,1));
        [~,idx] = max(abs(grat_xcross));
        %latency_grat = 180/pi*(-lags(idx) / (2*pi));
        phase_grat(jj) =-lags(idx) / Fs;
    %%%%%%%% Peak to peak amplitude %%%%%%%%%
        peak_bar(jj) = 180/pi*(max(y(:,1))-min(y(:,1)));
        peak_spot(jj) = 180/pi*(max(y1(:,1))-min(y1(:,1)));
        peak_grat(jj) = 180/pi*(max(y2(:,1))-min(y2(:,1)));
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%
    figure(hfig);
    haxes(ii)=subplot(4,2,ii);
    plot(t,y(:,1)/pi*180,'LineWidth',1.2, 'Color', 'blue');
    hold on
    plot(t1,y1(:,1)/pi*180+180,'LineWidth',1.2, 'Color', 'green');
    plot(t2,y2(:,1)/pi*180,'LineWidth',1.2, 'Color', 'red', 'LineStyle', '--');
    ylim([-60 150]);
%     plot(t0,theta_bar/pi*180,'black', 'LineWidth',1.2);
%     plot(t0,yb1/pi*180,'LineWidth',1.2, 'Color', 'blue');
%     hold on
%     plot(t0,ys1/pi*180,'LineWidth',1.2, 'Color', 'green');
    %plot(t0,yg1/pi*180,'LineWidth',1.2, 'Color', 'red', 'LineStyle', '--');
    plot(t0,theta_bar/pi*180,'black', 'LineWidth',1.2);
    if kk == 1
        legend('fly-bar','fly-spot', 'fly-grating', 'bar');
        title('Body angle'); xlabel('Time (s)');
    end
    set(0,'DefaultLegendAutoUpdate','off')
    ylabel('(deg)'); %set(gca,'TickDir','out','Box','off','YTick',-70:35:70);
    %set(haxes,'YLim',[-70 70]);
    set(haxes,'Box','off','TickDir','out');
    %set(haxes(1:2),'XTickLabel',[]);
    for i=1:length(haxes)
        y_lim=get(haxes(i),'YLim');
        hpatch=patch(haxes(i),[0.4 0.6 0.6 0.4],y_lim([1 1 2 2]),'r','FaceColor',ones(1,3)*.9,'EdgeColor','none');
        uistack(hpatch,'bottom');
    end
ii = ii+1;
end
end
figure(fib);subplot(3,1,1);
plot(ampl, latency_bar, 'blue'); hold on; plot(ampl, latency_spot, 'green');plot(ampl, latency_grat, 'red');
subplot(3,1,2);
plot(freq, peak_bar, 'blue'); hold on; plot(freq, peak_spot, 'green');plot(freq, peak_grat, 'red');
subplot(3,1,3);
plot(freq, phase_bar, 'blue'); hold on; plot(freq, phase_spot, 'green');plot(freq, phase_grat, 'red');

function dydt=sim_spot1(t,y, theta_bar, dot_theta_bar,t0)
theta_bar_now = interp1(t0,theta_bar,t);
dot_theta_bar_now = interp1(t0,dot_theta_bar,t);
I=6*10^-14;%moment of inertia
beta=10^-11;%drag coefficient
tau_spot=0.05;
A = -5e-11; p = 2*pi;
P=@(theta)A*4/p*abs(mod(theta-p/4,p)-p/2)-A;
V=@(theta)5e-12./(1+exp(-2*(3/8*pi-abs(theta))));

dydt=[y(2);
      (-beta*y(2)+y(3))/I;
      (-y(3)+P(theta_bar_now-y(1))+(dot_theta_bar_now-y(2))*V(theta_bar_now-y(1)))/tau_spot];
end

function dydt=sim_bar1(t,y, theta_bar, dot_theta_bar,t0)
theta_bar_now = interp1(t0,theta_bar,t);
dot_theta_bar_now = interp1(t0,dot_theta_bar,t);
I=6*10^-14;%moment of inertia
beta=10^-11;%drag coefficient
tau_bar=0.05;
P=@(theta)+10e-11*sin(theta); % pos to torque e-11
V=@(theta)1e-12*(20-2*theta^2); % vel to torque e-12

dydt=[y(2);
      (-beta*y(2)+y(3))/I;
      (-y(3)+P(theta_bar_now-y(1))+(dot_theta_bar_now-y(2))*V(theta_bar_now-y(1)))/tau_bar];
end

function dydt=sim_grat1(t,y, theta_bar, dot_theta_bar,t0)
dot_theta_grating_now=interp1(t0,dot_theta_bar,t);
I=6*10^-14;%moment of inertia
beta=10^-11;%drag coefficient
tau_opto=0.04;
k_opto = 25*10^-12; %velocity function % vel to torque e-10

dydt=[y(2);
      (-beta*y(2)+y(3))/I;
      (-y(3)+(dot_theta_grating_now-y(2))*k_opto)/tau_opto];
end