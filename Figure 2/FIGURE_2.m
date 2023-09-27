% FIGURE 2

% 2023/07/21 Angel Canelo & Anmo Kim

clear all; clc; close all;

t0=0:0.01:3.5;
hfig=figure();clf;set(gcf,'Color','w'); ii = 1;
fib = figure();clf;set(gcf,'Color','w'); barcolor = [0 0 1; 0 1 0; 1 0 0]; % Blue, Green, Red
for kk = 1:2
    if kk==1        
        theta_bar = 90*pi/180./(1+exp((-t0+0.5)*30));  % for sigmoid-like stimulus
    elseif kk==2
        theta_bar = sin(2*pi*t0);       % for sinusoidal stimulus
    end
dot_theta_bar = diff(theta_bar([1 1:end]))/(t0(2)-t0(1));
    
[t, y]=ode45(@(t,y) sim_bar1(t,y,theta_bar,dot_theta_bar,t0),t0, [0 0 0]');
[t1, y1]=ode45(@(t,y) sim_spot1(t,y,theta_bar,dot_theta_bar,t0),t0, [-pi 0 0]');
[t2, y2]=ode45(@(t,y) sim_grat1(t,y,theta_bar,dot_theta_bar,t0),t0, [0 0 0]');

%%%%%%%%%%%% LATENCIES %%%%%%%%%%%%%%%%%%
%theta_bar_now = interp1(t0,theta_bar,t);
% Calculate the baseline level of the signal
baseline = mean(theta_bar);
% Set the threshold
%threshold = 0 + 0.5*std(theta_bar);
threshold = 0.5*max(theta_bar);
% Detect the onsets
onsets = find(theta_bar >= threshold,1,'first');
t_onset = 0.4;%t0(onsets);
% bar
val50 = 0.5*max(y(:,1));    % 50% of peak value
ind_b = find(y(:,1)>=val50,1,'first');
t50bar = t(ind_b);
latency_bar = abs(t50bar - t_onset);
% spot
val50 = 0.5*max(y1(:,1)+pi);    % 50% of peak value
ind_b = find(y1(:,1)+pi>=val50,1,'first');
t50spot = t1(ind_b);
latency_spot = abs(t50spot - t_onset);% - 0.7; % 0.7 s offset
% grating
val50 = 0.5*max(y(:,1));%0.5*max(y2(:,1));    % 50% of peak value
ind_b = find(y2(:,1)>=val50,1,'first');
t50grat = t2(ind_b);
latency_grat = abs(t50grat - t_onset);
%%%%%%%%%%%% PHASE DIFFERENCE %%%%%%%%%%%
    if kk==2
        Fs = 1/0.01; % sampling frequency (f/time step)
        % bar
        [bar_xcross,lags] = xcorr(theta_bar(1:end), y(1:end,1));
        %stem(lags,phase_bar)
        [~,idx] = max(abs(bar_xcross));
        %latency_bar = 180/pi*(-lags(idx) / (2*pi));
        latency_bar = -lags(idx) / Fs;
        % spot
        [spot_xcross,lags] = xcorr(theta_bar(1:end), y1(1:end,1)+pi);
        [~,idx] = max(abs(spot_xcross));
        %latency_spot = 180/pi*(-lags(idx) / (2*pi));
        latency_spot = -lags(idx) / Fs;% - 0.2;   % 0.2 s offset
        % grating
        [grat_xcross,lags] = xcorr(theta_bar(1:end), y2(1:end,1));
        [~,idx] = max(abs(grat_xcross));
        %latency_grat = 180/pi*(-lags(idx) / (2*pi));
        latency_grat =-lags(idx) / Fs;
    %%%%%%%% Peak to peak amplitude %%%%%%%%%
        peak_bar = 180/pi*(max(y(:,1))-min(y(:,1)));
        peak_spot = 180/pi*(max(y1(:,1))-min(y1(:,1)));
        peak_grat = 180/pi*(max(y2(:,1))-min(y2(:,1)));
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%
    figure(hfig);
    haxes(ii)=subplot(2,1,ii);
    plot(t,y(:,1)/pi*180,'LineWidth',1.2, 'Color', 'blue');
    hold on
    plot(t1,y1(:,1)/pi*180+180,'LineWidth',1.2, 'Color', 'green');
    plot(t2,y2(:,1)/pi*180,'LineWidth',1.2, 'Color', 'red', 'LineStyle', '--');
%     plot(t0,theta_bar/pi*180,'black', 'LineWidth',1.2);
%     plot(t0,yb1/pi*180,'LineWidth',1.2, 'Color', 'blue');
%     hold on
%     plot(t0,ys1/pi*180,'LineWidth',1.2, 'Color', 'green');
    %plot(t0,yg1/pi*180,'LineWidth',1.2, 'Color', 'red', 'LineStyle', '--');
    plot(t0,theta_bar/pi*180,'black', 'LineWidth',1.2);
    if kk == 1
        legend('fly-bar','fly-spot', 'fly-grating', 'bar');
        title('Heading response'); xlabel('Time (s)');
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
    figure(fib);subplot(3,1,kk);
    plt_1 = bar([latency_bar latency_spot latency_grat], 'FaceColor', 'flat'); plt_1.CData = barcolor;
ii = ii+1;
end
subplot(3,1,kk+1);
plt_1 = bar([peak_bar peak_spot peak_grat], 'FaceColor', 'flat'); plt_1.CData = barcolor;
%print(fg3,'-painters','-depsc', 'figure_2_1.eps')

clear all; clc;
global t_gwn; global gwn; global dot_theta_grating;
t_gwn=[0:0.01:10]; %simulation time span, in sec
%gwn=randn(size(t_gwn))*10^-10;%gaussian white noise (arbitrary size)
theta_grating = linspace(-pi,pi/2,size(t_gwn,2));
dot_theta_grating = diff(theta_grating([1 1:end]))/(t_gwn(2)-t_gwn(1));
skip1 = [1 2 3 4]; skip2 = [5 6 7 8];
figure(3);clf;set(gcf,'Color','w'); hold on
color_plt = ['blue green red magenta'];
color_plt = split(color_plt);
for jj=1:4
for i=1:10
    gwn=randn(size(t_gwn))*10^-10;  %gaussian white noise (arbitrary size)
    options = ddeset('MaxStep',0.000125);
    y0=[-pi/1.5 0 0]; %initial condition (\psi=0, \dot\psi=0)
    if jj==1
    [t1, y1]=ode45(@simfly1,[0 t_gwn(end)], y0');
    elseif jj==2
    [t1, y1]=ode45(@simfly2,[0 t_gwn(end)], y0');
    elseif jj==3
    [t1, y1]=ode45(@simfly3,[0 t_gwn(end)], y0');
    elseif jj==4
    [t1, y1]=ode45(@simfly4,[0 t_gwn(end)], y0');
    end
    y{i} = interp1(t1, y1(:,1), t_gwn)';
    y_mod{i} = mod(y{i}+pi,2*pi)-pi;
    bCrossings=abs(diff(y_mod{i}))>1.9*pi;
    y_mod{i}(bCrossings)=nan;
end
y2 = cell2mat(y);
for j=1:size(t_gwn,2)
    hd_mean(j) = mean(y2(j,:), 2);
end

%plot the simulation result: psi(t) and its histogram
psi=mod(hd_mean+pi,2*pi)-pi;
bCrossings=abs(diff(psi))>1.9*pi;
psi(bCrossings)=nan;
subplot(2,4,skip1(jj));
hold on; for kk=1:10 p(kk) = plot(y_mod{kk},t_gwn, 'Color', uint8([17 17 17]), 'LineWidth',0.2); p(kk).Color(4) = 0.25; end
p(kk+1) = plot(psi,t_gwn, color_plt{jj},'LineWidth',3);
%p(kk+2) = plot(dot_theta_grating, t_gwn, 'r','LineWidth',2);
ylim=get(gca,'Ylim');
plot([0 0],ylim,'k--');
set(gca,'Box','off','Tickdir','out','FontSizE',12,'XLim',[-pi-.02 pi],'XTick',[-pi:pi/2:pi],...
    'XTickLabel',{'-\pi','-\pi/2','0','\pi/2','\pi'},'YDir','reverse','YTick',[0:5:10], 'YTickLabel',{0:5:10});
ylabel('time (s)');xlabel('heading (rad)');
title(['Mean response for 10 simulations']); legend([p(kk+1)],'mean response');
subplot(2,4,skip2(jj));
h = histfit(cell2mat(cellfun(@(x)x(:),y_mod(:),'un',0)),20,'kernel');
delete(h(1)); set(h(2),'color',color_plt{jj})
ylim=get(gca,'YLim');
hline = line([0 0],ylim); hline.Color = 'black'; hline.LineStyle = '--';
yt = get(gca, 'YTick');
set(gca,'FontSize',12,'Box','off','Tickdir','out','XLim',[-pi pi],'XTick',[-pi:pi/2:pi],...
    'XTickLabel',{'-\pi','-\pi/2','0','\pi/2','\pi'}, 'YTick', yt, 'YTickLabel', yt/(numel(y)*100));
xlabel('heading (rad)');
title('Probability density');
end

%Plotting position and velocity components
%%%%%
psi=linspace(-pi,pi,100);
%%%%%%%%%%%%%%% Bar %%%%%%%%%%%%%%%%%%%%%%%
pos_fcn_bar = 10*sin(psi);%position function
vel_fcn_bar = 20-2*psi.^2;
%%%%%%%%%%%%%%% Spot %%%%%%%%%%%%%%%%%%%%%%%
A = -5; x = psi; p = 2*pi;
pos_fcn_spot = A*4/p*abs(mod(x-p/4,p)-p/2)-A;
vel_fcn_spot = 5./(1+exp(-2*(3/8*pi-abs(psi))));%velocity function
%%%%%%%%%%%%%%% Grating %%%%%%%%%%%%%%%%%%%%
pos_fcn_grating = zeros(size(psi,2));%position function
vel_fcn_grating = 25*ones(size(psi,2));%position function
%%%%%
psi = 180/pi*psi;
figure(); clf;
subplot(2,1,1);
plot(psi,pos_fcn_spot,'green','LineWidth',2); hold on;
plot(psi,pos_fcn_bar,'blue','LineWidth',2);
plot(psi,pos_fcn_grating,'red','LineWidth',2);
axis tight;
%grid on
set(gca,'FontSize',12,'TickDir','out','Box','off');
set(gca,'YLim',[-10 10]);
legend('spot','bar','grating');
title(gca,'Position function');
%xlabel(gca,'Position of the object (rad)');
ylabel(gca,'Position func(deg)');
subplot(2,1,2);
plot(psi,vel_fcn_spot,'green','LineWidth',2); hold on;
plot(psi,vel_fcn_bar,'blue','LineWidth',2);
plot(psi,vel_fcn_grating,'red','LineWidth',2);
axis tight;
%grid on
set(gca,'FontSize',12,'TickDir','out','Box','off');
set(gca,'YLim',[0 30]);
legend('spot','bar','grating');
title(gca,'Velocity function');
xlabel(gca,'Pattern position(degree)');
ylabel(gca,'Velocity func (deg)');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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

function dydt=simfly1(t,y)
global t_gwn; global gwn;
Theta=0.6*10^-13;%kg m^2 -- moment of inertia (following Ristroph et al 2010)
k=1.0*10^-11;%kg m^2 s^-1 -- friction constant (following Ristroph et al 2010)
N_now=interp1(t_gwn,gwn,t);
P=@(theta)+10e-11*sin(theta); % pos to torque e-11
V=@(theta)1e-12*(20-2*theta^2); % vel to torque e-12
tau_bar=0.05;

dydt=[y(2);
      (-k*y(2)-y(3)+N_now)/Theta;
      (-y(3)+P(y(1))+(y(2))*V(y(1)))/tau_bar];
end

function dydt=simfly2(t,y)
global t_gwn; global gwn;
Theta=0.6*10^-13;%kg m^2 -- moment of inertia (following Ristriph et al 2010)
k=1.0*10^-11;%kg m^2 s^-1 -- friction constant (following Ristriph et al 2010)
N_now=interp1(t_gwn,gwn,t);
A = -5e-11; p = 2*pi;
P=@(theta)A*4/p*abs(mod(theta-p/4,p)-p/2)-A;
V=@(theta)5e-12./(1+exp(-2*(3/8*pi-abs(theta))));
tau_bar=0.05;

dydt=[y(2);
      (-k*y(2)-y(3)+N_now)/Theta;
      (-y(3)+P(y(1))+(y(2))*V(y(1)))/tau_bar];
end

function dydt=simfly3(t,y)
global t_gwn; global gwn;
Theta=0.6*10^-13;%kg m^2 -- moment of inertia (following Ristroph et al 2010)
k=1.0*10^-11;%kg m^2 s^-1 -- friction constant (following Ristroph et al 2010)
N_now=interp1(t_gwn,gwn,t);
k_opto = 25*10^-12; %velocity function % vel to torque e-10
tau_opto=0.04;

dydt=[y(2);
      (-k*y(2)-y(3)+N_now)/Theta;
      (-y(3)+y(2)*k_opto)/tau_opto];
end

function dydt=simfly4(t,y)
global t_gwn; global gwn;
Theta=0.6*10^-13;%kg m^2 -- moment of inertia (following Ristroph et al 2010)
k=1.0*10^-11;%kg m^2 s^-1 -- friction constant (following Ristroph et al 2010)
N_now=interp1(t_gwn,gwn,t);
k_opto = 25*10^-12; %velocity function % vel to torque e-10
tau_opto=0.04;

dydt=[y(2);
      (-k*y(2)-0*y(3)+N_now)/Theta;
      (-y(3)+y(2)*k_opto)/tau_opto];
end