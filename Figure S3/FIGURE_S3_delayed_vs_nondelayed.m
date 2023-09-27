% Figure S3
% Delayed vs Non delayed heading responses
% 2023/09/27 Angel Canelo & Anmo Kim

clear all; clc; close all;

t0=0:0.01:2.5;
hfig=figure();clf;set(gcf,'Color','w'); ii = 1;

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

[tt, yy]=ode45(@(t,y) sim_bar2(t,y,theta_bar,dot_theta_bar,t0),t0, [0 0 0]');
[tt1, yy1]=ode45(@(t,y) sim_spot2(t,y,theta_bar,dot_theta_bar,t0),t0, [-pi 0 0]');
[tt2, yy2]=ode45(@(t,y) sim_grat2(t,y,theta_bar,dot_theta_bar,t0),t0, [0 0 0]');

%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%
    haxes(ii)=subplot(2,1,ii);
    plot(t,y(:,1)/pi*180,'LineWidth',1.2, 'Color', 'blue');
    hold on
    plot(t1,y1(:,1)/pi*180+180,'LineWidth',1.2, 'Color', 'green');
    plot(t2,y2(:,1)/pi*180,'LineWidth',1.2, 'Color', 'red');
    
    plot(tt,yy(:,1)/pi*180,'LineWidth',1.2, 'Color', 'blue', 'LineStyle', '--');
    plot(tt1,yy1(:,1)/pi*180+180,'LineWidth',1.2, 'Color', 'green', 'LineStyle', '--');
    plot(tt2,yy2(:,1)/pi*180,'LineWidth',1.2, 'Color', 'red', 'LineStyle', '--');
    
    plot(t0,theta_bar/pi*180,'black', 'LineWidth',1.2);
    if kk == 1
        legend('fly-bar','fly-spot', 'fly-grating', 'fly-bar (non tau)','fly-spot (non tau)', 'fly-grating (non tau)', 'bar');
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
ii = ii+1;
end



function dydt=sim_spot1(t,y, theta_bar, dot_theta_bar,t0)
theta_bar_now = interp1(t0,theta_bar,t);
dot_theta_bar_now = interp1(t0,dot_theta_bar,t);
I=6*10^-14;%moment of inertia
beta=10^-11;%drag coefficient
tau_bar=0.05;
A = -5e-11; p = 2*pi;
P=@(theta)A*4/p*abs(mod(theta-p/4,p)-p/2)-A;
V=@(theta)5e-12./(1+exp(-2*(3/8*pi-abs(theta))));

dydt=[y(2);
      (-beta*y(2)+y(3))/I;
      (-y(3)+P(theta_bar_now-y(1))+(dot_theta_bar_now-y(2))*V(theta_bar_now-y(1)))/tau_bar];
end

function dydt=sim_bar1(t,y, theta_bar, dot_theta_bar,t0)
theta_bar_now = interp1(t0,theta_bar,t);
dot_theta_bar_now = interp1(t0,dot_theta_bar,t);
I=6*10^-14;%moment of inertia
beta=10^-11;%drag coefficient
tau_bar=0.05;
P=@(theta)+10e-11*sin(theta); % pos to torque e-11
V=@(theta)1e-12*(20-2*theta^2); % vel to torque e-12
% P=@(theta)+8*10^-11*sin(theta);%position function  % pos to torque e-11
% V=@(theta)+16*10^-12./(1+exp(-4*(3/5*pi-abs(theta))));%velocity function % vel to torque e-12

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
%k_opto = 5*10^-11; %velocity function % vel to torque e-10

dydt=[y(2);
      (-beta*y(2)+y(3))/I;
      (-y(3)+(dot_theta_grating_now-y(2))*k_opto)/tau_opto];
end

function dydt=sim_spot2(t,y, theta_bar, dot_theta_bar,t0)
theta_bar_now = interp1(t0,theta_bar,t);
dot_theta_bar_now = interp1(t0,dot_theta_bar,t);
I=6*10^-14;%moment of inertia
beta=10^-11;%drag coefficient
tau_bar=0.05;
A = -5e-11; p = 2*pi;
P=@(theta)A*4/p*abs(mod(theta-p/4,p)-p/2)-A;
V=@(theta)5e-12./(1+exp(-2*(3/8*pi-abs(theta))));

dydt=[y(2);
      (-beta*y(2)+P(theta_bar_now-y(1))+(dot_theta_bar_now-y(2))*V(theta_bar_now-y(1)))/I;
      (-y(3)+P(theta_bar_now-y(1))+(dot_theta_bar_now-y(2))*V(theta_bar_now-y(1)))/tau_bar];
end

function dydt=sim_bar2(t,y, theta_bar, dot_theta_bar,t0)
theta_bar_now = interp1(t0,theta_bar,t);
dot_theta_bar_now = interp1(t0,dot_theta_bar,t);
I=6*10^-14;%moment of inertia
beta=10^-11;%drag coefficient
tau_bar=0.05;
P=@(theta)+10e-11*sin(theta); % pos to torque e-11
V=@(theta)1e-12*(20-2*theta^2); % vel to torque e-12
% P=@(theta)+8*10^-11*sin(theta);%position function  % pos to torque e-11
% V=@(theta)+16*10^-12./(1+exp(-4*(3/5*pi-abs(theta))));%velocity function % vel to torque e-12

dydt=[y(2);
      (-beta*y(2)+P(theta_bar_now-y(1))+(dot_theta_bar_now-y(2))*V(theta_bar_now-y(1)))/I;
      (-y(3)+P(theta_bar_now-y(1))+(dot_theta_bar_now-y(2))*V(theta_bar_now-y(1)))/tau_bar];
end

function dydt=sim_grat2(t,y, theta_bar, dot_theta_bar,t0)
dot_theta_grating_now=interp1(t0,dot_theta_bar,t);
I=6*10^-14;%moment of inertia
beta=10^-11;%drag coefficient
tau_opto=0.04;
k_opto = 25*10^-12; %velocity function % vel to torque e-10
%k_opto = 5*10^-11; %velocity function % vel to torque e-10

dydt=[y(2);
      (-beta*y(2)+(dot_theta_grating_now-y(2))*k_opto)/I;
      (-y(3)+(dot_theta_grating_now-y(2))*k_opto)/tau_opto];
end