% FIGURE 4

% 2023/07/25 Angel Canelo & Anmo Kim

clear all; close all; clc;
global t0 theta_bar dot_theta_bar dot_theta_grating supp_global
t0=0:0.001:1.6; % 0:0.001:3.5;

for kk = 1:1
    if kk==1        
        theta_bar = 60*pi/180./(1+exp((-t0+0.5)*30));  % for sigmoid-like stimulus
        theta_grating = 0*pi/3./(1+exp((-t0+0.5)*30));  %30deg righwards
    elseif kk==2
        %theta_bar = 0*pi/3./(1+exp((-t0+0.5)*30));       % for sinusoidal stimulus
        theta_bar = 60*pi/180./(1+exp((-t0+2)*30));%-60*pi/180./(1+exp((-t0+2.25)*30));  % for sigmoid-like stimulus
        theta_grating = pi/3./(1+exp((-t0+0.5)*30))+pi/3./(1+exp((-t0+2)*30));  %30deg righwards
    end
dot_theta_bar = diff(theta_bar([1 1:end]))/(t0(2)-t0(1));
dot_theta_grating = diff(theta_grating([1 1:end]))/(t0(2)-t0(1));

[t, y]=ode45(@sim_additive1,[0 t0(end)], [0 0 0 0]');
[t2, y2]=ode45(@sim_ek1,[0 t0(end)], [0 0 0 0 0]');
[t1, y1]=ode45(@sim_bar1,[0 t0(end)], [0 0 0]');
[t3, y3]=ode45(@sim_grat1,[0 t0(end)], [0 0 0]');
[t4, y4]=ode45(@sim_all_or_non,[0 t0(end)], [0 0 0 0]');

t_supp = linspace(t4(1),t4(end),length(supp_global));
supp_global_now=interp1(t_supp,supp_global,t4);

%plotting
hfig=figure(kk);clf;set(gcf,'Color','w');
haxes(1)=subplot(231);
plot(t,y(:,1)/pi*180,'magenta','LineWidth',1.2); hold on
plot(t1,y1(:,1)/pi*180,'--b', 'LineWidth',1.2);
plot(t0,theta_bar/pi*180,'black','LineWidth',1.2)
plot(t0,theta_grating/pi*180,'--k','LineWidth',1.2);
title('fly angular position');ylabel('angle(deg)'); xlabel('time (s)'); 
haxes(2)=subplot(232);
plot(t2,y2(:,1)/pi*180, 'Color',[0.9290 0.6940 0.1250],'LineWidth',1.2); hold on
plot(t0,theta_bar/pi*180,'black','LineWidth',1.2)
plot(t0,theta_grating/pi*180,'--k','LineWidth',1.2);
haxes(3)=subplot(233);
%plot(t3,y3(:,1)/pi*180,'red','LineWidth',1.2);
plot(t4,y4(:,1)/pi*180,'Color', [200 100 16]/256,'LineWidth',1.2); hold on
plot(t0,theta_bar/pi*180,'black','LineWidth',1.2)
plot(t0,theta_grating/pi*180,'--k','LineWidth',1.2);
%legend('efference copy','additive','to bar','to grating','suppression','bar','grating'); set(0,'DefaultLegendAutoUpdate','off')
set(haxes(1:3),'YLim',[-10 90]);

haxes(4)=subplot(234); hold on;
plot(t,y(:,3),'magenta','LineWidth',1.2);
plot(t,y(:,4),'--m','LineWidth',1.2);
xlabel('time (s)'); ylabel('Torque (Nm)'); hold on
haxes(5)=subplot(235); hold on;
plot(t2,y2(:,3),'Color',[0.9290 0.6940 0.1250],'LineWidth',1.2);
plot(t2,y2(:,4),'Color',[0.9290 0.6940 0.1250],'LineStyle','--','LineWidth',1.2);
plot(t2,3e-13*(y2(:,5)/pi*180),'Color','black','LineWidth',1.2);
plot(t2,3e-13*(dot_theta_grating(1:1:end)/pi*180-y2(1:1:end,2)/pi*180),'Color','black','LineStyle','--','LineWidth',1.2);
haxes(6)=subplot(236); hold on;
plot(t4,y4(:,3),'Color', [200 100 16]/256,'LineWidth',1.2);
plot(t4,y4(:,4),'Color', [200 100 16]/256,'LineStyle','--','LineWidth',1.2);
plot(t4, supp_global_now(1:1:end),'black','LineWidth',1.2);
%legend('torque to bar (additive)','torque to grating (additive)','torque to bar (efference copy)','torque to grating (efference copy)', 'torque to bar (suppression)','torque to grating (suppression)');
set(haxes,'Box','off','TickDir','out');
set(haxes(4:6),'YLim',[-1.0e-10 1.5e-10]);
%set(haxes(4:6),'XTickLabel',[]);
for i=1:length(haxes)
    ylim=get(haxes(i),'YLim');
    hpatch=patch(haxes(i),[0.4 0.6 0.6 0.4],ylim([1 1 2 2]),'r','FaceColor',ones(1,3)*.9,'EdgeColor','none');
    %hpatch2=patch(haxes(i),[1.9 2.1 2.1 1.9],ylim([1 1 2 2]),'r','FaceColor',ones(1,3)*.9,'EdgeColor','none');
    uistack(hpatch,'bottom');
    %uistack(hpatch2,'bottom');
end
% print(hfig,'-painters','-depsc', 'fig42.eps')
supp_global = [];
end

function dydt=sim_bar1(t,y)
global t0 theta_bar dot_theta_bar;
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

function dydt=sim_grat1(t,y)
global t0 dot_theta_grating;
dot_theta_grating_now=interp1(t0,dot_theta_grating,t);
I=6*10^-14;%moment of inertia
beta=10^-11;%drag coefficient
tau_opto=0.04;
k_opto = 25*10^-12;

dydt=[y(2);
      (-beta*y(2)+y(3))/I;
      (-y(3)+(dot_theta_grating_now-y(2))*k_opto)/tau_opto];
end

function dydt=sim_additive1(t,y)
global t0 theta_bar dot_theta_bar dot_theta_grating;
theta_bar_now=interp1(t0,theta_bar,t);
dot_theta_bar_now=interp1(t0,dot_theta_bar,t);
dot_theta_grating_now=interp1(t0,dot_theta_grating,t);
I=6*10^-14;%moment of inertia
beta=10^-11;%drag coefficient
tau_bar=0.05;
tau_opto=0.04;
P=@(theta)+10e-11*sin(theta); % pos to torque e-11
V=@(theta)1e-12*(20-2*theta^2); % vel to torque e-12
k_opto = 25*10^-12;

dydt=[y(2);
      (-beta*y(2)+y(3)+y(4))/I;
      (-y(3)+P(theta_bar_now-y(1))+(dot_theta_bar_now-y(2))*V(theta_bar_now-y(1)))/tau_bar;
      (-y(4)+(dot_theta_grating_now-y(2))*k_opto)/tau_opto];
end

function dydt=sim_ek1(t,y)
global t0 theta_bar dot_theta_bar dot_theta_grating;
theta_bar_now=interp1(t0,theta_bar,t);
dot_theta_bar_now=interp1(t0,dot_theta_bar,t);
dot_theta_grating_now=interp1(t0,dot_theta_grating,t);
I=6*10^-14;%moment of inertia
beta=10^-11;%drag coefficient
tau_bar=0.05;
tau_opto=0.04;
P=@(theta)+10e-11*sin(theta); % pos to torque e-11
V=@(theta)1e-12*(20-2*theta^2); % vel to torque e-12
k_opto = 25*10^-12;

dydt=[y(2);
      (-beta*y(2)+y(3)+y(4))/I;
      (-y(3)+P(theta_bar_now-y(1))+(dot_theta_bar_now-y(2))*V(theta_bar_now-y(1)))/tau_bar;
      (-y(4)+(dot_theta_grating_now-y(2)+y(5))*k_opto)/tau_opto;
      (-beta*y(5)+y(3))/I;];
end

function dydt=sim_all_or_non(t,y)
global t0 theta_bar dot_theta_bar dot_theta_grating supp_global
theta_bar_now=interp1(t0,theta_bar,t);
dot_theta_bar_now=interp1(t0,dot_theta_bar,t);
dot_theta_grating_now=interp1(t0,dot_theta_grating,t);
I=6*10^-14;%moment of inertia
beta=10^-11;%drag coefficient
tau_bar=0.05;
tau_opto=0.04;
P=@(theta)+10e-11*sin(theta); % pos to torque e-11
V=@(theta)1e-12*(20-2*theta^2); % vel to torque e-12
if abs(y(3))>=1.5e-12
    k_opto = 0;
    supp = 1e-10;
else
    k_opto = 25*10^-12;
    supp = 0;
end

supp_global = [supp_global; supp];

dydt=[y(2);
      (-beta*y(2)+y(3)+y(4))/I;
      (-y(3)+P(theta_bar_now-y(1))+(dot_theta_bar_now-y(2))*V(theta_bar_now-y(1)))/tau_bar;
      (-y(4)+(dot_theta_grating_now-y(2))*k_opto)/tau_opto];
end