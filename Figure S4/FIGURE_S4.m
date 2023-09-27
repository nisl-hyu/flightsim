% FIGURE S4

% 2023/09/27 Angel Canelo & Anmo Kim

clear all; close all; clc;
global t0 theta_bar dot_theta_bar dot_theta_grating F_ext;
t0=0:0.0001:5;
tonset1 = 0.5; tonset2 = 1.5; tonset3 = 2.5; tonset4 = 3.5;
for kk = 1:1
    if kk==1      
        theta_bar = 60*pi/180./(1+exp((-t0+tonset1)*30))-60*pi/180./(1+exp((-t0+tonset2)*30))+60*pi/180./(1+exp((-t0+tonset3)*30));%-60*pi/180./(1+exp((-t0+2.25)*30));  % for sigmoid-like stimulus
        theta_grating = pi/3-pi/3./(1+exp((-t0+tonset2)*30))-pi/3./(1+exp((-t0+tonset3)*30))+pi/3./(1+exp((-t0+tonset4)*30));  %30deg righwards
        F_ext = zeros(1,length(t0));
    elseif kk==2
        theta_bar = 60*pi/180./(1+exp((-t0+1.25)*30))-60*pi/180./(1+exp((-t0+2.25)*30));  % for sigmoid-like stimulus
        theta_grating = 0*-60*pi/180./(1+exp((-t0+0.25)*30));%+30*pi/180./(1+exp((-t0+0.5)*30));  %30deg righwards
        F_ext = zeros(1,length(t0)); F_ext(t0>=0.25 & t0<=0.30) = -3e-10;
    end
dot_theta_bar = diff(theta_bar([1 1:end]))/(t0(2)-t0(1));
dot_theta_grating = diff(theta_grating([1 1:end]))/(t0(2)-t0(1));

[t2, y2]=ode45(@sim_ek1,[0 t0(end)], [0 0 0 0 0]');
[t4, y4]=ode45(@sim_allnone,[0 t0(end)], [0 0 0 0]');
[t1, y1]=ode45(@sim_bar1,[0 t0(end)], [0 0 0]');
[t3, y3]=ode45(@sim_grat1,[0 t0(end)], [0 0 0]');

%%% Latency
%threshold = 0.5*max(theta_bar);
% Detect the onsets
%onsets = find(theta_bar >= threshold,1,'first');
t_onset = tonset1-0.1; %t0(onsets);
meanpreEC1 = mean(y2(t2>=0.1 & t2<=0.4, 1));
meanpostEC1 = mean(y2(t2>=0.6 & t2<=0.9, 1));
amplEC1 = abs(meanpostEC1-meanpreEC1);
amplEC11 = 30*pi/180;%0.5*max(y2(t2>=0.4 & t2<=0.7, 1));    % 50% of peak value
ind_b = find(y2(t2>=0.4 & t2<=0.9, 1)>=amplEC11,1,'first');
ff = find(t2>=0.4,1,"first");
t50bar = t2(ff+ind_b);
latency_ec1 = abs(t50bar - t_onset);

meanpreallnone1 = mean(y4(t4>=0.1 & t4<=0.4, 1));
meanpostallnone1 = mean(y4(t4>=0.6 & t4<=0.9, 1));
amplallnone1 = abs(meanpostallnone1-meanpreallnone1);
amplallnone11 = 30*pi/180;%0.5*max(y4(t4>=0.4 & t4<=0.7, 1));    % 50% of peak value
ind_b = find(y4(t4>=0.4 & t4<=0.9, 1)>=amplallnone11,1,'first');
ff = find(t4>=0.4,1,"first");
t50bar = t4(ff+ind_b);
latency_allnone1 = abs(t50bar - t_onset);

t_onset2 = tonset2-0.1;
meanpreEC2 = mean(y2(t2>=tonset2-0.4 & t2<=tonset2-0.1, 1));
meanpostEC2 = mean(y2(t2>=tonset2+0.1 & t2<=tonset2+0.4, 1));
amplEC2 = abs(meanpostEC2-meanpreEC2);
amplEC22 = 30*pi/180;%1.5*min(y2(t2>=1.4 & t2<=1.7, 1));    % 50% of peak value
ind_b = find(y2(t2>=tonset2-0.1 & t2<=tonset2+0.2, 1)<=amplEC22,1,'first');
ff = find(t2>=tonset2-0.1,1,"first");
t50bar = t2(ff+ind_b);
latency_ec2 = abs(t50bar - t_onset2);

meanpreallnone2 = mean(y4(t4>=tonset2-0.4 & t4<=tonset2-0.1, 1));
meanpostallnone2 = mean(y4(t4>=tonset2+0.1 & t4<=tonset2+0.4, 1));
amplallnone2 = abs(meanpostallnone2-meanpreallnone2);
amplallnone22 = 30*pi/180;%1.5*min(y4(t4>=1.4 & t4<=1.7, 1));    % 50% of peak value
ind_b = find(y4(t4>=tonset2-0.1 & t4<=tonset2+0.2, 1)<=amplallnone22,1,'first');
ff = find(t4>=tonset2-0.1,1,"first");
t50bar = t4(ff+ind_b);
latency_allnone2 = abs(t50bar - t_onset2);

t_onset3 = tonset3-0.1;
meanpreEC3 = mean(y2(t2>=tonset3-0.4 & t2<=tonset3-0.1, 1));
meanpostEC3 = mean(y2(t2>=tonset3+0.1 & t2<=tonset3+0.4, 1));
amplEC3 = abs(meanpostEC3-meanpreEC3);
amplEC33 = 30*pi/180;%0.5*max(y2(t2>=2.4 & t2<=2.7, 1));    % 50% of peak value
ind_b = find(y2(t2>=tonset3-0.1 & t2<=tonset3+0.2, 1)>=amplEC33,1,'first');
ff = find(t2>=tonset3-0.1,1,"first");
t50bar = t2(ff+ind_b);
latency_ec3 = abs(t50bar - t_onset3);

meanpreallnone3 = mean(y4(t4>=tonset3-0.4 & t4<=tonset3-0.1, 1));
meanpostallnone3 = mean(y4(t4>=tonset3+0.1 & t4<=tonset3+0.4, 1));
amplallnone3 = abs(meanpostallnone3-meanpreallnone3);
amplallnone33 = 30*pi/180;%0.5*max(y4(t4>=2.4 & t4<=2.7, 1));    % 50% of peak value
ind_b = find(y4(t4>=tonset3-0.1 & t4<=tonset3+0.4, 1)>=amplallnone33,1,'first');
ff = find(t4>=tonset3-0.1,1,"first");
t50bar = t4(ff+ind_b);
latency_allnone3 = abs(t50bar - t_onset3);
%%%%%%%%%%%%%

%plotting
hfig=figure(kk);clf;set(gcf,'Color','w');
haxes(1)=subplot(211);
plot(t2,y2(:,1)/pi*180, 'Color',[0.9290 0.6940 0.1250],'LineWidth',1.2); hold on
plot(t4,y4(:,1)/pi*180, 'Color',[0.6350 0.0780 0.1840],'LineWidth',1.2);
% plot(t1,y1(:,1)/pi*180,'--b', 'LineWidth',1.2);
% plot(t3,y3(:,1)/pi*180,'red','LineWidth',1.2);
plot(t0,theta_bar/pi*180,'black','LineWidth',1.2)
plot(t0,theta_grating/pi*180,'Color',[.7 .7 .7],'LineWidth',1.2);
legend('graded','all-or-none'); set(0,'DefaultLegendAutoUpdate','off')
title('fly angular position');ylabel('Body angle(deg)'); xlabel('time (s)'); set(haxes(1),'YLim',[-60 180]);
haxes(2)=subplot(212); hold on;
plot(t2,y2(:,3),'Color',[0.9290 0.6940 0.1250],'LineWidth',1.2);
plot(t2,y2(:,4),'Color',[0.9290 0.6940 0.1250],'LineStyle','--','LineWidth',1.2);
plot(t4,y4(:,3), 'Color',[0.6350 0.0780 0.1840],'LineWidth',1.2);
plot(t4,y4(:,4), 'Color',[0.6350 0.0780 0.1840],'LineStyle','--','LineWidth',1.2);
plot(t0,F_ext,'Color', 'black','LineWidth',1.2);
legend('torque to bar (graded)','torque to grating (graded)','torque to bar (all-or-none)','torque to grating (all-or-none)');
title('Torque response');xlabel('time (s)'); ylabel('(Nm)'); hold on
set(haxes(1),'Ylim', [-60 100]);
set(haxes(2),'Ylim', [-1e-10 1e-10]);
set(haxes,'Box','off','TickDir','out');
set(haxes(1:1),'XTickLabel',[]);
for i=1:length(haxes)
    ylim=get(haxes(i),'YLim');
    hpatch=patch(haxes(i),[0.4 0.6 0.6 0.4],ylim([1 1 2 2]),'r','FaceColor',ones(1,3)*.9,'EdgeColor','none');
    hpatch2=patch(haxes(i),[1.4 1.6 1.6 1.4],ylim([1 1 2 2]),'r','FaceColor',ones(1,3)*.9,'EdgeColor','none');
    hpatch3=patch(haxes(i),[2.4 2.6 2.6 2.4],ylim([1 1 2 2]),'r','FaceColor',ones(1,3)*.9,'EdgeColor','none');
    hpatch4=patch(haxes(i),[3.4 3.6 3.6 3.4],ylim([1 1 2 2]),'r','FaceColor',ones(1,3)*.9,'EdgeColor','none');
    uistack(hpatch,'bottom');
    uistack(hpatch2,'bottom');
    uistack(hpatch3,'bottom');
    uistack(hpatch4,'bottom');
end

fib = figure();clf;set(gcf,'Color','w'); barcolor = [0.9290 0.6940 0.1250; 0.6350 0.0780 0.1840]; % Blue, Green, Red
figure(fib); subplot(6,1,1);
plt_1 = bar([amplEC1*180/pi amplallnone1*180/pi], 'FaceColor', 'flat'); plt_1.CData = barcolor;
subplot(6,1,2);
plt_2 = bar([latency_ec1 latency_allnone1], 'FaceColor', 'flat'); plt_2.CData = barcolor;
subplot(6,1,3);
plt_3 = bar([amplEC2*180/pi amplallnone2*180/pi], 'FaceColor', 'flat'); plt_3.CData = barcolor;
subplot(6,1,4);
plt_4 = bar([latency_ec2 latency_allnone2], 'FaceColor', 'flat'); plt_4.CData = barcolor;
subplot(6,1,5);
plt_5 = bar([amplEC3*180/pi amplallnone3*180/pi], 'FaceColor', 'flat'); plt_5.CData = barcolor;
subplot(6,1,6);
plt_6 = bar([latency_ec3 latency_allnone3], 'FaceColor', 'flat'); plt_6.CData = barcolor;

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

function dydt=sim_ek1(t,y)
global t0 theta_bar dot_theta_bar dot_theta_grating F_ext;
theta_bar_now=interp1(t0,theta_bar,t);
dot_theta_bar_now=interp1(t0,dot_theta_bar,t);
dot_theta_grating_now=interp1(t0,dot_theta_grating,t);
F_ext_now=interp1(t0,F_ext,t);
I=6*10^-14;%moment of inertia
beta=10^-11;%drag coefficient
tau_bar=0.05;
tau_opto=0.04;
P=@(theta)+10e-11*sin(theta); % pos to torque e-11
V=@(theta)1e-12*(20-2*theta^2); % vel to torque e-12
k_opto = 25*10^-12;

dydt=[y(2);
      (-beta*y(2)+y(3)+y(4)+F_ext_now)/I;
      (-y(3)+P(theta_bar_now-y(1))+(dot_theta_bar_now-y(2))*V(theta_bar_now-y(1)))/tau_bar;
      (-y(4)+(dot_theta_grating_now+y(5)-y(2))*k_opto)/tau_opto;
      (-beta*y(5)+y(3))/I;];
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

function dydt=sim_allnone(t,y)
global t0 theta_bar dot_theta_bar dot_theta_grating F_ext;
theta_bar_now=interp1(t0,theta_bar,t);
dot_theta_bar_now=interp1(t0,dot_theta_bar,t);
dot_theta_grating_now=interp1(t0,dot_theta_grating,t);
F_ext_now=interp1(t0,F_ext,t);
I=6*10^-14;%moment of inertia
beta=10^-11;%drag coefficient
tau_bar=0.05;
tau_opto=0.04;
P=@(theta)+10e-11*sin(theta); % pos to torque e-11
V=@(theta)1e-12*(20-2*theta^2); % vel to torque e-12
k_opto = 25*10^-12;
if abs(y(3))>1.5e-11
    k_opto = 0;
else
    k_opto = 25*10^-12;
end

dydt=[y(2);
      (-beta*y(2)+y(3)+y(4)+F_ext_now)/I;
      (-y(3)+P(theta_bar_now-y(1))+(dot_theta_bar_now-y(2))*V(theta_bar_now-y(1)))/tau_bar;
      (-y(4)+(dot_theta_grating_now-y(2))*k_opto)/tau_opto];
end