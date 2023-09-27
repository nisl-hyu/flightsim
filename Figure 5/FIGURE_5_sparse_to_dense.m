% FIGURE 5
% Sparse to dense background
% 2023/09/27 Angel Canelo & Anmo Kim

clear all; clc; close all;
h = 0.0002; % step size
t0=0:h:2.2;
theta_grating = 0*linspace(-pi,pi/2,size(t0,2));
dot_theta_grating = diff(theta_grating([1 1:end]))/(t0(2)-t0(1));
ds = 1;     % For downsampling plots
% epplt1 = 200; epplt2 = 500; epplt3 = 800;
epplt1 = 200; epplt2 = 410; epplt3 = 500;

theta_bar = pi/3./(1+exp((-t0+0.5)*30));
dot_theta_bar = diff(theta_bar([1 1:end]))/(t0(2)-t0(1));

ep_vec = 0:0.01:20;
%k_opto = 0.5e-1+((7e-1)./(1+exp((-ep_vec+5)*2)));
k_opto = 0.5e-1*ones(1,length(ep_vec)); k_opto(ep_vec>4) = 7e-1;
n_ep = 0:length(ep_vec); c_e_plot = 0; base_e = 0;

% MLP weights
global wh wo ro_global
units = 5; inputs = 1; outputs = 1;
wh = 0.1*rand(units,inputs);
wo = 0.025*rand(outputs,units);
alpha = 0.00008; count = 0;
[t, y0] = ode45(@(t,y) sim_ideal(t,y,theta_bar,dot_theta_bar,t0),t0, [0 0 0]');
for ep=1:length(ep_vec)
    [t01, y01] = ode45(@(t,y) sim_effc(t,y,theta_bar,dot_theta_bar,...
        dot_theta_grating,y0(:,1),t0,k_opto(ep),alpha),t0, [0 0 0 0 0]');
    [t02, y02] = ode45(@(t,y) sim_fixed(t,y,theta_bar,dot_theta_bar,...
        dot_theta_grating,y0(:,1),t0,k_opto(ep)),t0, [0 0 0 0 0]');

    %%%%%%%%%%%% LATENCIES %%%%%%%%%%%%%%%%%%
    %theta_bar_now = interp1(t0,theta_bar,t);
    % Calculate the baseline level of the signal
    % baseline = mean(theta_bar);
    % % Set the threshold
    % %threshold = 0 + 0.5*std(theta_bar);
    % threshold = 0.5*max(theta_bar);
    % % Detect the onsets
    % onsets = find(theta_bar >= threshold,1,'first');
    t_onset = 0.4;%t0(onsets);
    % bar
    val50 = 0.5*max(y01(:,1));    % 50% of peak value
    ind_b = find(y01(:,1)>=val50,1,'first');
    t50bar = t(ind_b);
    latency_bar(ep) = abs(t50bar - t_onset);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    effc{ep} = y01(:,5)';
    ye1_buff{ep} = y01(:,1)';
    k_opto_buff(ep) = mean(ro_global); %mean(y01(:,6)');
    MSE(ep) = mean((180/pi*y0(:,1)'-180/pi*y02(:,1)').^2);
    MSE_exp(ep) = mean((180/pi*y0(:,1)'-180/pi*y01(:,1)').^2);
    ro_global = [];
    disp(['Iteration ',num2str(ep),' of ', num2str(length(ep_vec))]);
end
% Plotting
cm = hsv(4);
hfig=figure();clf;set(gcf,'Color','w'); haxes(1) = subplot(4,1,1); hold on
h2=plot(t01(1:ds:end), ye1_buff{epplt1}(1:ds:end)*180/pi,'Color',cm(1,:),'LineStyle','-');
h3=plot(t01(1:ds:end), ye1_buff{epplt2}(1:ds:end)*180/pi,'Color',cm(2,:));
h4=plot(t01(1:ds:end), ye1_buff{epplt3}(1:ds:end)*180/pi,'Color',cm(3,:));
h5=plot(t01(1:ds:end), ye1_buff{end}(1:ds:end)*180/pi,'Color',cm(4,:),'LineStyle','--');
h6=plot(t0(1:ds:end),theta_grating(1:ds:end)*180/pi,'Color',[0, 0, 0, 0.5]);
h1=plot(t0(1:ds:end),theta_bar(1:ds:end)*180/pi, 'black');
legend([h1(1),h6(1),h2(1),h3(1),h4(1),h5(1)],'bar','grating',['ep ',num2str(epplt1)], ['ep ',num2str(epplt2)], ['ep ',num2str(epplt3)], ['ep ',num2str(length(ye1_buff)-1)]); xlabel('Time (s)')
set(haxes(1), 'XLim',[0 2.2], 'XTick',0:0.2:2.2); ylabel('Body angle (degree)');
haxes(2) = subplot(4,1,2); plot(n_ep(1:end-1),k_opto,'black'); hold on;
plot(n_ep(1:end-1),k_opto_buff,'b')
legend('Vgrtng per ep','Vgrtng learned by effc'); ylabel('Vgrtng (Nm)')
haxes(3) = subplot(4,1,3); plot(n_ep(1:end-1),MSE, 'Color',[0.9290 0.6940 0.1250]); hold on;
plot(n_ep(1:end-1),MSE_exp, 'Color',[0.4940 0.1840 0.5560]);
legend('MSE fixed', 'MSE auto-tuned'); ylabel('MSE (degree)');
haxes(4) = subplot(4,1,4);
plot(n_ep(1:end-1),latency_bar,'b');
xlabel('Iteration#'); ylabel('50% onset latency (ms)');
set(haxes(:),'TickDir','out','Box','off');
%print(hfig,'-painters','-depsc', 'fig_5.eps')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dydt=sim_ideal(t,y, theta_bar, dot_theta_bar,t0)
theta_bar_now=interp1(t0,theta_bar,t);
dot_theta_bar_now=interp1(t0,dot_theta_bar,t);
beta=10^-1;%drag coefficient
I = 6*10^-4;%moment of inertia
tau_bar=0.05;
P=@(theta)+10e-1*sin(theta); % pos to torque e-11
V=@(theta)1e-2*(20-2*theta^2); % vel to torque e-12
dydt=[y(2);
      (-beta*y(2)+y(3))/I;
      (-y(3)+P(theta_bar_now-y(1))+(dot_theta_bar_now-y(2))*V(theta_bar_now-y(1)))/tau_bar];
end

function dydt=sim_effc(t,y, theta_bar, dot_theta_bar,dot_theta_grating,y0,t0,Vgrat, alpha)
global wh wo ro_global
theta_bar_now=interp1(t0,theta_bar,t);
dot_theta_bar_now=interp1(t0,dot_theta_bar,t);
dot_theta_grating_now=interp1(t0,dot_theta_grating,t);
y0_now=interp1(t0,y0,t);
beta=10^-1;%drag coefficient
I = 6*10^-4;%moment of inertia
tau_bar=0.05;
tau_opto=0.04;
P=@(theta)+10e-1*sin(theta); % pos to torque e-11
V=@(theta)1e-2*(20-2*theta^2); % vel to torque e-12

ri = y(3);
rh = 1./(1+exp(-wh*ri)); %wh*ri;    % applying sigmoid to hidden layer
ro = wo*rh;

do = (y0_now-y(1))^2;      % error
dh = (rh.*(1-rh)).*(wo'*do);    %wo'*do     % sigmoid in hidden layer
wo = wo + alpha*(do*rh');
wh = wh + alpha*(dh*ri);

ro_global = [ro_global; ro];

dydt=[y(2);
      (-beta*y(2)+y(3)+y(4))/I;
      (-y(3)+P(theta_bar_now-y(1))+(dot_theta_bar_now-y(2))*V(theta_bar_now-y(1)))/tau_bar;
      (-y(4)+(dot_theta_grating_now-y(2)+y(5))*Vgrat)/tau_opto;
      (-beta*y(5)+y(3)*ro)/I];
end

function dydt=sim_fixed(t,y, theta_bar, dot_theta_bar,dot_theta_grating,y0,t0,Vgrat)
theta_bar_now=interp1(t0,theta_bar,t);
dot_theta_bar_now=interp1(t0,dot_theta_bar,t);
dot_theta_grating_now=interp1(t0,dot_theta_grating,t);
y0_now=interp1(t0,y0,t);
beta=10^-1;%drag coefficient
I = 6*10^-4;%moment of inertia
tau_bar=0.05;
tau_opto=0.04;
P=@(theta)+10e-1*sin(theta); % pos to torque e-11
V=@(theta)1e-2*(20-2*theta^2); % vel to torque e-12

dydt=[y(2);
      (-beta*y(2)+y(3)+y(4))/I;
      (-y(3)+P(theta_bar_now-y(1))+(dot_theta_bar_now-y(2))*V(theta_bar_now-y(1)))/tau_bar;
      (-y(4)+(dot_theta_grating_now-y(2)+y(5))*Vgrat)/tau_opto;
      (-beta*y(5)+y(3)*0.5e-1)/I];
end