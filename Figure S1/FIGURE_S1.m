% FIGURE S1
% 2023/07/21  Angel Canelo & Anmo Kim

clear all; clc; close all;

global t0 theta_bar dot_theta_bar alpha tau_fatigue theta_bar_ccw dot_theta_bar_ccw;
t0=-1:0.001:11;
alpha_vec = linspace(0.4,1,7);
tau_fatigue_vec = linspace(11000,17000,7);
theta_bar = linspace(-pi,pi,size(t0,2));
dot_theta_bar = diff(theta_bar([1 1:end]))/(t0(2)-t0(1));

theta_bar_ccw = linspace(pi,-pi,size(t0,2));
dot_theta_bar_ccw = diff(theta_bar_ccw([1 1:end]))/(t0(2)-t0(1));
dir = [1 -1];

%%%%%%%%%%%%
pos_fcn_bar = 10*sin(theta_bar);%position function
vel_fcn_bar = (20-2*theta_bar.^2);
A = -5; x = theta_bar; p = 2*pi;
pos_fcn_spot = A*4/p*abs(mod(x-p/4,p)-p/2)-A;
vel_fcn_spot = 5./(1+exp(-2*(3/8*pi-abs(theta_bar))));%velocity function

pos_fcn_grat = zeros(1,length(theta_bar));%position function
vel_fcn_grat = 0*ones(1,length(theta_bar));%velocity function
tau_rise = 500;
input = 25*ones(1,length(theta_bar)); input(t0>=10) = 0;
    for nn=1:length(theta_bar)-1
        if t0(nn)>=0  % onset time
            vel_fcn_grat(nn+1) = vel_fcn_grat(nn) + 1/tau_rise*(-vel_fcn_grat(nn)+input(nn));
        end
    end
%%%%%%%%%%%%

clc;
load('..\data\POSVEL_data\Averaged_wba');
fill = 1;
ind = [3 2 1];
colr = {[0 1 0] [0 0 1] [1 0 0]};
for kkk=1:length(alpha_vec)
    alpha = alpha_vec(kkk);
for nnn=1:length(tau_fatigue_vec)
    tau_fatigue = tau_fatigue_vec(nnn);
for select=1:3
dt=0.001;
[t, torque0]=simFly1(theta_bar, dt, select, dir(1));
[t_ccw, torque0_ccw]=simFly1(theta_bar_ccw, dt, select, dir(2));

%%%%%%%%%%%% PERFORMANCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%torque_exp = 135/5*pi/180*wba_full{select};
torque_exp = exp_wba{select};
t_exp = linspace(-1,11,size(torque_exp,2));
torque_exp_now = interp1(t_exp,torque_exp,t);
    if fill<=3
        torque_exp_save(fill,:) = interp1(t_exp,torque_exp,t);
        fill = fill + 1;
    end
err = sqrt(mean((torque_exp_now - torque0).^2)); %RMSE
if select == 1
    spot_torque_save(kkk,nnn,:) = torque0;
    spot_torque_save_ccw(kkk,nnn,:) = torque0_ccw;
    spot_opt_param(kkk,nnn) = err;
elseif select == 2
    bar_torque_save(kkk,nnn,:) = torque0;
    bar_torque_save_ccw(kkk,nnn,:) = torque0_ccw;
    bar_opt_param(kkk,nnn) = err;
else % grating is computed in a reduced range
    opt_param(kkk,nnn) = err;
    %torque_exp_save = torque_exp_now;
    torque_save(kkk,nnn,:) = torque0;
    torque_save_ccw(kkk,nnn,:) = torque0_ccw;
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
end
end
end
%plotting

[M,I] = min(opt_param(:));
[row,col]=ind2sub(size(opt_param),I);

% Grating
figure();
subplot(2,1,1)
load('..\data\POSVEL_data\position_velocity_components_grating');
plot(t,squeeze((torque_save(row,col,:)+flip(torque_save_ccw(row,col,:)))/2),'Color','red','LineWidth',2); hold on;
plot(t, pos_fcn_grat,'--','Color','red','LineWidth',2);
legend('fatigue','non fatigue')
ylabel('Position response (degree)')
title('Components to grating')
subplot(2,1,2)
plot(t,squeeze((torque_save(row,col,:)-flip(torque_save_ccw(row,col,:)))/2),'Color','red','LineWidth',2); hold on
plot(t, vel_fcn_grat,'--','Color','red','LineWidth',2);
legend('fatigue','non fatigue')
ylabel('Velocity response (degree)')
xlabel('Time (s)')

% Bar
figure();
subplot(2,1,1)
load('..\data\POSVEL_data\position_velocity_components_bar');
t_exp = linspace(0,10,243);
bar_pos = interp1(t_exp,pos_tm,t);
bar_vel = interp1(t_exp,vel_tm,t);
plot(t,squeeze((bar_torque_save(row,col,:)-bar_torque_save_ccw(row,col,:))/2),'Color','blue','LineWidth',2); hold on;
plot(t, pos_fcn_bar,'--','Color','blue','LineWidth',2);
legend('fatigue','non fatigue')
ylabel('Position response (degree)')
title('Components to bar')
subplot(2,1,2)
plot(t,squeeze((bar_torque_save(row,col,:)+bar_torque_save_ccw(row,col,:))/2),'Color','blue','LineWidth',2); hold on;
plot(t, vel_fcn_bar,'--','Color','blue','LineWidth',2);
legend('fatigue','non fatigue')
ylabel('Velocity response (degree)')
xlabel('Time (s)')

% Spot
figure();
subplot(2,1,1)
load('..\data\POSVEL_data\position_velocity_components_spot');
spot_pos = interp1(t_exp,pos_tm,t);
spot_vel = interp1(t_exp,vel_tm,t);
plot(t,squeeze((spot_torque_save(row,col,:)-spot_torque_save_ccw(row,col,:))/2),'Color','green','LineWidth',2); hold on;
plot(t, pos_fcn_spot,'--','Color','green','LineWidth',2);
legend('fatigue','non fatigue')
ylabel('Position response (degree)')
title('Components to spot')
subplot(2,1,2)
plot(t,squeeze((spot_torque_save(row,col,:)+spot_torque_save_ccw(row,col,:))/2),'Color','green','LineWidth',2); hold on;
plot(t, vel_fcn_spot,'--','Color','green','LineWidth',2);
legend('fatigue','non fatigue')
ylabel('Velocity response (degree)')
xlabel('Time (s)')

figure();
subplot(3,1,1)
plot(t, squeeze(torque_save(row,col,:)),'Color','red','LineWidth',2); hold on
plot(t, vel_fcn_grat,'--','Color','red','LineWidth',2);
plot(t, torque_exp_save(3,:),'Color','black','LineWidth',2);
legend('fatigue','non fatigue')
ylabel('Wcw (degree)')
title('Wing response to each pattern')
subplot(3,1,2)
plot(t, squeeze(bar_torque_save(row,col,:)),'Color','blue','LineWidth',2); hold on
plot(t, pos_fcn_bar+vel_fcn_bar,'--','Color','blue','LineWidth',2);
plot(t, torque_exp_save(2,:),'Color','black','LineWidth',2);
legend('fatigue','non fatigue')
ylabel('Wcw (degree)')
subplot(3,1,3)
plot(t, squeeze(spot_torque_save(row,col,:)),'Color','green','LineWidth',2); hold on
plot(t, pos_fcn_spot+vel_fcn_spot,'--','Color','green','LineWidth',2);
plot(t, torque_exp_save(1,:),'Color','black','LineWidth',2);
legend('fatigue','non fatigue')
ylabel('Wcw (degree)')
xlabel('Time (s)')

xvalues = tau_fatigue_vec;
yvalues = alpha_vec;
fg = figure();
h = heatmap(xvalues,yvalues,round(opt_param,3), 'Colormap',jet); colorbar
h.Title = 'RMSE (deg) experimental vs predicted torque';
h.XLabel = 'tau fatigue (ms)';
h.YLabel = 'alpha (a.u.)';
%print(fg,'-painters','-depsc', 'matrix_tau_alpha.eps')

function [t0, output]=simFly1(psi, dt, select,dir)
global dot_theta_bar tau_fatigue alpha
%psi: pattern position in radian
%output: torque that a fly is generating 
t0=-1:0.001:11;
if select == 2
%%%%%%%%%%%%%%% Bar %%%%%%%%%%%%%%%%%%%%%%%
pos_fcn = 10*sin(psi);%position function
vel_fcn = 20-2*psi.^2;%velocity function
output2=pos_fcn+vel_fcn;
output = 0*ones(1,length(psi));
A=zeros(1,length(psi));
tau_rise = 500;
    for nn=1:length(psi)-1
            output(nn+1) = output(nn) + 1/tau_rise*(-output(nn)+output2(nn)-A(nn));
            A(nn+1) = A(nn) + 1/tau_fatigue*(-A(nn)+alpha*output2(nn));
    end
%%%%%%%%%%%%%%% Spot %%%%%%%%%%%%%%%%%%%%%%%
elseif select == 1
%pos_fcn = -8 * sawtooth(1*(psi+pi/2),1/2);%position function
A = -5; x = psi; p = 2*pi;
pos_fcn = A*4/p*abs(mod(x-p/4,p)-p/2)-A;
vel_fcn = 5./(1+exp(-2*(3/8*pi-abs(psi)))); %10-1*psi.^2;%velocity function
output2=pos_fcn+vel_fcn;
output = 0*ones(1,length(psi));
A=zeros(1,length(psi));
tau_rise = 500;
    for nn=1:length(psi)-1
            output(nn+1) = output(nn) + 1/tau_rise*(-output(nn)+output2(nn)-A(nn));
            A(nn+1) = A(nn) + 1/tau_fatigue*(-A(nn)+alpha*output2(nn));
    end
%%%%%%%%%%%%%%% Grating %%%%%%%%%%%%%%%%%%%
elseif select == 3
pos_fcn = zeros(size(psi,1));%position function
vel_fcn = 0*ones(1,length(psi));%velocity function
tau_rise = 500;
A=zeros(1,length(psi));
input = dir*25*ones(1,length(psi)); input(t0>=10) = 0;
    for nn=1:length(psi)-1
        if t0(nn)>=0  % onset time
            vel_fcn(nn+1) = vel_fcn(nn) + 1/tau_rise*(-vel_fcn(nn)+input(nn)-A(nn));
            A(nn+1) = A(nn) + 1/tau_fatigue*(-A(nn)+alpha*vel_fcn(nn));
        end
        %A(nn+1) = A(nn) + 1/5000*(-A(nn)+0.4*70);
    end
    %figure; plot(A)
    output=pos_fcn+(vel_fcn);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
end