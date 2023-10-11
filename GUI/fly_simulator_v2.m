function varargout = fly_simulator_v2(varargin)
% GUI
% 2023/07/26 Angel Canelo & Anmo Kim

% FLY_SIMULATOR_V2 MATLAB code for fly_simulator_v2.fig
%      FLY_SIMULATOR_V2, by itself, creates a new FLY_SIMULATOR_V2 or raises the existing
%      singleton*.
%
%      H = FLY_SIMULATOR_V2 returns the handle to a new FLY_SIMULATOR_V2 or the handle to
%      the existing singleton*.
%
%      FLY_SIMULATOR_V2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FLY_SIMULATOR_V2.M with the given input arguments.
%
%      FLY_SIMULATOR_V2('Property','Value',...) creates a new FLY_SIMULATOR_V2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before fly_simulator_v2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to fly_simulator_v2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help fly_simulator_v2

% Last Modified by GUIDE v2.5 10-Oct-2023 11:01:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fly_simulator_v2_OpeningFcn, ...
                   'gui_OutputFcn',  @fly_simulator_v2_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before fly_simulator_v2 is made visible.
function fly_simulator_v2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to fly_simulator_v2 (see VARARGIN)
clc;
% Choose default command line output for fly_simulator_v2
handles.output = hObject;
axes(handles.simulator);
% Update handles structure
guidata(hObject, handles);
set(handles.simulator,'visible','off')
axis(handles.simulator, 'equal'); axis([-3 3 -3 3]); P1 = [0 0];
P4_tr = viscircles(P1,2,'color','green','LineWidth', 20); hold on
line1 = line([-2.5 -2.3],[0 0], 'color','black','LineWidth', 2);
line2 = line([2.5 2.3],[0 0], 'color','black','LineWidth', 2);
line3 = line([0 0],[-2.5 -2.3], 'color','black','LineWidth', 2);
line4 = line([0 0],[2.5 2.3], 'color','black','LineWidth', 2);
% saveas(gca, './gui_figures/simulator_gui.png', 'png');
set(handles.pushbutton1,'Enable', 'off')
set(handles.ampl,'Enable', 'off')
set(handles.ampl2,'Enable', 'off')
set(handles.onset,'Enable', 'off')
% UIWAIT makes fly_simulator wait for user response (see UIRESUME)
% uiwait(handles.figure1);



% --- Outputs from this function are returned to the command line.
function varargout = fly_simulator_v2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 %load('./generated data/gui_variables.mat');
 % open video file
 set(handles.popupmenu1,'Enable', 'off')
 set(handles.pushbutton1,'Enable', 'off')
 set(handles.ampl,'Enable', 'off')
 set(handles.ampl2,'Enable', 'off')
 set(handles.onset,'Enable', 'off')
 axes(handles.simulator);
 global theta_fly theta_obj theta_fly2 theta_obj2 theta_grat amp amp2
    str=get(handles.popupmenu1,'String');
    val=get(handles.popupmenu1,'Value');
    patt=str{val};
    amp = str2double(get(handles.ampl, 'String'));
    amp2 = str2double(get(handles.ampl2, 'String'));
    if get(handles.popupmenu2,'Value')==2
        fly_simulator = VideoWriter(sprintf('./gui_figures/simulator_%s.mp4', patt),'MPEG-4'); %open video file
        fly_simulator.FrameRate = 20;  %can adjust this, 5 - 10 works well for me
        open(fly_simulator)
    end
    P1 = [0 0]; axis(handles.simulator, 'equal'); axis([-3 3 -3 3]);
    offset = pi/2; %theta_grat = zeros(length(theta_obj));
    onset = str2double(get(handles.onset, 'String'));
    t_onset = onset+0.1;
    if strcmp(patt, 'Grating') || strcmp(patt, 'Addition-based') || strcmp(patt,'Graded vs All-or-none') || strcmp(patt,'All-or-none efference copy')
        h = 0.001; t0=0:h:t_onset+1; check = 0 : 10 : length(t0);
    else
        h = 0.01; t0=0:h:t_onset+1; check = 0 : 1 : length(t0);
    end
    theta_obj = pi/180*amp./(1+exp((-t0+t_onset)*30)); dot_theta_bar = diff(theta_obj([1 1:end]))/(t0(2)-t0(1));
    theta_grating = pi/180*amp2./(1+exp((-t0+t_onset)*30));
    theta_grat = theta_grating;
    dot_theta_grating = diff(theta_grating([1 1:end]))/(t0(2)-t0(1));
    theta_fly(1) = 0; ye2(1) = 0; ye3(1) = 0; ye4(1) = 0; ye5(1) = 0;
    theta_fly2(1) = 0; yed2(1) = 0; yed3(1) = 0; yed4(1) = 0;
    beta=10^-11; I = 6*10^-14; tau_bar=0.05; tau_opto=0.04; k_opto = 25*10^-12; %5e-11;

    P=@(theta)+10e-11*sin(theta); % pos to torque e-11
    V=@(theta)1e-12*(20-2*theta^2); % vel to torque e-12

    A = -5e-11; p = 2*pi;
    Ps=@(theta)A*4/p*abs(mod(theta-p/4,p)-p/2)-A;
    Vs=@(theta)5e-12./(1+exp(-2*(3/8*pi-abs(theta))));
    ye3k1 = 0; ye3k2 = 0;ye3k3 = 0;
for th=1:length(theta_obj)
    textLabel = sprintf('Time = %0.2f s', t0(th));
    set(handles.time, 'String', textLabel);
    switch patt
        case 'Bar'
            theta_fly(th+1) = theta_fly(th) + h*ye2(th);
            ye2(th+1) = ye2(th) + h*(-beta*ye2(th)+ye3(th))/I;
            ye3(th+1) = ye3(th) + h*(-ye3(th)+P(theta_obj(th)-theta_fly(th))+(dot_theta_bar(th)-ye2(th))*V(theta_obj(th)-theta_fly(th)))/tau_bar;
        case 'Spot'
            theta_fly(th+1) = theta_fly(th) + h*ye2(th);
            ye2(th+1) = ye2(th) + h*(-beta*ye2(th)+ye3(th))/I;
            ye3(th+1) = ye3(th) + h*(-ye3(th)+Ps(theta_obj(th)-theta_fly(th))+(dot_theta_bar(th)-ye2(th))*Vs(theta_obj(th)-theta_fly(th)))/tau_bar;
        case 'Grating'
            theta_fly(th+1) = theta_fly(th) + h*ye2(th);
            ye2(th+1) = ye2(th) + h*(-beta*ye2(th)+ye3(th))/I;
            ye3(th+1) = ye3(th) + h*(-ye3(th)+(dot_theta_grating(th)-ye2(th))*k_opto)/tau_opto;
        case 'Addition-based'
            theta_fly(th+1) = theta_fly(th) + h*ye2(th);
            ye2(th+1) = ye2(th) + h*(-beta*ye2(th)+ye3(th)+ye4(th))/I;
            ye3(th+1) = ye3(th) + h*(-ye3(th)+P(theta_obj(th)-theta_fly(th))+(dot_theta_bar(th)-ye2(th))*V(theta_obj(th)-theta_fly(th)))/tau_bar;
            ye4(th+1) = ye4(th) + h*(-ye4(th)+(dot_theta_grating(th)-ye2(th))*k_opto)/tau_opto;
        case 'Graded efference copy'
            theta_fly(th+1) = theta_fly(th) + h*ye2(th);
            ye2(th+1) = ye2(th) + h*(-beta*ye2(th)+ye3(th)+ye4(th))/I;
            ye3(th+1) = ye3(th) + h*(-ye3(th)+P(theta_obj(th)-theta_fly(th))+(dot_theta_bar(th)-ye2(th))*V(theta_obj(th)-theta_fly(th)))/tau_bar;
            ye4(th+1) = ye4(th) + h*(-ye4(th)+(dot_theta_grating(th)-ye2(th))*k_opto+ye5(th))/tau_opto;
            ye5(th+1) = ye5(th) + h*(-beta*ye5(th)+ye3(th)*k_opto)/I;
        case 'All-or-none efference copy'
            theta_fly(th+1) = theta_fly(th) + h*ye2(th);
            ye2(th+1) = ye2(th) + h*(-beta*ye2(th)+ye3(th)+ye4(th))/I;
            ye3(th+1) = ye3(th) + h*(-ye3(th)+P(theta_obj(th)-theta_fly(th))+(dot_theta_bar(th)-ye2(th))*V(theta_obj(th)-theta_fly(th)))/tau_bar;
            if abs(ye3(th+1))>=1.5e-11
                k_opto = 0;
            else
                k_opto = 25*10^-12;
            end
            ye4(th+1) = ye4(th) + h*(-ye4(th)+(dot_theta_grating(th)-ye2(th))*k_opto)/tau_opto;
        case 'Graded vs All-or-none'
            theta_fly(th+1) = theta_fly(th) + h*ye2(th);
            ye2(th+1) = ye2(th) + h*(-beta*ye2(th)+ye3(th)+ye4(th))/I;
            ye3(th+1) = ye3(th) + h*(-ye3(th)+P(theta_obj(th)-theta_fly(th))+(dot_theta_bar(th)-ye2(th))*V(theta_obj(th)-theta_fly(th)))/tau_bar;
            ye4(th+1) = ye4(th) + h*(-ye4(th)+(dot_theta_grating(th)-ye2(th))*k_opto+ye5(th))/tau_opto;
            ye5(th+1) = ye5(th) + h*(-beta*ye5(th)+ye3(th)*k_opto)/I;
            
            theta_fly2(th+1) = theta_fly2(th) + h*yed2(th);
            yed2(th+1) = yed2(th) + h*(-beta*yed2(th)+yed3(th)+yed4(th))/I;
            yed3(th+1) = yed3(th) + h*(-yed3(th)+P(theta_obj(th)-theta_fly2(th))+(dot_theta_bar(th)-yed2(th))*V(theta_obj(th)-theta_fly2(th)))/tau_bar;
            if abs(yed3(th+1))>=1.5e-11
                k_opto2 = 0;
            else
                k_opto2 = 25*10^-12;
            end
            yed4(th+1) = yed4(th) + h*(-yed4(th)+(dot_theta_grating(th)-yed2(th))*k_opto2)/tau_opto;
    end        
    P2 = 0.8*[cos(offset-theta_fly(th)) sin(offset-theta_fly(th))];
    P4_tr = viscircles(P1,2,'color','green','LineWidth', 20);
    if isequal(patt,'Bar')       
        P4 = 2.25*[cos(offset-theta_obj(th)) sin(offset-theta_obj(th))];
        P5 = 1.75*[cos(offset-theta_obj(th)) sin(offset-theta_obj(th))];
        P6 = 0.4*[cos(offset-theta_fly(th)) sin(offset-theta_fly(th))];
        obj = line([P1(1) P4(1)],[P1(2) P4(2)], 'color','black','LineWidth', 10);
        hide = line([P1(1) P5(1)],[P1(2) P5(2)], 'color','white','LineWidth', 10);
        %fly = line([P1(1) P2(1)],[P1(2) P2(2)],'color','blue','LineWidth', 5);
        %P1_circ = viscircles(P1,0.07,'color','black');
        %P2_circ = viscircles(P2,0.1);
        %P2_circ = viscircles(P6,0.07,'color','black'); hold on
%         t = linspace(0,2*pi) ;
%         a = 0.5 ; b = 0.1 ; x = a*cos(t); y = b*sin(t) ;
%         R = Rot2D(pi/2-theta_fly(th)); XY_rotated = R*[x;y];
%         ellip = plot(XY_rotated(1,:),XY_rotated(2,:),'blue');
        hold on
        [fly1_handles fly1_handles2]=plot_fly(handles.simulator,theta_fly(th),ye3(th),[0 0 1],1.5);
        if th==1
            % saveas(gca, sprintf('./gui_figures/simulator_gui_%s.png', patt), 'png');
        end
        axes(handles.trace);
        plot(t0(1:th),theta_fly(1:th)*180/pi,'blue'); hold on; plot(t0(1:th),theta_obj(1:th)*180/pi,'black'); box off;
        xlim(handles.trace,[t0(1),t0(end)]); ylim(handles.trace,[min([0 amp]),max([0 amp])]); set(gca,'XTickLabel',[]);
        title('Bar      {\color[rgb]{0, 0, 1}Fly}');
        ylabel('Angle (deg)');
        axes(handles.torque);
        % plot(t0(1:th),ye3(1:th),'blue');
        plot(t0(1:th),ye3(1:th),'red');
        xlim(handles.torque,[t0(1),t0(end)]); ylim(handles.torque,[min([0 amp/abs(amp)*7e-11 ye3(1:th)]),max([0 amp/abs(amp)*7e-11 ye3(1:th)])]); box off;
        title('{\color[rgb]{1, 0, 0}Wings}');
        xlabel('Time (s)'); ylabel('Torque (Nm)');
        axes(handles.simulator);
        if get(handles.popupmenu2,'Value')==2
            frame = getframe(gcf); %get frame
            writeVideo(fly_simulator, frame);
        end
        pause(0.001); delete(obj);%delete(P1_circ);delete(P2_circ); %delete(fly)
        delete(hide); delete(fly1_handles);delete(fly1_handles2);%delete(ellip);
    elseif isequal(patt,'Spot')       
        P4 = 2.10*[cos(offset-theta_obj(th)) sin(offset-theta_obj(th))];
        P5 = 1.73*[cos(offset-theta_obj(th)) sin(offset-theta_obj(th))];
        P6 = 1.90*[cos(offset-theta_obj(th)) sin(offset-theta_obj(th))];
        P7 = 0.4*[cos(offset-theta_fly(th)) sin(offset-theta_fly(th))];
        obj = line([P1(1) P4(1)],[P1(2) P4(2)], 'color','black','LineWidth', 7);
        make_sp = line([P1(1) P6(1)],[P1(2) P6(2)], 'color','green','LineWidth', 10);
        hide = line([P1(1) P5(1)],[P1(2) P5(2)], 'color','white','LineWidth', 10);
        %fly = line([P1(1) P2(1)],[P1(2) P2(2)],'color','green','LineWidth', 5);
%         P1_circ = viscircles(P1,0.07,'color','black');
%         P2_circ = viscircles(P7,0.07,'color','black'); hold on
%         t = linspace(0,2*pi) ;
%         a = 0.5 ; b = 0.1 ; x = a*cos(t); y = b*sin(t) ;
%         R = Rot2D(pi/2-theta_fly(th)); XY_rotated = R*[x;y];
%         ellip = plot(XY_rotated(1,:),XY_rotated(2,:),'green');
        hold on
        [fly1_handles fly1_handles2]=plot_fly(handles.simulator,theta_fly(th),ye3(th),[0 1 0],1.5);
        if th==1
            % saveas(gca, sprintf('./gui_figures/simulator_gui_%s.png', patt), 'png');
        end
        axes(handles.trace);
        plot(t0(1:th),theta_fly(1:th)*180/pi,'green'); hold on; plot(t0(1:th),theta_obj(1:th)*180/pi,'black');box off;
        xlim(handles.trace,[t0(1),t0(end)]); ylim(handles.trace,[min([amp-amp/abs(amp)*180 amp]),max([amp-amp/abs(amp)*180 amp])]); set(gca,'XTickLabel',[]);
        title('Spot      {\color[rgb]{0, 1, 0}Fly}');
        ylabel('Angle (deg)');
        axes(handles.torque);
        % plot(t0(1:th),ye3(1:th),'green'); xlabel('Time (s)'); ylabel('Torque (Nm)');
        plot(t0(1:th),ye3(1:th),'red'); xlabel('Time (s)'); ylabel('Torque (Nm)');
        xlim(handles.torque,[t0(1),t0(end)]); ylim(handles.torque,[min([0 amp/abs(amp)*5e-11 ye3(1:th)]),max([0 amp/abs(amp)*5e-11 ye3(1:th)])]);box off;
        title('{\color[rgb]{1, 0, 0}Wings}');
        axes(handles.simulator);
        if get(handles.popupmenu2,'Value')==2
            frame = getframe(gcf); %get frame
            writeVideo(fly_simulator, frame);
        end
        pause(0.001);delete(obj);delete(fly1_handles);%delete(P1_circ);delete(ellip);%delete(fly);
        delete(fly1_handles2);
        delete(hide); delete(make_sp); %delete(P2_circ);
    elseif isequal(patt,'Grating')
        if any(check==th)
        P4 = 2.25*[cos(offset-theta_grat(th)) sin(offset-theta_grat(th))];
        P5 = 1.75*[cos(offset-theta_grat(th)) sin(offset-theta_grat(th))];
        P6 = 2.25*[cos(2*offset-theta_grat(th)) sin(2*offset-theta_grat(th))];
        P7 = 1.75*[cos(2*offset-theta_grat(th)) sin(2*offset-theta_grat(th))];
        P8 = 2.25*[cos(3*offset-theta_grat(th)) sin(3*offset-theta_grat(th))];
        P9 = 1.75*[cos(3*offset-theta_grat(th)) sin(3*offset-theta_grat(th))];
        P10 = 2.25*[cos(4*offset-theta_grat(th)) sin(4*offset-theta_grat(th))];
        P11 = 1.75*[cos(4*offset-theta_grat(th)) sin(4*offset-theta_grat(th))];
        
        P12 = 2.25*[cos(1.5*offset-theta_grat(th)) sin(1.5*offset-theta_grat(th))];
        P13 = 1.75*[cos(1.5*offset-theta_grat(th)) sin(1.5*offset-theta_grat(th))];
        P14 = 2.25*[cos(2.5*offset-theta_grat(th)) sin(2.5*offset-theta_grat(th))];
        P15 = 1.75*[cos(2.5*offset-theta_grat(th)) sin(2.5*offset-theta_grat(th))];
        P16 = 2.25*[cos(3.5*offset-theta_grat(th)) sin(3.5*offset-theta_grat(th))];
        P17 = 1.75*[cos(3.5*offset-theta_grat(th)) sin(3.5*offset-theta_grat(th))];
        P18 = 2.25*[cos(4.5*offset-theta_grat(th)) sin(4.5*offset-theta_grat(th))];
        P19 = 1.75*[cos(4.5*offset-theta_grat(th)) sin(4.5*offset-theta_grat(th))];
        
        obj = line([P1(1) P4(1)],[P1(2) P4(2)], 'color','black','LineWidth', 10);
        hide = line([P1(1) P5(1)],[P1(2) P5(2)], 'color','white','LineWidth', 10);
        obj1 = line([P1(1) P6(1)],[P1(2) P6(2)], 'color','black','LineWidth', 10);
        hide1 = line([P1(1) P7(1)],[P1(2) P7(2)], 'color','white','LineWidth', 10);
        obj2 = line([P1(1) P8(1)],[P1(2) P8(2)], 'color','black','LineWidth', 10);
        hide2 = line([P1(1) P9(1)],[P1(2) P9(2)], 'color','white','LineWidth', 10);
        obj3 = line([P1(1) P10(1)],[P1(2) P10(2)], 'color','black','LineWidth', 10);
        hide3 = line([P1(1) P11(1)],[P1(2) P11(2)], 'color','white','LineWidth', 10);
        
        obj4 = line([P1(1) P12(1)],[P1(2) P12(2)], 'color','black','LineWidth', 10);
        hide4 = line([P1(1) P13(1)],[P1(2) P13(2)], 'color','white','LineWidth', 10);
        obj5 = line([P1(1) P14(1)],[P1(2) P14(2)], 'color','black','LineWidth', 10);
        hide5 = line([P1(1) P15(1)],[P1(2) P15(2)], 'color','white','LineWidth', 10);
        obj6 = line([P1(1) P16(1)],[P1(2) P16(2)], 'color','black','LineWidth', 10);
        hide6 = line([P1(1) P17(1)],[P1(2) P17(2)], 'color','white','LineWidth', 10);
        obj7 = line([P1(1) P18(1)],[P1(2) P18(2)], 'color','black','LineWidth', 10);
        hide7 = line([P1(1) P19(1)],[P1(2) P19(2)], 'color','white','LineWidth', 10);        
        hold on
        [fly1_handles fly1_handles2]=plot_fly(handles.simulator,theta_fly(th),ye3(th),[1 0 0],1.5);
        if th==1
            % saveas(gca, sprintf('./gui_figures/simulator_gui_%s.png', patt), 'png');
        end
        axes(handles.trace);
        plot(t0(1:th),theta_fly(1:th)*180/pi,'red'); hold on; plot(t0(1:th),theta_grat(1:th)*180/pi,'black');box off;
        xlim(handles.trace,[t0(1),t0(end)]); ylim(handles.trace,[min([0 amp2]),max([0 amp2])]); set(gca,'XTickLabel',[]);
        title('Grating      {\color[rgb]{1, 0, 0}Fly}');
        ylabel('Angle (deg)');
        axes(handles.torque);
        plot(t0(1:th),ye3(1:th),'red');xlabel('Time (s)');box off; ylabel('Torque (Nm)');
        xlim(handles.torque,[t0(1),t0(end)]); ylim(handles.torque,[min([0 amp2/abs(amp2)*7e-11 ye3(1:th)]),max([0 amp2/abs(amp2)*7e-11 ye3(1:th)])])
        title('{\color[rgb]{1, 0, 0}Wings}');
        axes(handles.simulator);         
        if get(handles.popupmenu2,'Value')==2
            frame = getframe(gcf); %get frame
            writeVideo(fly_simulator, frame);
        end
        pause(0.001);
        delete(obj);%delete(P1_circ);%delete(fly);
        delete(hide); delete(fly1_handles);%delete(P2_circ);delete(ellip);
        delete(fly1_handles2);
        delete(obj1);delete(hide1);delete(obj2);delete(hide2);delete(obj3);delete(hide3);
        
        delete(obj4);delete(hide4);delete(obj5);delete(hide5);delete(obj6);delete(hide6);
        delete(obj7);delete(hide7); 
        end
    elseif isequal(patt,'Addition-based') || isequal(patt,'Graded efference copy') || isequal(patt,'All-or-none efference copy')
        if any(check==th)
        %P5_tr = viscircles(P1,2.5,'color','green','LineWidth', 20);
        %%%%%%%%%% Grating %%%%%%%%%%%%%%%%%%%%%%
        P4 = 2.25*[cos(offset-theta_grat(th)) sin(offset-theta_grat(th))];  % original value 2.75
        P5 = 1.75*[cos(offset-theta_grat(th)) sin(offset-theta_grat(th))];
        P6 = 2.25*[cos(2*offset-theta_grat(th)) sin(2*offset-theta_grat(th))];
        P7 = 1.75*[cos(2*offset-theta_grat(th)) sin(2*offset-theta_grat(th))];
        P8 = 2.25*[cos(3*offset-theta_grat(th)) sin(3*offset-theta_grat(th))];
        P9 = 1.75*[cos(3*offset-theta_grat(th)) sin(3*offset-theta_grat(th))];
        P10 = 2.25*[cos(4*offset-theta_grat(th)) sin(4*offset-theta_grat(th))];
        P11 = 1.75*[cos(4*offset-theta_grat(th)) sin(4*offset-theta_grat(th))];
        
        P12 = 2.25*[cos(1.5*offset-theta_grat(th)) sin(1.5*offset-theta_grat(th))];
        P13 = 1.75*[cos(1.5*offset-theta_grat(th)) sin(1.5*offset-theta_grat(th))];
        P14 = 2.25*[cos(2.5*offset-theta_grat(th)) sin(2.5*offset-theta_grat(th))];
        P15 = 1.75*[cos(2.5*offset-theta_grat(th)) sin(2.5*offset-theta_grat(th))];
        P16 = 2.25*[cos(3.5*offset-theta_grat(th)) sin(3.5*offset-theta_grat(th))];
        P17 = 1.75*[cos(3.5*offset-theta_grat(th)) sin(3.5*offset-theta_grat(th))];
        P18 = 2.25*[cos(4.5*offset-theta_grat(th)) sin(4.5*offset-theta_grat(th))];
        P19 = 1.75*[cos(4.5*offset-theta_grat(th)) sin(4.5*offset-theta_grat(th))];
        
        obj = line([P1(1) P4(1)],[P1(2) P4(2)], 'color',ones(1,3)*.6,'LineWidth', 10);
        hide = line([P1(1) P5(1)],[P1(2) P5(2)], 'color','white','LineWidth', 10);
        obj1 = line([P1(1) P6(1)],[P1(2) P6(2)], 'color',ones(1,3)*.6,'LineWidth', 10);
        hide1 = line([P1(1) P7(1)],[P1(2) P7(2)], 'color','white','LineWidth', 10);
        obj2 = line([P1(1) P8(1)],[P1(2) P8(2)], 'color',ones(1,3)*.6,'LineWidth', 10);
        hide2 = line([P1(1) P9(1)],[P1(2) P9(2)], 'color','white','LineWidth', 10);
        obj3 = line([P1(1) P10(1)],[P1(2) P10(2)], 'color',ones(1,3)*.6,'LineWidth', 10);
        hide3 = line([P1(1) P11(1)],[P1(2) P11(2)], 'color','white','LineWidth', 10);
        
        obj4 = line([P1(1) P12(1)],[P1(2) P12(2)], 'color',ones(1,3)*.6,'LineWidth', 10);
        hide4 = line([P1(1) P13(1)],[P1(2) P13(2)], 'color','white','LineWidth', 10);
        obj5 = line([P1(1) P14(1)],[P1(2) P14(2)], 'color',ones(1,3)*.6,'LineWidth', 10);
        hide5 = line([P1(1) P15(1)],[P1(2) P15(2)], 'color','white','LineWidth', 10);
        obj6 = line([P1(1) P16(1)],[P1(2) P16(2)], 'color',ones(1,3)*.6,'LineWidth', 10);
        hide6 = line([P1(1) P17(1)],[P1(2) P17(2)], 'color','white','LineWidth', 10);
        obj7 = line([P1(1) P18(1)],[P1(2) P18(2)], 'color',ones(1,3)*.6,'LineWidth', 10);
        hide7 = line([P1(1) P19(1)],[P1(2) P19(2)], 'color','white','LineWidth', 10);
        %P4_tr = viscircles(P1,2,'color','green','LineWidth', 20);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        P_bar = 2.25*[cos(offset-theta_obj(th)) sin(offset-theta_obj(th))];
        P_bar_h = 1.75*[cos(offset-theta_obj(th)) sin(offset-theta_obj(th))];
        obj_bar = line([P1(1) P_bar(1)],[P1(2) P_bar(2)], 'color','black','LineWidth', 10);
        hide_bar = line([P1(1) P_bar_h(1)],[P1(2) P_bar_h(2)], 'color','white','LineWidth', 10);
        hold on
        if isequal(patt,'Addition-based')
            [fly1_handles fly1_handles2]=plot_fly(handles.simulator,theta_fly(th),ye3(th),[1 0 1],1.5);
            if th==1
                % saveas(gca, sprintf('./gui_figures/simulator_gui_%s.png', patt), 'png');
            end
            axes(handles.trace);
            plot(t0(1:th),theta_fly(1:th)*180/pi,'magenta'); hold on; plot(t0(1:th),theta_grat(1:th)*180/pi,'color',ones(1,3)*.6);
            plot(t0(1:th),theta_obj(1:th)*180/pi,'black');        
            box off;
            xlim(handles.trace,[t0(1),t0(end)]); ylim(handles.trace,[min([0 -amp]),max([0 amp])]); set(gca,'XTickLabel',[]);
            title('Bar   {\color[rgb]{0.6, 0.6, 0.6}Grating}   {\color[rgb]{1, 0, 1}Fly}');
            ylabel('Angle (deg)');
            axes(handles.torque);
            plot(t0(1:th),ye3(1:th),'red'); hold on
            plot(t0(1:th),ye4(1:th),'--r');
            xlabel('Time (s)');box off; ylabel('Torque (Nm)');
            xlim(handles.torque,[t0(1),t0(end)]); ylim(handles.torque,[min([0 -amp/abs(amp)*10e-11 ye3(1:th)]),max([0 amp/abs(amp)*10e-11 ye3(1:th)])])
            title('{\color[rgb]{1, 0, 0}Wings}');
            axes(handles.simulator);
        elseif isequal(patt,'Graded efference copy')
            [fly1_handles fly1_handles2]=plot_fly(handles.simulator,theta_fly(th),ye3(th),[0.9290, 0.6940, 0.1250],1.5);
            if th==1
                % saveas(gca, sprintf('./gui_figures/simulator_gui_%s.png', patt), 'png');
            end
            axes(handles.trace);
            plot(t0(1:th),theta_fly(1:th)*180/pi,'color',[0.9290, 0.6940, 0.1250]); hold on; plot(t0(1:th),theta_grat(1:th)*180/pi,'color',ones(1,3)*.6);
            plot(t0(1:th),theta_obj(1:th)*180/pi,'black');        
            xlim(handles.trace,[t0(1),t0(end)]); ylim(handles.trace,[min([0 -amp]),max([0 amp])]); set(gca,'XTickLabel',[]);box off;
            title('Bar   {\color[rgb]{0.6, 0.6, 0.6}Grating}   {\color[rgb]{0.9290, 0.6940, 0.1250}Fly}');
            ylabel('Angle (deg)');
            axes(handles.torque);   
            plot(t0(1:th),ye3(1:th),'color','red'); hold on
            plot(t0(1:th),ye4(1:th),'--r');
            xlabel('Time (s)');box off; ylabel('Torque (Nm)');
            xlim(handles.torque,[t0(1),t0(end)]); ylim(handles.torque,[min([0 -amp/abs(amp)*7e-11 ye3(1:th)]),max([0 amp/abs(amp)*7e-11 ye3(1:th)])])
            title('{\color[rgb]{1, 0, 0}Wings}');
            axes(handles.simulator);
        elseif isequal(patt,'All-or-none efference copy')
            [fly1_handles fly1_handles2]=plot_fly(handles.simulator,theta_fly(th),ye3(th),[200 100 16]/256,1.5);
            if th==1
                % saveas(gca, sprintf('./gui_figures/simulator_gui_%s.png', patt), 'png');
            end
            axes(handles.trace);
            plot(t0(1:th),theta_fly(1:th)*180/pi,'color',[200 100 16]/256); hold on; plot(t0(1:th),theta_grat(1:th)*180/pi,'color',ones(1,3)*.6);
            plot(t0(1:th),theta_obj(1:th)*180/pi,'black');        
            xlim(handles.trace,[t0(1),t0(end)]); ylim(handles.trace,[min([0 -amp]),max([0 amp])]); set(gca,'XTickLabel',[]);box off;
            title('Bar   {\color[rgb]{0.6, 0.6, 0.6}Grating}   {\color[rgb]{0.7812, 0.3906, 0.0625}Fly}');
            ylabel('Angle (deg)');
            axes(handles.torque);   
            plot(t0(1:th),ye3(1:th),'color','red'); hold on
            plot(t0(1:th),ye4(1:th),'--r');
            xlabel('Time (s)');box off; ylabel('Torque (Nm)');
            xlim(handles.torque,[t0(1),t0(end)]); ylim(handles.torque,[min([0 -amp/abs(amp)*7e-11 ye3(1:th)]),max([0 amp/abs(amp)*7e-11 ye3(1:th)])])
            title('{\color[rgb]{1, 0, 0}Wings}');
            axes(handles.simulator);
        end
        if get(handles.popupmenu2,'Value')==2
            frame = getframe(gcf); %get frame
            writeVideo(fly_simulator, frame);
        end
        pause(0.001); delete(obj_bar); delete(fly1_handles);%delete(fly);
        delete(fly1_handles2);
        delete(hide_bar); %delete(P1_circ); %delete(P2_circ); delete(ellip);
        
        delete(obj);delete(hide);
        delete(obj1);delete(hide1);delete(obj2);delete(hide2);delete(obj3);delete(hide3);
        
        delete(obj4);delete(hide4);delete(obj5);delete(hide5);delete(obj6);delete(hide6);
        delete(obj7);delete(hide7);  
        end           
    elseif isequal(patt,'Graded vs All-or-none')
        if any(check==th)
        %P5_tr = viscircles(P1,2.5,'color','green','LineWidth', 20);
        %%%%%%%%%% Grating %%%%%%%%%%%%%%%%%%%%%%
        P4 = 2.25*[cos(offset-theta_grat(th)) sin(offset-theta_grat(th))];  % original value 2.75
        P5 = 1.75*[cos(offset-theta_grat(th)) sin(offset-theta_grat(th))];
        P6 = 2.25*[cos(2*offset-theta_grat(th)) sin(2*offset-theta_grat(th))];
        P7 = 1.75*[cos(2*offset-theta_grat(th)) sin(2*offset-theta_grat(th))];
        P8 = 2.25*[cos(3*offset-theta_grat(th)) sin(3*offset-theta_grat(th))];
        P9 = 1.75*[cos(3*offset-theta_grat(th)) sin(3*offset-theta_grat(th))];
        P10 = 2.25*[cos(4*offset-theta_grat(th)) sin(4*offset-theta_grat(th))];
        P11 = 1.75*[cos(4*offset-theta_grat(th)) sin(4*offset-theta_grat(th))];
        
        P12 = 2.25*[cos(1.5*offset-theta_grat(th)) sin(1.5*offset-theta_grat(th))];
        P13 = 1.75*[cos(1.5*offset-theta_grat(th)) sin(1.5*offset-theta_grat(th))];
        P14 = 2.25*[cos(2.5*offset-theta_grat(th)) sin(2.5*offset-theta_grat(th))];
        P15 = 1.75*[cos(2.5*offset-theta_grat(th)) sin(2.5*offset-theta_grat(th))];
        P16 = 2.25*[cos(3.5*offset-theta_grat(th)) sin(3.5*offset-theta_grat(th))];
        P17 = 1.75*[cos(3.5*offset-theta_grat(th)) sin(3.5*offset-theta_grat(th))];
        P18 = 2.25*[cos(4.5*offset-theta_grat(th)) sin(4.5*offset-theta_grat(th))];
        P19 = 1.75*[cos(4.5*offset-theta_grat(th)) sin(4.5*offset-theta_grat(th))];
        
        obj = line([P1(1) P4(1)],[P1(2) P4(2)], 'color',ones(1,3)*.6,'LineWidth', 10);
        hide = line([P1(1) P5(1)],[P1(2) P5(2)], 'color','white','LineWidth', 10);
        obj1 = line([P1(1) P6(1)],[P1(2) P6(2)], 'color',ones(1,3)*.6,'LineWidth', 10);
        hide1 = line([P1(1) P7(1)],[P1(2) P7(2)], 'color','white','LineWidth', 10);
        obj2 = line([P1(1) P8(1)],[P1(2) P8(2)], 'color',ones(1,3)*.6,'LineWidth', 10);
        hide2 = line([P1(1) P9(1)],[P1(2) P9(2)], 'color','white','LineWidth', 10);
        obj3 = line([P1(1) P10(1)],[P1(2) P10(2)], 'color',ones(1,3)*.6,'LineWidth', 10);
        hide3 = line([P1(1) P11(1)],[P1(2) P11(2)], 'color','white','LineWidth', 10);
        
        obj4 = line([P1(1) P12(1)],[P1(2) P12(2)], 'color',ones(1,3)*.6,'LineWidth', 10);
        hide4 = line([P1(1) P13(1)],[P1(2) P13(2)], 'color','white','LineWidth', 10);
        obj5 = line([P1(1) P14(1)],[P1(2) P14(2)], 'color',ones(1,3)*.6,'LineWidth', 10);
        hide5 = line([P1(1) P15(1)],[P1(2) P15(2)], 'color','white','LineWidth', 10);
        obj6 = line([P1(1) P16(1)],[P1(2) P16(2)], 'color',ones(1,3)*.6,'LineWidth', 10);
        hide6 = line([P1(1) P17(1)],[P1(2) P17(2)], 'color','white','LineWidth', 10);
        obj7 = line([P1(1) P18(1)],[P1(2) P18(2)], 'color',ones(1,3)*.6,'LineWidth', 10);
        hide7 = line([P1(1) P19(1)],[P1(2) P19(2)], 'color','white','LineWidth', 10);
        %P4_tr = viscircles(P1,2,'color','green','LineWidth', 20);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        P_bar = 2.25*[cos(offset-theta_obj(th)) sin(offset-theta_obj(th))];
        P_bar_h = 1.75*[cos(offset-theta_obj(th)) sin(offset-theta_obj(th))];
        obj_bar = line([P1(1) P_bar(1)],[P1(2) P_bar(2)], 'color','black','LineWidth', 10);
        hide_bar = line([P1(1) P_bar_h(1)],[P1(2) P_bar_h(2)], 'color','white','LineWidth', 10);
        P2_2 = 0.8*[cos(offset-theta_fly2(th)) sin(offset-theta_fly2(th))];
        hold on
        [fly1_handles fly1_handles2]=plot_fly(handles.simulator,theta_fly(th),ye3(th),[0.9290, 0.6940, 0.1250],1.5);
        [fly2_handles fly2_handles2]=plot_fly(handles.simulator,theta_fly2(th),yed3(th),[200 100 16]/256,1.4);
        if th==1
            % saveas(gca, sprintf('./gui_figures/simulator_gui_%s.png', patt), 'png');
        end
        axes(handles.trace);
        plot(t0(1:th),theta_fly(1:th)*180/pi,'color',[0.9290, 0.6940, 0.1250]); hold on; 
        plot(t0(1:th),theta_fly2(1:th)*180/pi,'color',[200 100 16]/256); plot(t0(1:th),theta_grat(1:th)*180/pi,'color',ones(1,3)*.6);
        plot(t0(1:th),theta_obj(1:th)*180/pi,'black');        
        box off;
        xlim(handles.trace,[t0(1),t0(end)]); ylim(handles.trace,[min([0 -amp]),max([0 amp])]); set(gca,'XTickLabel',[]);
        title('Bar   {\color[rgb]{0.6, 0.6, 0.6}Grating}   {\color[rgb]{0.9290, 0.6940, 0.1250}Fly (Graded)}   {\color[rgb]{0.7812, 0.3906, 0.0625}Fly (All-or-none)}');
        ylabel('Angle (deg)');
        axes(handles.torque);
        % plot(t0(1:th),ye3(1:th),'color',[0.9290, 0.6940, 0.1250]); hold on
        plot(t0(1:th),ye3(1:th),'color','red'); hold on
        plot(t0(1:th),ye4(1:th),'--r');
        % plot(t0(1:th),yed3(1:th),'yellow');xlabel('Time (s)');box off;ylabel('Torque (Nm)');
        plot(t0(1:th),yed3(1:th),'color',[200 100 16]/256);
        plot(t0(1:th),yed4(1:th),'color',[200 100 16]/256,'LineStyle','--');
        xlabel('Time (s)');box off;ylabel('Torque (Nm)');
        xlim(handles.torque,[t0(1),t0(end)]); ylim(handles.torque,[min([0 -amp/abs(amp)*10e-11 yed3(1:th)]),max([0 amp/abs(amp)*10e-11 yed3(1:th)])])
        title('{\color[rgb]{1, 0, 0}Wings (Graded)}      {\color[rgb]{0.7812, 0.3906, 0.0625}Wings (All-or-none)}');
        axes(handles.simulator);
        if get(handles.popupmenu2,'Value')==2
            frame = getframe(gcf); %get frame
            writeVideo(fly_simulator, frame);
        end
        pause(0.001); %delete(P1_circ);
        %delete(fly_eff); delete(fly_add);
        delete(fly1_handles); delete(fly2_handles);
        delete(fly1_handles2); delete(fly2_handles2);
        delete(obj_bar);
        delete(hide_bar); %delete(P2_circ); delete(ellip);
        
        delete(obj);delete(hide);
        delete(obj1);delete(hide1);delete(obj2);delete(hide2);delete(obj3);delete(hide3);
        
        delete(obj4);delete(hide4);delete(obj5);delete(hide5);delete(obj6);delete(hide6);
        delete(obj7);delete(hide7);       
        end
    end
end
set(handles.popupmenu1,'Enable', 'on')
set(handles.pushbutton1,'Enable', 'on')
set(handles.ampl,'Enable', 'on')
set(handles.ampl2,'Enable', 'on')
set(handles.onset,'Enable', 'on')
close(fly_simulator_v2)
 


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
restex = {'Response to bar' 'Response to spot' 'Response to grating' 'Response Addition model'...
                'Response Graded efference copy' 'Response All-or-none efference copy' 'Response Graded vs All-or-none'};
contents = get(handles.popupmenu1,'String');
ind = get(handles.popupmenu1,'Value');
set(handles.response,'String',restex{ind});
switch contents{ind}
    case 'Bar'
        patt = 'Bar';
%         data = load('./generated data/gui_variables_bar.mat');
%         theta_fly = data.theta_fly_bar; theta_obj = data.theta_obj_bar;
    case 'Spot'
        patt = 'Spot';
%         data = load('./generated data/gui_variables_spot.mat');
%         theta_fly = data.theta_fly_spot; theta_obj = data.theta_obj_spot;
    case 'Grating'
        patt = 'Grating';
%         data = load('./generated data/gui_variables_grat.mat');
%         theta_fly = data.theta_fly_grat; theta_obj = data.theta_obj_grat;
    case 'Addition-based'
        patt = 'Addition-based';
%         data = load('./generated data/gui_variables_add.mat');
%         theta_fly = data.theta_fly_add; theta_obj = data.theta_obj_add;
    case 'Graded efference copy'
        patt = 'Graded efference copy';
%         data = load('./generated data/gui_variables_eff.mat');
%         theta_fly = data.theta_fly_eff; theta_obj = data.theta_obj_eff;
    case 'All-or-none efference copy'
        patt = 'All-or-none efference copy';
%         data = load('./generated data/gui_variables_eff.mat');
%         theta_fly = data.theta_fly_eff; theta_obj = data.theta_obj_eff;
    case 'Graded vs All-or-none'
        patt = 'Graded vs All-or-none';
%         data = load('./generated data/gui_variables_add.mat');
%         theta_fly = data.theta_fly_add; theta_obj = data.theta_obj_add;
%         data2 = load('./generated data/gui_variables_eff.mat');
%         theta_fly2 = data2.theta_fly_eff; theta_obj2 = data2.theta_obj_eff;
end
set(handles.ampl,'Enable', 'on')


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function trace_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

function ampl_KeyPressFcn(hObject, eventdata, handles)
 % hObject    handle to pushbutton1 (see GCBO)
 % eventdata  reserved - to be defined in a future version of MATLAB
 % handles    structure with handles and user data (see GUIDATA)
 % get the value from the text box and convert to a number
 set(handles.ampl2,'Enable', 'on')
 %set(handles.pushbutton1,'Enable', 'on')
 
 
function f1_ = f1(t,y1)
    f1_= y1;
function f2_ = f2(t,y2,y3,y4)
    I = 6e-14;
    beta = 1e-11;
    f2_ = 1 / I * (-beta * y2 + y3 + y4);
function f3_ = f3(t,y2,y3, dot_theta_grat)
    kgrat = 5e-11;
    f3_ = 1 / 0.04 * (-y3 + kgrat * (dot_theta_grat - y2));
    
function R = Rot2D(angle)
    R = [cos(angle), -sin(angle); 
        sin(angle), cos(angle)];
     
% --- Executes during object creation, after setting all properties.
function time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

function [fly_handles fly_handles2] = plot_fly(haxes,angle,wb,color,line_width)
theta=-pi:0.01:pi;
angle = -180/pi*angle;
wb = -1e10*wb*180/pi;
bodyc = {[222 184 135]/255 [205 133  63]/255 [160  82  45]/255 [1 0 0] [1 0 0] [0 0 0] [0 0 0] [0 0 0] [0 0 0] [0 0 0] [0 0 0] [0 0 0]};
%calculate x, y for the head
x0=0.01*40*cosd(angle+90);y0=0.01*40*sind(angle+90);
rx=0.01*20;
ry=0.01*15;
x{1} = x0 + rx*cos(theta)*cosd(angle) - ry*sin(theta)*sind(angle);
y{1} = y0 + rx*cos(theta)*sind(angle) + ry*sin(theta)*cosd(angle);

%calculate x,y for the thorax
x0=0.01*7*cosd(angle+90);y0=0.01*7*sind(angle+90);
rx=0.01*20;
ry=0.01*20;
%calculate x and y points
x{2} = x0 + rx*cos(theta)*cosd(angle) - ry*sin(theta)*sind(angle);
y{2} = y0 + rx*cos(theta)*sind(angle) + ry*sin(theta)*cosd(angle);

%calculate x,y for the abdomen
x0=0.01*(-37*cosd(angle+90));y0=0.01*(-37*sind(angle+90));
rx=0.01*20;
ry=0.01*30;
%calculate x and y points
x{3} = x0 + rx*cos(theta)*cosd(angle) - ry*sin(theta)*sind(angle);
y{3} = y0 + rx*cos(theta)*sind(angle) + ry*sin(theta)*cosd(angle);

%calculate x,y for the left eye
x0=0.01*42*cosd(angle+90+15);y0=0.01*42*sind(angle+90+15);
rx=0.01*5;
ry=0.01*5;
%calculate x and y points
x{4} = x0 + rx*cos(theta)*cosd(angle) - ry*sin(theta)*sind(angle);
y{4} = y0 + rx*cos(theta)*sind(angle) + ry*sin(theta)*cosd(angle);

%calculate x,y for the right eye
x0=0.01*42*cosd(angle+90-15);y0=0.01*42*sind(angle+90-15);
rx=0.01*5;
ry=0.01*5;
%calculate x and y points
x{5} = x0 + rx*cos(theta)*cosd(angle) - ry*sin(theta)*sind(angle);
y{5} = y0 + rx*cos(theta)*sind(angle) + ry*sin(theta)*cosd(angle);

%calculate x,y for the left wing
x{6}=[0.01*25*cosd(angle+90+70+wb) 0.01*100*cosd(angle+90+70+wb)];y{6}=[0.01*25*sind(angle+90+70+wb) 0.01*100*sind(angle+90+70+wb)];

%calculate x,y for the right wing
x{7}=[0.01*25*cosd(angle+90-70+wb) 0.01*100*cosd(angle+90-70+wb)];y{7}=[0.01*25*sind(angle+90-70+wb) 0.01*100*sind(angle+90-70+wb)];

%calculate x,y for the left wing
x{8}=[0.01*25*cosd(-angle+90+70-wb) 0.01*100*cosd(angle-90-70+wb)];y{8}=[0.01*25*sind(angle-90-70+wb) 0.01*100*sind(angle-90-70+wb)];

%calculate x,y for the right wing
x{9}=[0.01*25*cosd(-angle+90-70-wb) 0.01*100*cosd(-angle+90-70-wb)];y{9}=[0.01*25*sind(-angle-90-70-wb) 0.01*100*sind(-angle-90-70-wb)];

%closing wing perimeter
x{10} = [0.01*100*cosd(angle+90+70+wb) 0.01*100*cosd(angle-90-70+wb)]; y{10} = [0.01*100*sind(angle+90+70+wb) 0.01*100*sind(angle-90-70+wb)];
x{11} = [0.01*100*cosd(angle+90-70+wb) 0.01*100*cosd(-angle+90-70-wb)]; y{11} = [0.01*100*sind(angle+90-70+wb) 0.01*100*sind(-angle-90-70-wb)];

% Reference line
x{12} = [0.01*100*cosd(angle+90) 0.01*100*cosd(-angle+90)];
y{12} = [0.01*100*sind(angle+90) 0.01*100*sind(-angle-90)];

% if(nargin<5)
    for j=1:length(x)
        if j==12
            fly_handles(j) = plot(x{j},y{j},'Color',color,'LineWIdth',line_width,'LineStyle','--');
            %fly_handles2(j) = patch(x{j},y{j},bodyc{j});
        elseif j>=6 && j<=11 && line_width ~= 1.4
            fly_handles(j) = plot(x{j},y{j},'Color','red','LineWIdth',line_width);
        elseif j>=6 && j<=11 && line_width == 1.4
            fly_handles(j) = plot(x{j},y{j},'Color',color,'LineWIdth',1.5);
        else
            fly_handles(j) = plot(x{j},y{j},'Color',color,'LineWIdth',line_width);
            fly_handles2(j) = patch(x{j},y{j},bodyc{j});
        end
    end


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ampl2_Callback(hObject, eventdata, handles)
% hObject    handle to ampl2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ampl2 as text
%        str2double(get(hObject,'String')) returns contents of ampl2 as a double


% --- Executes during object creation, after setting all properties.
function ampl2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ampl2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on ampl2 and none of its controls.
function ampl2_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to ampl2 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
set(handles.onset,'Enable', 'on')
%set(handles.pushbutton1,'Enable', 'on')



function onset_Callback(hObject, eventdata, handles)
% hObject    handle to onset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of onset as text
%        str2double(get(hObject,'String')) returns contents of onset as a double


% --- Executes during object creation, after setting all properties.
function onset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to onset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on onset and none of its controls.
function onset_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to onset (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
set(handles.pushbutton1,'Enable', 'on')
