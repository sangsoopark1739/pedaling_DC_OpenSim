%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main file: run IPOPT to solve optimal control pedaling problems using DC
% Dependencis:
% pedaling_DC_ConFun_Ipopt.m
% pedaling_DC_ObjFun_Ipopt.m
% pedaling_DC_ObjFunGrad_Tracking.m
% pedaling_DC_ConFun_Ipopt.m
% pedaling_DC_ConJac_Sparse_Ipopt.m
% pedaling_DC_ConJacStruct_Ipopt.m
% pedaling_DC_Callback.m
% WriteOpenSimStatesFile.m
% WriteOpenSimControlFile.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
close all

global auxdata
global FacG FacJ JacGrps JacStruct
global states_all controls_all
global trackinginfo
global osimState osimModel

% The number of nodes (Density of the optimal solution)
N = 31;
% Offset angle between the pedal angle in the model and experimental pedaling data in degree
offset = 13; 

% data in 360 data points/cycle, 0-359 degrees
% Mean left pedal force (AP and V), pedal angle, EMG data of fifteen recreational cyclists (Park and Caldwell, In Press)
load([pwd,'\Experimental_Data\TwoL_Tracking_Data.mat']);

% Time for completion of one crank revolution 
tFinal = 1.914;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% Prepare the experimental data
% Organzie the experimental data 
% 1) Interpolate into the number of nodes (Left leg)
old_crankAng = (0:1:360)';
new_crankAng = (0:360/(N-1):360)';
pedalAng = interp1(old_crankAng, Two_L.Ang_Two_L(:,1), new_crankAng);
pedalAng_std = interp1(old_crankAng, Two_L.Ang_Two_L(:,2), new_crankAng);
L_pedalAng = (pedalAng+offset).*pi/180; % degrees to radians
L_pedalAng_std = pedalAng_std.*pi/180; % degrees to radians

L_pedalForceAP = interp1(old_crankAng, Two_L.APF_Two_L(:,1), new_crankAng);
L_pedalForceAP_std = interp1(old_crankAng, Two_L.APF_Two_L(:,2), new_crankAng);
L_pedalForceV = interp1(old_crankAng, Two_L.VF_Two_L(:,1), new_crankAng);
L_pedalForceV_std = interp1(old_crankAng, Two_L.VF_Two_L(:,2), new_crankAng);

% 2) Produce pedal forces and angle patterns for the right leg by shifting 180 degrees of those patterns in the left leg
R_pedalAng = [L_pedalAng(round(N/2):N,1);L_pedalAng(2:round(N/2),1)];
R_pedalAng_std = [L_pedalAng_std(round(N/2):N,1);L_pedalAng_std(2:round(N/2),1)];
R_pedalForceAP = [L_pedalForceAP(round(N/2):N,1);L_pedalForceAP(2:round(N/2),1)];
R_pedalForceAP_std = [L_pedalForceAP_std(round(N/2):N,1);L_pedalForceAP_std(2:round(N/2),1)];
R_pedalForceV = [L_pedalForceV(round(N/2):N,1);L_pedalForceV(2:round(N/2),1)];
R_pedalForceV_std = [L_pedalForceV_std(round(N/2):N,1);L_pedalForceV_std(2:round(N/2),1)];

% 3) 'trankinginfo' contains the experimental pedal force and angle patterns
trackinginfo.meanForce = [L_pedalForceAP, L_pedalForceV, R_pedalForceAP, R_pedalForceV];
trackinginfo.stdForce = [L_pedalForceAP_std,L_pedalForceV_std,R_pedalForceAP_std,R_pedalForceV_std];
trackinginfo.meanAngle = [L_pedalAng, R_pedalAng];
trackinginfo.stdAngle = [L_pedalAng_std, R_pedalAng_std];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize FacG, FacJ and JacGrps (used by numjac)?
FacG    = [];
FacJ    = [];
JacGrps = [];

% Import the OpenSim modeling classes
import org.opensim.modeling.*

% Read in the osim model
nameofModel = [pwd,'\Model\Pedaling_Model.osim'];
osimModel = Model(nameofModel); 

% Initialize the model (this builds the system and initialize the state)
osimState = osimModel.initSystem();

% Get the number of states, coordinates, muscles and controls from the model;
% in this case the number of controls equals the number of muscles
Nstates       = osimModel.getNumStateVariables();
Ncontrols     = osimModel.getNumControls();
Ncoord        = osimModel.getNumCoordinates(); 
model_muscles = osimModel.getMuscles();
Nmuscles      = model_muscles.getSize();

% Auxiliary data to be passed to the optimizer
auxdata.model      = osimModel;
auxdata.N          = N;
auxdata.Nstates    = Nstates;
auxdata.Ncontrols  = Ncontrols; 
auxdata.Nmuscles   = Nmuscles;
auxdata.Ncoord     = Ncoord;

% Get the names of the states from the model
states_all = cell(Nstates,1);
for i = 1:Nstates
   states_all(i,1) = cell(osimModel.getStateVariableNames().getitem(i-1));
end

controls_all = cell(Ncontrols,1);
for i = 1:Ncontrols
   currentMuscle = model_muscles.get(i-1);
   controls_all(i,1) = cell(currentMuscle.getName());
end

h = tFinal/(N-1);      % time interval between nodes
dc_time = tFinal*linspace(0,1,N)';

auxdata.h          = h;
auxdata.time       = dc_time;
auxdata.tFinal     = tFinal;

%%%%%%%%%%%%%%%%%% IC states
[file_input, pathname] = uigetfile({'*.sto', 'OpenSim States Files (*.sto)'}, ...
                         'Select the initial states file','MultiSelect', 'off');
temp_s = importdata(strcat(pathname,file_input)); % import states data
old_time_s = temp_s.data(:,1);         % time 
old_data_s = temp_s.data(:,2:end);     
old_text_s = temp_s.textdata(7,2:end); 

% Arrange the initial guess by nodes and states
x0_temp = zeros(N,Nstates);  % pre-allocate space
chk_counter =0;
for j = 1:size(states_all,1) 
    for k = 1:size(old_data_s,2)
        if strcmp(old_text_s(k),states_all(j)) == 1
            % interpolate initial guess to the defined temporal grid
            x0_temp(:,j) = interp1(old_time_s,old_data_s(:,k),dc_time);
            chk_counter = chk_counter+1;
        end
    end
end
if (chk_counter ~= size(states_all,1))
    disp('Inconsistent number of states for inital guess!');
    disp('Check the input file...');
    return
end

% Arrange the initial guess into a column vector
% [ [pos(t_0) ... pos(t_N)], [vel(t_0) ... vel(t_N)], ... etc]
for i = 1:Nstates
    X0(N*(i-1)+1:N*i,:) = x0_temp(:,i);
end

% Get the names of the controls/muscles from the model (same in this case)

% Load the file that contain the initial guess for the controls (excitations)
[file_input, pathname] = uigetfile({'*.sto', 'OpenSim Controls (excitation) Files (*.sto)'}, ...
                                   'Select the initial controls file','MultiSelect', 'off');
temp_i = importdata(strcat(pathname,file_input)); % import controls data (1st column is time)
old_time_i = temp_i.data(:,1);         % time 
old_data_i = temp_i.data(:,2:end);     % the controls, 
old_text_i = temp_i.textdata(7,2:end); % controls names, 

% Arrange the initial guess by nodes and controls
u0_temp = zeros(N,Ncontrols);
for j = 1:size(controls_all,1)
    for k = 1:size(old_data_i,2)
        if strcmp(old_text_i(k),controls_all(j)) == 1
            % interpolate to the temporal grid
            u0_temp(:,j) = interp1(old_time_i,old_data_i(:,k),dc_time);
        end
    end
end

% Append the initial guess for the controls to the end of X0
for i = 1:Ncontrols
    X0(Nstates*N + N*(i-1)+1 : Nstates*N + N*i,:) = u0_temp(:,i);
end

% Check: make sure both files (states and controls) are consistent
if (old_time_i(end) ~= old_time_s(end))
    disp('Time stamp of states and controls for initial guess do not match!');
    disp('Check the input files...');
    return
end

% Append the initial guess of tFinal to the end of X0
X0(end+1,1) = tFinal;

% Load the MAT file that contains the sparsity structure of the constraint
% Jaoobian and the groups of columns that do not share any non-zero rows;
[file_input, pathname] = uigetfile( ...
{'*.mat', 'Jacobian MAT File (*.mat)'},'Select MAT file','MultiSelect', 'off');
Jacobian  = load(strcat(pathname,file_input));
JacStruct = Jacobian.Jstruct;
figure('Name','Jacobian Matrix','NumberTitle','off')
spy(JacStruct)

% Number of unknowns and constraints
NX             = size(X0,1);
C              = pedaling_DC_ConFun_Ipopt(X0);
Nconst         = size(C,1);
auxdata.NX     = NX;
auxdata.Nconst = Nconst;

% Bounds on the optimization parameters
Pos_LB1(1:N) =  -2.0944;            Pos_UB1(1:N) = 2.0944;      % Hip R
Pos_LB2(1:N) =  0;                  Pos_UB2(1:N) = 2.0944;      % Knee R
Pos_LB3(1:N) = -pi/2;               Pos_UB3(1:N) = pi/2;        % Ankle R

Pos_LB4(1:N) =  -2.0944;            Pos_UB4(1:N) = 2.0944;      % Hip L
Pos_LB5(1:N) =  0;                  Pos_UB5(1:N) = 2.0944;      % Knee L
Pos_LB6(1:N) = -pi/2;               Pos_UB6(1:N) = pi/2;        % Ankle L

Pos_LB7(1:N) = 0;                   Pos_UB7(1:N) = 6*pi;        % Crank

Pos_LB = [Pos_LB1,Pos_LB2,Pos_LB3,Pos_LB4,Pos_LB5,Pos_LB6,Pos_LB7];
Pos_UB = [Pos_UB1,Pos_UB2,Pos_UB3,Pos_UB4,Pos_UB5,Pos_UB6,Pos_UB7];

Vel_LB1(1:6*N) =  -10*pi;           Vel_UB1(1:6*N) = 10*pi;   % Hip Knee, Ank ang vel
Vel_LB2(1:N) = 0;                   Vel_UB2(1:N) = 10*pi;

Act_LB(1:Nmuscles*N) = 0.011;        Act_UB(1:Nmuscles*N) = 0.99;
Fib_LB(1:Nmuscles*N) = 0.011;        Fib_UB(1:Nmuscles*N) = 0.99;       

for i = 1:Nmuscles % Activation and fiber length are arranged in alternate manner
    Act_Fib_LB(2*N*(i-1)+1:2*N*(i)) = [Act_LB(N*(i-1)+1:N*(i)) Fib_LB(N*(i-1)+1:N*(i))];
    Act_Fib_UB(2*N*(i-1)+1:2*N*(i)) = [Act_UB(N*(i-1)+1:N*(i)) Fib_UB(N*(i-1)+1:N*(i))];
end

%---
Con_LB(1:Nmuscles*N) = 0.011;        Con_UB(1:Nmuscles*N) = 0.99; %constraints on controls

% Set LB and UB both to fFinal to have a fixed final time (IPOPT handles this
% very nicely)
Time_LB = tFinal;    Time_UB = tFinal; % pedal cycle time

lb = [Pos_LB Vel_LB1 Vel_LB2 Act_Fib_LB Con_LB Time_LB]';
ub = [Pos_UB Vel_UB1 Vel_UB2 Act_Fib_UB Con_UB Time_UB]';

% Functions and options for IPOPT
funcs.objective         = @pedaling_DC_ObjFun_Ipopt;
funcs.gradient          = @pedaling_DC_ObjFunGrad_Tracking;
funcs.constraints       = @pedaling_DC_ConFun_Ipopt;
funcs.jacobian          = @pedaling_DC_ConJac_Sparse_Ipopt;
funcs.jacobianstructure = @pedaling_DC_ConJacStruct_Ipopt;
funcs.iterfunc          = @pedaling_DC_Callback;

options.cl = zeros(Nconst,1);
options.cu = zeros(Nconst,1);
options.lb = lb;
options.ub = ub;

options.ipopt.hessian_approximation = 'limited-memory';
options.ipopt.recalc_y_feas_tol = 1e-02;

options.ipopt.max_cpu_time    = 43200;
options.ipopt.max_iter        = 1000;
options.ipopt.tol             = 1e02;
options.ipopt.dual_inf_tol    = 1e02;
options.ipopt.constr_viol_tol = 1e-03;
options.ipopt.compl_inf_tol   = 1e-03;
options.ipopt.acceptable_tol  = 1e-03;
options.ipopt.acceptable_iter = 10;

options.ipopt.bound_frac = 0.001;  %0.0001;
options.ipopt.bound_push = 0.001;  %0.0001;
options.ipopt.ma57_pre_alloc        = 2.0;

options.ipopt.print_level = 5;

% IPOPT
[Xopt, info] = ipopt(X0,funcs,options);

% Save Lagrange multiplier info for possible re-start (warm start)
save info.mat info Xopt

X_state_opt = zeros(N,Nstates); %pre-allocate size
for i = 1:Nstates
    X_state_opt(:,i) = Xopt(N*(i-1)+1:N*i,1);
end

X_controls_opt = zeros(N,Ncontrols); %pre-allocate size
for i = 1:Ncontrols
    X_controls_opt(:,i) = Xopt(Nstates*N + N*(i-1)+1:Nstates*N + N*i,1);
end

%X_time_opt = zeros(N,1); %pre-allocate size
X_time_opt = Xopt(end,1)*linspace(0,1,N)';

% Create data structure for the states file
StatesData = struct();
StatesData.name = [char(osimModel.getName()), '_Optimal_States_', char(date)];
StatesData.nRows = size(X_state_opt, 1);
StatesData.nColumns = Nstates+1; %All the states + time
StatesData.inDegrees = false;
StatesData.labels = cell(1,StatesData.nColumns); 
StatesData.labels{1}= 'time';
for j = 2:1:StatesData.nColumns
   StatesData.labels{j} = char(states_all(j-1));
end
StatesData.data = [X_time_opt, X_state_opt];
WriteOpenSimStatesFile(StatesData)

% Create data structure for the controls file
ControlData = struct();
ControlData.name = [char(osimModel.getName()), '_Optimal_Controls_', char(date)];
ControlData.nRows = size(X_controls_opt, 1);
ControlData.nColumns = Ncontrols+1; %All the controls + time
ControlData.inDegrees = false;
ControlData.labels = cell(1,ControlData.nColumns); 
ControlData.labels{1}= 'time';
for j = 2:1:ControlData.nColumns
   ControlData.labels{j} = char(controls_all(j-1));
end
ControlData.data = [X_time_opt, X_controls_opt];
WriteOpenSimControlFile(ControlData)

save optimalRes.mat StatesData ControlData
