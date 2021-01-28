function b = pedaling_DC_Callback(t, f, currentData)
%
% This callback function can optionally be called to write out
% temp states and control files at the end of every iteration
% to monitor the ongoing progress of an optimization.
%
% This current version can cause the optimization to crash if the current
% iterate is really bad. This can probably be fixed by doing some size checking
% and not writing out the temp files for that iteration.
%
% Authors: Leng-Feng Lee & Brian Umberger, UMass Amherst
%

global auxdata
global states_all controls_all

% extra parameters from auxdata
%dc_time   = auxdata.time;
N         = auxdata.N;
Nstates   = auxdata.Nstates;
Ncontrols = auxdata.Ncontrols;

x = currentData.x;

% %-------------write temp State file----------------------
% for i = 1:Nstates,
%     x_state(:,i) = x(N*(i-1)+1:N*i,1); %column: state; row: time ptn
% end

% Create the time column from current value of tFinal
dc_time = x(end,1)*linspace(0,1,N)';

%-------------write temp state file----------------------
for i = 1:Nstates
    x_state(:,i) = x(N*(i-1)+1:N*i,1); %column: state; row: time ptn
end
OutputData = struct();
OutputData.name = ['pedaling_DC_Temp_States'];
OutputData.nRows = size(x_state, 1);
OutputData.nColumns = Nstates+1; %All the states + time
OutputData.inDegrees = false;
OutputData.labels = cell(1,OutputData.nColumns);
OutputData.labels{1}= 'time';
for j = 2:1:OutputData.nColumns
   OutputData.labels{j} = char(states_all(j-1));
end
OutputData.data = [dc_time, x_state]; %time and corresponding optimized states
WriteOpenSimStatesFile(OutputData);

%-------------write temp control file----------------------
for i = 1:Ncontrols
    x_controls(:,i) = x(Nstates*N + N*(i-1)+1:Nstates*N + N*i,1); %column: controls; row: time ptn
end
ControlData = struct();
ControlData.name = [ 'pedaling_DC_Temp_Control'];
ControlData.nRows = size(x_controls, 1);
ControlData.nColumns = Ncontrols+1; %All the inputs + time
ControlData.inDegrees = false;
ControlData.labels = cell(1,ControlData.nColumns);
ControlData.labels{1}= 'time';
for j = 2:1:ControlData.nColumns
   ControlData.labels{j} = char(controls_all(j-1));
end
ControlData.data = [dc_time, x_controls]; %time and corresponding optimized states
WriteOpenSimControlFile(ControlData);
%---------------------------------------------

b = true;

end
