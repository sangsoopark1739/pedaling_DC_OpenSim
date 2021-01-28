function f = pedaling_DC_ObjFun_Ipopt(X)
%
% Computes the objective function value to solve the optimal control pedaling problem
%	 X = current set of optimization parameters
%   f = objective function value

% Dependencies: 
% computeOpenSimModelPedalForce.m
% computeOpenSimModelPedalAngle_both.m

% By setting specific flags below the obj func value may contain terms
% associated with tracking kinematics and/or GRFs, minimizing muscle activations
% or other quantities the user specifies
%
% Authors: Brian Umberger & Leng-Feng Lee, UMass Amherst
%          Updated by Sangsoo Park & Brian Umberger, UMass Amherst

    global auxdata
    global trackinginfo

    pause(0.000001)
    
    % Auxiliary data that are passed in
    osimModel = auxdata.model;
    N         = auxdata.N;
    Nmuscles  = auxdata.Nmuscles;
    Ncoord    = auxdata.Ncoord;
    Nstates   = auxdata.Nstates;
    Ncontrols = auxdata.Ncontrols; 
    tFinal    = auxdata.tFinal;

    % Import the OpenSim modeling classes
    import org.opensim.modeling.*

    % Check to see if model state is initialized by checking size
    if(osimModel.getWorkingState().getNY() == 0)
       osimState = osimModel.initSystem();
    else
       osimState = osimModel.updWorkingState(); 
    end


    % Create two variables "states" and "controls" from X
    % X = [jnt states (11 pos, 11 vel), mus states (18 act, 18 fib), 18 controls] at N time points
    states = zeros(N,Nstates); %pre-allocate size
    for i = 1:Nstates
        states(:,i) = X(N*(i-1)+1:N*i,1); %column: state; row: time ptn
    end

    controls = zeros(N,Ncontrols); %pre-allocate size
    for i = 1:Ncontrols
        controls(:,i) = X(Nstates*N + N*(i-1)+1:Nstates*N + N*i,1); %column: controls; row: time ptn
    end
    
    % Create the time column from current value of tFinal
    tFinal  = X(end,1); 
    dc_time = tFinal*linspace(0,1,N);
    h       = X(end,1)/(N-1);   % time interval between nodes

    
    % Cal model pedal forces, Compute pedal forces from Spring 
    Mod_Force = zeros(N,4); %XYZ
    Mod_Angle = zeros(N,2); % Left/Right pedal angles from the model
    for i = 1:N
         temp_Force = computeOpenSimModelPedalForce(states(i,:)',dc_time(i),osimModel,osimState)'; % compute pedal forces
         Mod_Force(i,:) = [temp_Force(1:2);temp_Force(4:5)];
         Mod_Angle(i,1:2) = computeOpenSimModelPedalAngle_both(states(i,:)',dc_time(i),osimModel,osimState)'; % compute pedal angles
    end
    
    %
    trackinginfo.Mod_Force = Mod_Force;% left-right
    trackinginfo.Mod_Angle = Mod_Angle; % left-right
    Exp_meanForce = trackinginfo.meanForce;
    Exp_stdForce = trackinginfo.stdForce;
    Exp_meanAng = trackinginfo.meanAngle;
    Exp_stdAng = trackinginfo.stdAngle;
    
    %
    force_err_sum = 0;
    angle_err_sum = 0;
    
    % Summed, squared Ang tracking error
    numForceTrackVar = 4;
    numAngleTrackVar = 2;
	
	% Forces that you would like to track, normalized by std -> unitless (van den Bogert - 2011)
    for i = 1:numForceTrackVar
      for j = 1:N
         force_err_sum = force_err_sum + ((Exp_meanForce(j,i) - Mod_Force(j,i))/Exp_stdForce(j,i))^2;
      end
    end
    
	% Tracking the pedal angles, normalized by std -> unitless (van den Bogert - 2011)
    for i = 1:numAngleTrackVar % both pedal angles
        for j=1:N          
          angle_err_sum = angle_err_sum + ((Exp_meanAng(j,i) - Mod_Angle(j,i))/Exp_stdAng(j,i))^2;
        end
    end
    
	% Muscle activations
	Mus_Act = zeros(N,Nmuscles); %pre-allocate size
    j = 0;
    for i = (2*Ncoord+1):2:(Nstates) % muscle activation and fiber length are arranged in alternate manner
        j = j + 1;
        Mus_Act(:,j) = X(N*(i-1)+1:N*i,1); %column: muscle; row: nodes (time steps)
    end
	
	SumIntegAct = 0;
    for i = 1:size(Mus_Act,2)
       SumIntegAct = SumIntegAct + (h * trapz(Mus_Act(:,i).^3));
    end

	% 
	totalAct = (SumIntegAct/Nmuscles)/N;
	totalTrackingErr = ((force_err_sum + angle_err_sum)/(numForceTrackVar+numAngleTrackVar))/N;
    	
    % Weightings
    ActWeight = 1;
    TrackWeight = 1;
	
	% Objective function
    f = ActWeight*SumIntegAct + TrackWeight*totalTrackingErr;
end









