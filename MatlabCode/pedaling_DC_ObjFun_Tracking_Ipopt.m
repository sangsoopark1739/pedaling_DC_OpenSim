function f = pedaling_DC_ObjFun_Tracking_Ipopt(X)
%
% Computes the objective function value for the 2D Walking model simulations
%	 X = current set of optimization parameters
%   f = objective function value
%
% By setting specific flags below the obj func value may contain terms
% associated with tracking kinematics and/or GRFs, minimizing muscle activations
% or other quantities the user specifies
%
% Authors: Brian Umberger & Leng-Feng Lee, UMass Amherst
% 

    global auxdata
    global trackinginfo
   
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
    
    %
    Xnum     = trackinginfo.Xnum;
    Nodenum  = trackinginfo.Nodenum;
    
    Exp_meanForce = trackinginfo.meanForce;
    Exp_stdForce = trackinginfo.stdForce;    
    Exp_meanAng = trackinginfo.meanAngle;
    Exp_stdAng = trackinginfo.stdAngle;
    
    Mod_Force = trackinginfo.Mod_Force;
    Mod_Angle = trackinginfo.Mod_Angle;
    
    
    if Xnum <= N*2*Ncoord
      
      temp_Force = computeOpenSimModelPedalForce(states(Nodenum,:)',dc_time(Nodenum),osimModel,osimState)';
      Mod_Force(Nodenum,1:2) = temp_Force(1:2);
      Mod_Force(Nodenum,3:4) = temp_Force(4:5);
      
      temp_Angle = computeOpenSimModelPedalAngle_both(states(Nodenum,:)',dc_time(Nodenum),osimModel,osimState)'; 
      Mod_Angle(Nodenum,1:2) = temp_Angle;
      
      Nodenum = Nodenum + 1;   
      if Nodenum > N
          Nodenum = 1;
      end
      trackinginfo.Nodenum = Nodenum;
    end
    
    %%
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
    
	% Angles that you would like to track
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
   
    f = ActWeight*SumIntegAct + TrackWeight*totalTrackingErr;
end









