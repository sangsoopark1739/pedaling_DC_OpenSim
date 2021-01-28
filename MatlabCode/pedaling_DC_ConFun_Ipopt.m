function ceq = pedaling_DC_ConFun_Ipopt(X)
    %
    % Computes the nonlinear equality and inequality contraints
    %	X      = current set of optimization parameters
    %  ceq    = vector of equality constraints
    % dependencies: computeOpenSimModelXdot_parallel.m
    % Authors: Brian Umberger & Leng-Feng Lee, UMass Amherst
    %          Updated by Sangsoo Park & Brian Umberger, UMass Amherst

    global auxdata

    % Import the OpenSim modeling classes
    import org.opensim.modeling.*

    % extract the nesessary auxiliary data
    N         = auxdata.N;
    Nstates   = auxdata.Nstates;
    Nmuscles  = auxdata.Nmuscles;
    Ncontrols = auxdata.Ncontrols;
    osimModel = auxdata.model;
    Ncoord    = auxdata.Ncoord;

    % Check to see if model state is initialized by checking size
    if(osimModel.getWorkingState().getNY() == 0)
       osimState = osimModel.initSystem();
    else
       osimState = osimModel.updWorkingState(); 
    end

    % First, we need to extract all of the states and controls by node (time step)
    % from X. This will result in two matrices, states and controls, that have a
    % number of rows equal to the number of nodes and a number of columns
    % equal to the number of states/controls, as appropriate.

    states=zeros(N,Nstates); %pre-allocate size
    for i = 1:Nstates
        states(:,i) = X(N*(i-1)+1:N*i,1); %column: state; row: nodes (time steps)
    end

    controls = zeros(N,Ncontrols); %pre-allocate size
    for i = 1:Ncontrols
        controls(:,i) = X(Nstates*N + N*(i-1)+1:Nstates*N + N*i,1); %column: controls; row: nodes (time steps)
    end
    
    % Create the time column from current value of tFinal
    tFinal  = X(end,1);
    dc_time = tFinal*linspace(0,1,N);
    h       = X(end,1)/(N-1);   % time interval between nodes
    
    % ------------ Scaling factors for the dyanmics & continuity constraints -----------

    DynConstScaling(1:7)   = 1.0;  % velocities,
    DynConstScaling(8:14)  = 0.01; % accelerations, 
    DynConstScaling(15:50) = 0.1;  % muscle dynamics, 9 muscles * 2

    DynConstScaling = DynConstScaling';
    
    % ----------- Compute the equality constraints, ceq --------------------------------
    % Compute the constraint violation using backward Euler method
    states_dot = zeros(N-1,Nstates);
    for i = 1:Nstates
       for j = 1:N-1
          states_dot(j,i) = (states(j+1,i)-states(j,i))/h;
       end
    end

    % Get state derivatives from OpenSim
    x_dot = zeros(N-1,Nstates);
    for i = 1:N-1
       x_dot(i,:) = computeOpenSimModelXdot_parallel(states(i+1,:)',controls(i+1,:)',dc_time(i+1))';
    end

    % Evaluate the constraint violoation at every time step except the last node, N
    ceq_temp= zeros(N-1,Nstates);
    for i = 1:Nstates
        ceq_temp(:,i) = states_dot(1:end,i) - x_dot(1:end,i);
        ceq_temp(:,i) = ceq_temp(:,i) * DynConstScaling(i); % scale the dynamics constraints
    end

    % Re-arrange the constraint violations in one long vector
    for i = 1:Nstates
        ceq((N-1)*(i-1)+1:(N-1)*i,:) = ceq_temp(:,i);
    end
    

    % Additional task constraints
    % All states except crank angle are same at beginning and end
    for i = 1:6
       ceq(end+1,1)= states(end,i) - states(1,i); % joint angles
    end
    
    ceq(end+1,1)= states(end,7) - states(1,7) - 2*pi; % one crank revolution at 30 rpm;
    
    for i = 8:Nstates
       ceq(end+1,1)= states(end,i) - states(1,i); % the rest of the states
    end
    
    % controls at end must be same as controls in beginning
    for i = 1:Ncontrols
       ceq(end+1,1)= controls(end,i) - controls(1,i);
    end
    
    
    ceq(end+1,1) = states(1,7);  % inital crank angle = 0

 end
 