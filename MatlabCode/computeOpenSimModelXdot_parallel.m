function x_dot = computeOpenSimModelXdot_parallel(states_local,controls_local,t)
% This function sets the model to a particular state, applies the controls
% (excitations) and returns the derivaties of the state variables.
%
% Author: Brian Umberger, Vinh Nguyen, UMass Amherst
%
% Note: This function draws from the OpenSimPlantFunction.m function by Daniel
% Jacobs (Stanford), which is part of the OpenSim dynamic walking example


global osimModel osimState;

% Import the OpenSim modeling classes
import org.opensim.modeling.*

% Get some variables  
numVar = osimState.getNY();
nControls = osimModel.getNumControls();

% Pass in time
osimState.setTime(t);

% Set state values using a vector
y = Vector(numVar, 0.0);
for i = 0:numVar-1
   y.set(i, states_local(i+1,1));
end
osimState.setY(y);

% Create model control vector
modelControls = Vector(nControls, 0.0);

% Pass control values into a control vector
for i = 0:nControls-1
    modelControls.set(i, controls_local(i+1,1));
end

% Update the state velocity calculations
osimModel.computeStateVariableDerivatives(osimState);

% Pass control vector into state
osimModel.setControls(osimState,modelControls);

% Update the derivative calculations in the state variable
osimModel.computeStateVariableDerivatives(osimState);

% Set output variable to the new state
 x_dot = zeros(numVar,1);
 VecDot = osimState.getYDot(); % using a vector is faster
 
for i = 0:1:numVar-1
    x_dot(i+1,1) = VecDot.get(i);
end

end
