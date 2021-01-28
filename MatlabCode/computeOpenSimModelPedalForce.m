function PedalForce = computeOpenSimModelPedalForce(states,t,osimModel,osimState)
% This function sets the model to a particular state, returns the three
% components of the GRFs applied to both feet.
%
% Author: Brian Umberger, UMass Amherst
%

% Import the OpenSim modeling classes
import org.opensim.modeling.*

% Update state with current values  
osimState.setTime(t);
numVar = osimState.getNY();
for i = 0:numVar-1
   osimState.updY().set(i, states(i+1,1));
end

% Update the state velocity calculations
osimModel.computeStateVariableDerivatives(osimState);

% Right GRFs

rPedalForce_Array = ArrayDouble();
rPedalForce_Array = osimModel.getForceSet().get('Spring_r').getRecordValues(osimState);
rAPPF = 1 * rPedalForce_Array.getitem(0);
rVPF  = 1 * rPedalForce_Array.getitem(1);
rMLPF = -1 * rPedalForce_Array.getitem(2);

% Left GRFs
lPedalForce_Array = ArrayDouble();
lPedalForce_Array = osimModel.getForceSet().get('Spring_l').getRecordValues(osimState);
lAPPF = 1 * lPedalForce_Array.getitem(0);
lVPF  = 1 * lPedalForce_Array.getitem(1);
lMLPF = -1 * lPedalForce_Array.getitem(2);

% Put GRFs in a variable to return to the calling fucntion
PedalForce = [lAPPF lVPF lMLPF rAPPF rVPF rMLPF];


end
