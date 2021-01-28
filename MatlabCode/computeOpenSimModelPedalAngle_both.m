
function PedalAng = computeOpenSimModelPedalAngle_both(states,t,osimModel,osimState)

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
    
    toes_L_pos = zeros(1,3);
    calcn_L_pos = zeros(1,3);
    toes_R_pos = zeros(1,3);
    calcn_R_pos = zeros(1,3);
    
    toes_L_pos = computeOpenSimModelPosition(states,t,osimModel,osimState,'toes_l', 'ground')';
    calcn_L_pos = computeOpenSimModelPosition(states,t,osimModel,osimState, 'calcn_l', 'ground')';
    
    toes_R_pos = computeOpenSimModelPosition(states,t,osimModel,osimState,'toes_r', 'ground')';
    calcn_R_pos = computeOpenSimModelPosition(states,t,osimModel,osimState, 'calcn_r', 'ground')';   
    
    Vec_calcn_toes_L = toes_L_pos - calcn_L_pos; % in GCS
    Vec_calcn_toes_R = toes_R_pos - calcn_R_pos; % in GCS
    Vec_X_unit = [1,0,0];
        
    Y_L = Vec_calcn_toes_L;
    Y_R = Vec_calcn_toes_R;
    X = Vec_X_unit;
        
    mag_Vec_Y_L = sqrt(Y_L(1)^2 + Y_L(2)^2);
    mag_Vec_Y_R = sqrt(Y_R(1)^2 + Y_R(2)^2);
    mag_Vec_X = sqrt(X(1)^2 + X(2)^2);
    
    dir_L = cross(Y_L,X); % counterclockwise (+) 
    dir_R = cross(Y_R,X);
    
    tempAng_L = acos(dot(X(1:2),Y_L(1:2))/(mag_Vec_X*mag_Vec_Y_L))*180/pi;
    tempAng_R = acos(dot(X(1:2),Y_R(1:2))/(mag_Vec_X*mag_Vec_Y_R))*180/pi;
    
    if dir_L(3) > 0
        pedal_L = tempAng_L;
    else
        pedal_L = -tempAng_L;
    end
    
    if dir_R(3) > 0
        pedal_R = tempAng_R;
    else
        pedal_R = -tempAng_R;
    end
    
    L_PedalAng = pedal_L * pi / 180; % degrees to radians
    R_PedalAng = pedal_R * pi / 180; % degrees to radians
    
    PedalAng = [L_PedalAng,R_PedalAng]; % Should check which one for which col
    
end