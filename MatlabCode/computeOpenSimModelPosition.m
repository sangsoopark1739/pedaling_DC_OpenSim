

function pos = computeOpenSimModelPosition(states,t,osimModel,osimState, name1, name2)
    
    import org.opensim.modeling.*

    % Update state with current values  
    osimState.setTime(t);
    numVar = osimState.getNY();
    for i = 0:numVar-1
       osimState.updY().set(i, states(i+1));
    end


    % Update the state velocity calculations
    osimModel.computeStateVariableDerivatives(osimState);

    simbodyEngine = osimModel.getSimbodyEngine();

    globalposition = Vec3(0.00,  0.00, 0.00);
    localposition = Vec3(0.00,  0.00, 0.00);

    simbodyEngine.transformPosition(osimState, osimModel.getBodySet().get(name1), localposition, osimModel.getBodySet().get(name2), globalposition);
    
    px = globalposition.get(0);
    py = globalposition.get(1);
    pz = globalposition.get(2);
    
    pos = [px;py;pz];
