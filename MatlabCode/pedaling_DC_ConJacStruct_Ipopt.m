function JacobianStructure = pedaling_DC_ConJacStruct_Ipopt
%
% Returns the sparsity structure of the constraint Jacobian for IPOPT; it
% can either be assumed to be dense, or the actual sparsity structure can be
% loaded from a MAT file and stored as a global variable
%

%-----------------------------------------------------------------------

% Uncomment this section to use a sparse Jacobian
global JacStruct

JacobianStructure = JacStruct;
JacobianStructure = sparse(JacobianStructure);

end


