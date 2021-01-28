function J = pedaling_DC_ConJac_Sparse_Ipopt(X)
%
% Function to compute the Jacobian of the constaints for IPOPT using
% finite forward differences; step size is selected uniquely for each
% variable; a grouping strategy based on the sparsity pattern is used
% to calculate the full Jacobian with as few evaluations of the full
% constraint vector as possible.
%
% This version is called when the sparsity pattern of the Jacobian is known.
%
%   X = current set of optimization parameters
%   J = Jacobian of constraints with respect to X


global FacJ JacStruct JacGrps

t          = 0; % dummy variable; not used but needed by numjac
fX         = pedaling_DC_ConFun_Ipopt(X);
thresh     = 1e0*ones(size(X));
vectorized = 0;
func       = @pedaling_DC_ConFun_NumJac;

S = JacStruct;

[dfdx,FacJ,JacGrps] = numjac(func,t,X,fX,thresh,FacJ,vectorized,S,JacGrps);

J = dfdx;

    
end

