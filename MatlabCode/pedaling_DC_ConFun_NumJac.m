 function ceq = pedaling_DC_ConFun_NumJac(t,X)
%
% This is just a wrapper function that allows the use of Matlab
% built-in numjac.m function for numerical derivatives;
% numjac expects t and X to be passed in while IPOPT expects only
% X (and possibly auxdata) to be passed in.


ceq = pedaling_DC_ConFun_Ipopt(X);


end
