function f = pedaling_DC_ObjFun_NumJac(t,X)
%
% This is just a wrapper function that allows the use of Matlab's
% built-in numjac.m function for numerical derivatives;
% numjac expects t and X to be passed in while IPOPT expects only
% X (and possibly auxdata) to be passed in.

f = pedaling_DC_ObjFun_Ipopt(X);

end
