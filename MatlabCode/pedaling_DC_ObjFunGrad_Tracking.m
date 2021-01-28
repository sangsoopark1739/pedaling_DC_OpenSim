function G = pedaling_DC_ObjFunGrad_Tracking(X)
%
% Function to compute the gradient of the objective function for IPOPT
% using finite forward differences; 
%   X = current set of optimization parameters
%   G = gradient of ObjFun with respect to X
%

global trackinginfo

n   = length(X);
G   = zeros(n,1);
fx  = pedaling_DC_ObjFun_Ipopt(X);
del = sqrt(eps);
xperturb = X;

trackinginfo.Nodenum = 1;

for i = 1:n

   trackinginfo.Xnum = i;

   xperturb(i) = xperturb(i) + del;
   fperturb = pedaling_DC_ObjFun_Tracking_Ipopt(xperturb);
   G(i) = (fperturb-fx) / del;
   xperturb(i) = X(i);
   
end


end





