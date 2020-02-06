function [ydot] = linearDynamics(y,F)
params = systemParams();
l = params.l; mo = params.mo; ml = params.ml; g = params.g; 
A = [0 0 1 0; 0 0 0 1 ; 0 g*ml/mo 0 0; 0 -(g*ml + g*mo)/(l*mo) 0 0];

B = [0 0 1/mo -1/(l*mo)]';

ydot = A*y + B*F;

end



