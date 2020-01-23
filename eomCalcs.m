function [xddot thetaddot] = eomCalcs(y,F)
%%% This function just takes the current state and calculates the two
%%% second derivatives (x double dot and theta double dot)
params = systemParams();
l = params.l; mo = params.mo; ml = params.ml; g = params.g; 
x = y(1); theta = y(2); xdot = y(3); thetadot = y(4);

xddot = (l*ml*sin(theta)*thetadot^2 + F + g*ml*cos(theta)*sin(theta))...
    /(ml + mo - ml*cos(theta)^2);

thetaddot =  -(l*ml*cos(theta)*sin(theta)*thetadot^2 + F*cos(theta) + ...
    g*ml*sin(theta) + g*mo*sin(theta))/(l*(ml + mo - ml*cos(theta)^2));

end