function [yNew] = eulerStep(yCurrent,ydotCurrent,dt)
x = yCurrent(1); theta = yCurrent(2); xdot = yCurrent(3);
thetadot = yCurrent(4); xddot = ydotCurrent(3); 
thetaddot = ydotCurrent(4);

xnew = x + xdot*dt;
thetanew = theta + thetadot*dt;

xdotnew = xdot + xddot*dt;
thetadotnew = thetadot + thetaddot*dt;

yNew = [xnew thetanew xdotnew thetadotnew]';

end