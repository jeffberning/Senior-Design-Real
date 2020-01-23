function [yNew] = RK4Step(yCurrent,ydotCurrent,dt,F)
x = yCurrent(1); theta = yCurrent(2); xdot = yCurrent(3);
thetadot = yCurrent(4); xddot = ydotCurrent(3); 
thetaddot = ydotCurrent(4);

[xddot1 thetaddot1] = eomCalcs(yCurrent,F);
k1 = [xdot, thetadot, xddot1, thetaddot1]';
[xddot2 thetaddot2] = eomCalcs(yCurrent+dt/2*k1,F);
k2 = [xdot+dt/2*k1(3), thetadot+dt/2*k1(4),xddot2,thetaddot2]';
[xddot3 thetaddot3] = eomCalcs(yCurrent+dt/2*k2,F);
k3 = [xdot+dt/2*k2(3), thetadot+dt/2*k2(4), xddot3,thetaddot3]';
[xddot4 thetaddot4] = eomCalcs(yCurrent+dt*k3,F);
k4 = [xdot+dt*k3(3) thetadot+dt*k3(4), xddot4, thetaddot4]';


yNew = yCurrent + dt/6*(k1 + 2*k2 + 2*k3 + k4);
% xnew = x + xdot*dt;
% thetanew = theta + thetadot*dt;
% 
% xdotnew = xdot + xddot*dt;
% thetadotnew = thetadot + thetaddot*dt;
% 
% yNew = [xnew thetanew xdotnew thetadotnew]';

end