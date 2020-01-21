function [xnew,thetanew,xdotnew,thetadotnew] = eulerStep(x,theta,xdot,...
    thetadot,xddot,thetaddot,dt)


xnew = x + xdot*dt;
thetanew = theta + thetadot*dt;

xdotnew = xdot + xddot*dt;
thetadotnew = thetadot + thetaddot*dt;

end