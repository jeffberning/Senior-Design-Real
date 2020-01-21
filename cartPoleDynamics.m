% All Final Units in SI
% This code will use a Euler integration approach to solve the full
% dynamics of the system
close all
clear all
clc
%% Solve for Dynamics (Already solved using D'Alambert/Lagrange)
syms mo ml l theta g F thetadot xdot real
H = [mo+ml ml*l*cos(theta); ml*l*cos(theta) ml*l^2];
Hinv = inv(H);
Hinv = simplify(Hinv);
C = [0 -ml*l*thetadot*sin(theta); 0 0];
Tau = [0;ml*g*l*sin(theta)];
U = [F;0];
qddot = simplify(Hinv*(U-C*[xdot;thetadot] - Tau));

%% Simulation
params = systemParams();
t0 = 0;
t_end = 10; % 10 secs, can be changed
dt = 0.01; % 0.01 secs, also can be changed


% let state be y = [x theta xdot thetadot]';

x0 = 0; % initial horizontal position
xdot0 = 0; % initial horizontal velocity
theta0 = 90*pi/180; % initial angle
thetadot0 = 0; % initial angular velocity
F0 = 0; % initial force
% Calculate initial second derivatives
[xddot0,thetaddot0] = eomCalcs(x0,theta0,xdot0,thetadot0,F0);


%%% Initialize all the vectors (create empty vectors of the correct length
%%% for everything except time
tVec = linspace(t0,t_end,(t_end-t0)/dt+1); % Time Vector
xVec = zeros(1,length(tVec));
xdotVec = zeros(1,length(tVec));
thetaVec = zeros(1,length(tVec));
thetadotVec = zeros(1,length(tVec));
xddotVec = zeros(1,length(tVec));
thetaddotVec = zeros(1,length(tVec));
Fvec = zeros(1,length(tVec));
% arbitrary control policy of F = -10*xdot (damping)


xVec(1) = x0; xdotVec(1) = xdot0; xddotVec(1) = xddot0;
thetaVec(1) = theta0; thetadotVec(1) = thetadot0; 
thetaddotVec(1) = thetaddot0;

for iCnt = 2:length(tVec)
    % Euler integration to get new values
    [xnew,thetanew,xdotnew,thetadotnew] = eulerStep(xVec(iCnt-1),...
        thetaVec(iCnt-1),xdotVec(iCnt-1),thetadotVec(iCnt-1),...
        xddotVec(iCnt-1),thetaddotVec(iCnt-1),dt);
    Fvec(iCnt) = -30*xdotnew;
    [xddotnew,thetaddotnew] =  eomCalcs(xnew,thetanew,xdotnew,...
        thetadotnew,Fvec(iCnt));
    
    xVec(iCnt) = xnew; xdotVec(iCnt) = xdotnew; xddotVec(iCnt) = xddotnew;
    thetaVec(iCnt) = thetanew; thetadotVec(iCnt) = thetadotnew;
    thetaddotVec(iCnt) = thetaddotnew;
    
end




%% Plotting (Animation)

figure(1)
hold on
axis([-2 2 -2 2])
l = params.l;
yCart = 0;
xCart0 = xVec(1);
cart = plot(xCart0,yCart,'kx');

xLoad0 = xCart0 + l*sin(thetaVec(1));
yLoad0 = yCart - l*cos(thetaVec(1));
Load = plot(xLoad0,yLoad0,'rx');
for jCnt = 1:length(tVec) 
    xCart = xVec(jCnt);
    xLoad = xCart + l*sin(thetaVec(jCnt));
    yLoad = yCart - l*cos(thetaVec(jCnt));
    set(cart,'XData',xCart);
    set(cart,'YData',yCart);
    set(Load,'XData',xLoad);
    set(Load,'YData',yLoad);
    drawnow
    pause(0.00001) 
end







