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


%% Load Control
uVec = load('optimalControlFirst.mat');
uVec = uVec.u_opt;
%% Simulation
params = systemParams();
t0 = 0;
t_end = 50; % 10 secs, can be changed
%dt = t_end/30;
dt = 0.01; % 0.01 secs, also can be changed


% let state be y = [x theta xdot thetadot]';
in2met = 0.0254;
%x0 = -12*in2met; % initial horizontal position
x0 = 0;
xdot0 = 0; % initial horizontal velocity
theta0 = pi/2; % initial angle
thetadot0 = 0; % initial angular velocity
y0 = [x0 theta0 xdot0 thetadot0]';
F0 = 0; % initial force
% Calculate initial second derivatives
[xddot0,thetaddot0] = eomCalcs(y0,F0);


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

%%% LINEAR DYNAMICS MATRICES
l = params.l; mo = params.mo; ml = params.ml; g = params.g; 
A = [0 0 1 0; 0 0 0 1 ; 0 g*ml/mo 0 0; 0 -(g*ml + g*mo)/(l*mo) 0 0];
B = [0 0 1/mo -1/(l*mo)]';
Q = [1000 0 0 0; 0 10 0 0; 0 0 1000 0; 0 0 0 10];
R = .001;
K = lqr(A,B,Q,R,0);
for iCnt = 2:length(tVec)
    % Euler integration to get new values
    yCurrent = [xVec(iCnt-1) thetaVec(iCnt-1) xdotVec(iCnt-1) thetadotVec(iCnt-1)]';
    ydotCurrent = [xdotVec(iCnt-1) thetadotVec(iCnt-1) xddotVec(iCnt-1) thetaddotVec(iCnt-1)]';
    %[yNew] = eulerStep(yCurrent,ydotCurrent,dt);
    %Fvec(iCnt) = -30*yCurrent(3);
   %Fvec(iCnt) = uVec(iCnt);
    %Fvec(iCnt) = -K*yCurrent;
    [yNew] = RK4Step(yCurrent,ydotCurrent,dt,Fvec(iCnt));
    
    [xddotnew,thetaddotnew] =  eomCalcs(yNew,Fvec(iCnt));
    
    xVec(iCnt) = yNew(1); xdotVec(iCnt) = yNew(3);
    xddotVec(iCnt) = xddotnew;
    thetaVec(iCnt) = yNew(2); thetadotVec(iCnt) = yNew(4);
    thetaddotVec(iCnt) = thetaddotnew;
    
end




%% Plotting (Animation)
animate = 1;
if animate == 1
figure(1)
hold on
axis([-2 2 -2 2])
l = params.l;
yCart = 0;
xCart0 = xVec(1);
cart = plot(xCart0,yCart,'kx');
boxBottom = plot(xCart0 + [-.25 .25],yCart + [-.15 -.15],'k-','LineWidth',2);
boxTop = plot(xCart0 + [-.25 .25],yCart + [.15 .15],'k-','LineWidth',2);
boxLeft = plot(xCart0 + [-.25 -.25],yCart + [-.15 .15],'k-','LineWidth',2);
boxRight = plot(xCart0 + [.25 .25],yCart + [-.15 .15],'k-','LineWidth',2);
xLoad0 = xCart0 + l*sin(thetaVec(1));
yLoad0 = yCart - l*cos(thetaVec(1));
Load = plot(xLoad0,yLoad0,'r.','MarkerSize',30);
pole = plot([xCart0 xLoad0],[yCart yLoad0],'k-','LineWidth',2);
textTime = text(-2, 2.25, 'Time: 0.0 (s)','fontsize',16);
for jCnt = 1:length(tVec) 
    xCart = xVec(jCnt);
    xLoad = xCart + l*sin(thetaVec(jCnt));
    yLoad = yCart - l*cos(thetaVec(jCnt));
    set(cart,'XData',xCart);
    set(cart,'YData',yCart);
    set(Load,'XData',xLoad);
    set(Load,'YData',yLoad);
    set(pole,'XData',[xCart xLoad])
    set(pole,'YData',[yCart,yLoad])
    set(boxBottom,'XData',xCart + [-.25 .25])
    set(boxBottom,'YData',yCart + [-.15 -.15])
    set(boxTop,'XData',xCart + [-.25 .25])
    set(boxTop,'YData',yCart + [.15 .15])
    set(boxLeft,'XData',xCart + [-.25 -.25])
    set(boxLeft,'YData',yCart + [-.15 .15])
    set(boxRight,'XData',xCart + [.25 .25])
    set(boxRight,'YData',yCart + [-.15 .15])
    set(textTime,'string',['Time: ',num2str(tVec(jCnt)),' (s)'])
    drawnow
    %pause(0.0001) 
end

end




figure(2)
hold on
plot(tVec,xVec/in2met)
xlabel('Time (s)')
ylabel('Position (in)')
figure(3)
plot(tVec,thetaVec*180/pi)

xlabel('Time (s)')
ylabel('Angle (Degrees)')

figure(4)

plot(tVec,Fvec)
xlabel('Time (s)')
ylabel('Force (N)')








