clear all
clc
close all

syms x xdot theta thetadot F real
syms l ml g mo real

y = [x theta xdot thetadot]';
yEq = [0 0 0 0]';
FEq = 0;

ydot(1) = y(3); ydot(2) = y(4);

ydot(3) = (l*ml*sin(theta)*thetadot^2 + F + g*ml*cos(theta)*sin(theta))...
    /(ml + mo - ml*cos(theta)^2);

ydot(4) = -(l*ml*cos(theta)*sin(theta)*thetadot^2 + F*cos(theta) + ...
    g*ml*sin(theta) + g*mo*sin(theta))/(l*(ml + mo - ml*cos(theta)^2));
ydot = ydot';

A = simplify(jacobian(ydot,y))';


B = simplify(jacobian(ydot,F));

AVal = subs(A,[x theta xdot thetadot]',yEq);
AVal = subs(AVal,F,FEq);

BVal = subs(B,[x theta xdot thetadot]',yEq);
BVal = subs(BVal,F,FEq);






