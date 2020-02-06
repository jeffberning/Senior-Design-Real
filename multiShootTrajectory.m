% multiShootTrajectory

clear
clc
close all
p = genpath('C:\Matlab_programs'); % This creates a path to the folder where CasADi lives
addpath(p) % This adds the previous path to the working directory
import casadi.* % loads CasADi


N = 30; % Number of control intervals

params = systemParams();
l = params.l; mo = params.mo; ml = params.ml; g = params.g; 

xDes = [12*.0254 0 0 0]'; % desired [x,theta,xdot,thetadot]

% Declare model Variables
t = MX.sym('t',1); % Timing variable (this one goes from zero to 1, not actual time)
x = MX.sym('x',4); % STATE [x theta xdot thetadot]
x0 = MX.sym('x0',4); % Initial State
tf = MX.sym('tf',1); % Time Scaling Variable
u = MX.sym('u',1); % Control Force
time = t*tf; % this is actual time

[xdot] = trajDynamics(x,u); % get dynamics [xdot thetadot xddot thetaddot]'
%[xdot] = linearDynamics(x,u);
%L = sqrt((x(1)-xDes(1))^2)%*tf;% + (x(2)-xDes(2))^2 + (x(3)-xDes(3))^2 + (x(4)-xDes(4))^2);%*tf;
L = tf;
f = Function('f',{x,u,tf},{xdot,L}); % creates a function that takes in the state, control, tf, and gives dynamics and cost

%%% INTEGRATOR (RK but different formulation)

X0 = MX.sym('X0',4);
T0 = MX.sym('t0');
M = 5; % # of RK steps per shooting segment
DT = 1/(N*M); % time step but for RK 
X = X0;
Q = 0; % Initial Cost
t0 = T0;
[X,Q] = rungeKutta4traj(M,DT,X,Q,t0,tf,f,u); % RK gives back new x and new cost
% function that takes in x, u, tf and gives back new x and q
F = Function('F',{T0,X0,u,tf},{X,Q},{'t0','x0','u','tf'},{'xf','qf'});

% SET UP STUFF
problem.vars = {}; % create empty structure array to initialize
problem.vars_init = []; % array of initial guesses
problem.vars_lb = []; % lower bounds
problem.vars_ub = []; % upper bounds
problem.cost = 0; % initial cost
problem.constraints = {}; % initialize constraints
problem.constraints_lb = []; % bounds on constraints
problem.constraints_ub = []; % bounds on constraints

u_inds = []; % indices
x_inds = [];
other_inds = [];
p_inds = [];
tf_inds = [];
t_inds = [];
F_inds = [];

% bounds and guesses for the time
tf_guess = 10;
tf_lb = 1;
tf_ub = 25;

[problem,tf,tf_inds(end+1,:)] = add_var(problem,'tf',1,tf_guess,tf_lb,tf_ub);

% initialize the state
xinit = [0 0 0 0]';
xinit_lb = xinit;
xinit_ub = xinit;

% intialize the control, but this one is just a guess
uinit = 0; % just a guess
uinit_lb = -10;
uinit_ub = 10;
cost = 0; % initialize the cost
[problem, Xk, x_inds(end+1,:)] = add_var(problem,'X0',4,xinit,xinit_lb,xinit_ub);
[problem,Uk,u_inds(end+1,:)] = add_var(problem,'U0',1,uinit,uinit_lb,uinit_ub);


[problem,Xk,Xk_end,x_inds,Uk,u_inds,cost] = multiShoot(problem,Xk,x_inds,Uk,u_inds,tf,F,N,cost);



problem.cost = cost;

prob = struct('f',problem.cost,'x',vertcat(problem.vars{:}),'g',vertcat(problem.constraints{:}));
s_opts = struct;
s_opts.ipopt.max_iter = 1000;
solver = nlpsol('solver','ipopt',prob,s_opts);

tic
sol = solver('x0',problem.vars_init,'lbx',problem.vars_lb,'ubx',problem.vars_ub,...
    'lbg',problem.constraints_lb,'ubg',problem.constraints_ub);
toc
w_opt = full(sol.x);

tf_opt = reshape(w_opt(tf_inds(:)),size(tf_inds));
x_opt = reshape(w_opt(x_inds(:)),N+1,length(x));
u_opt = reshape(w_opt(u_inds(:)),size(u_inds));

figure(1)
plot(1:N+1,x_opt(:,1)./0.0254)
xlabel('Step')
ylabel('Position (in)')

figure(2)
plot(1:N+1,u_opt)
xlabel('Step')
ylabel('Force (N)')

figure(2)
function [problem, var, var_inds] = add_var(problem, name, len, init, lb, ub)
    import casadi.*
    var = MX.sym(name, len);
    problem.vars = {problem.vars{:}, var};
    problem.vars_lb = [problem.vars_lb;lb];
    problem.vars_ub = [problem.vars_ub;ub];
    problem.vars_init = [problem.vars_init; init];
    
    ind_end = length(problem.vars_lb);
    ind_begin = ind_end-len+1;
    var_inds = ind_begin: ind_end;
end

function [problem] = add_constraint(problem, constr, lb, ub)
    problem.constraints = {problem.constraints{:},constr};
    problem.constraints_lb = [problem.constraints_lb; lb];
    problem.constraints_ub = [problem.constraints_ub; ub];
end

