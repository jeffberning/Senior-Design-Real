function [X,Q] = rungeKutta4traj(M,DT,X,Q,t0,tf,f,u)

for j = 1:M
    [k1,k1_q] = f(X,u,tf);
    [k2,k2_q] = f(X + DT/2*k1,u, tf);
    [k3,k3_q] = f(X + DT/2*k2,u,tf);
    [k4,k4_q] = f(X + DT*k3,u,tf);
    X = X + DT/6*(k1 + 2*k2 + 2*k3 + k4);
    Q = Q + DT/6*(k1_q + 2*k2_q + 2*k3_q + k4_q);
    t0 = t0 + DT;
end


end