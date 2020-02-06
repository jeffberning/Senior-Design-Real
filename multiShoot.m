function [problem,Xk,Xk_end,x_inds,Uk,u_inds,cost] = multiShoot(problem,Xk,x_inds,Uk,u_inds,tf,F,N,cost)

for k = 0:1:N-1
    Fk = F('t0',k*1/(N),'x0',Xk,'tf',tf,'u',Uk);
    Xk_end = Fk.xf;
    cost = cost + Fk.qf;
    [problem,Xk,x_inds(end+1,:)] = add_var(problem,['X_' num2str(k+1)],4,[0 0 0 0]',[-inf -inf -inf -inf]',[inf inf inf inf]');
    [problem,Uk,u_inds(end+1,:)] = add_var(problem,['U_' num2str(k+1)],1,0,-25,25);
    [problem] = add_constraint(problem,Xk(2), -0.4362, 0.4363);
    [problem] = add_constraint(problem,Xk_end-Xk, [0 0 0 0]',[0 0 0 0]'); 
    
    if k == (N-1)
         [problem] = add_constraint(problem, Xk_end , [(12-.1)*.0254, .0087, 0, 0]', [(12+.1)*.0254,.0087 , 0, 0]');
     end
end


end