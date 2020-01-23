function params = systemParams()
%%% This function just creates a variable for params that contains all the
%%% relevant parameters. Having this function should reduce the clutter in
%%% the main codes

in2met = 0.0254;
lb2kg = 0.453592;
wireLength = 30.25; %inches
loadWeight = 5; % lbf (assume on earth)
cartWeight = 2; %lbf (assume on earth) just initial guess
params.l = wireLength*in2met;
params.ml = loadWeight*lb2kg;
params.g = 9.81;
params.mo = cartWeight*lb2kg;



end 