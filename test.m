clc;
close all;
clear;
format long g;
addpath ../tstFcn;
%% 

costFunc = @rosen;              
dim      = 50;
varMin   =-30;
varMax   = 30;

% costFunc = @Sphere;           
% dim     = 50;
% varMin   =-100;
% varMax   = 100;

res = pso_alg(costFunc, dim, varMin, varMax, 200, 200);
plot(res);
res(end)