% ------------------------------------------------------------------------------
% Author: Narayan P. Subramaniyam
% Affiliation: MET Faculty, Tampere University
% Email: narayan.subramaniyam@tuni.fi
%
% Description:
% This MATLAB code is developed as part of my research. Feel free to reuse 
% or modify this code, provided that you give proper attribution by citing 
% the associated paper:
% 
% [Insert Full Citation of Your Paper Here]
%
% License:
% This code is licensed under a Creative Commons Attribution 4.0 International License.
% You are free to share and adapt the material for any purpose, even commercially,
% under the following terms:
% 1. You must give appropriate credit, provide a link to the license, and indicate 
%    if changes were made.
% 2. You must cite the original paper if you use this code in your work.
%
% For more details on the license, visit:
% https://creativecommons.org/licenses/by/4.0/
% ------------------------------------------------------------------------------

% The code returns first-order, second-order and total-order sobol indices

clear all
addpath(strcat(pwd, '/func'))
f=@nmm_jr_SA;
dim=13;
n=20000;
lb = [0  0  0  1  25  6.5  0  0  0  0  2 0.5 0.3];
ub = [10 50 200 1000 140 110 1350 1080 340 340 9 7.5 0.8];
res=sobol_nmm(f,dim,n, 'sampler', 'lhs', 'lb', lb,'ub', ub);