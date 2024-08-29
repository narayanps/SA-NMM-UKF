function run_jr_sim(j)
% ------------------------------------------------------------------------------
% Author: Narayan P. Subramaniyam
% Affiliation: MET Faculty, Tampere University
% Email: narayan.subramaniyam@tuni.fi
%
% Description:
% This MATLAB code is developed as part of my research. Feel free to reuse 
% or modify this code, provided that you give proper attribution by citing 
% the associated paper. This code estimates parameters of NMM model of any
% required combination as inpt by variable j. Please check C and list for
% the combination and values. List contains the parameter combinations 1-2,
% 1-5, 1-6 pertaining to A, B, mu, sigma, a, b, C1, C2, C3, C4, v0, e0, r. Check the paper 
% 
%
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

load list
load C
combs = list(j,:);
vals = C(j,:);
addpath(strcat(pwd,'/func'))
SNR= [1 3 5 10];
T=100000;
T_ignore=10000;
n_s=6;
n_p=length(combs);
p = 90 + 30.*randn(1,T+T_ignore);
H=zeros(1,n_s+n_p);
H(2)=1; H(3)=-1; 
state_0=zeros(n_s,1);
state=zeros(T,n_s);
y=zeros(1,T);
f=@nmm_jr_param_est;
f_nmm=@nmm_jr_ukf_gen;
h=@h_meas;
dt=0.001;
X = [6; 70; 90; 30; 100; 50 ; 135; 108 ; 33.75; 33.75; 6; 2.5; 0.56];
X_vals = repmat(X, [1 T+T_ignore]);
X_vals(list(j,1),50001:end) = C(j,1);
X_vals(list(j,2),50001:end) = C(j,2);   
state(1,:) = f(X_vals(:,1), state_0, p(1));
y(1) = h(H(1:n_s), state(1,:)');
for i=2:T+T_ignore
    state(i,:) = f(X_vals(:,i), state(i-1,:)',p(i));
    y(i) = h(H(1:n_s), state(i,:)');
end
y=y-mean(y);
sig_var = var(y);
R = sig_var./SNR(4);
noise = chol(R)*randn(1,T+T_ignore);
y=y+noise;
y=y(T_ignore+1:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Q0=1e-10.*eye(n_s+n_p);
Q0(7,7)=1e-2;
R0=1;
x_hat0=zeros(n_s+n_p,1);
for k=1:length(vals)
x_hat0(n_s+k,1) = vals(1,k) - vals(1,k)*0.1;
Q0(n_s+k,n_s+k)=1e-2;
end
P_xx0=100*eye(n_s+n_p);
ukf_params.alpha=1;
ukf_params.beta=0;
ukf_params.kappa=0;%3-(n_s+n_p);
[x_hat_f,P_xx_f,x_hat_s,P_xx_s,~, ~]=uks_em_nmm_gen(ukf_params,f_nmm,h,x_hat0,P_xx0, y, H, Q0, R0, p, list(j,:), n_s);

%SAVE RES
save(sprintf('/x_f_%d.mat',j), 'x_hat_f')
save(sprintf('x_s_%d.mat',j), 'x_hat_s')
save(sprintf('P_f_%d.mat',j), 'P_xx_f')
save(sprintf('P_s_%d.mat',j), 'P_xx_s')