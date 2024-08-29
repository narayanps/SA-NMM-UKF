function [x_hat_f,P_xx_f,x_hat_s,P_xx_s] =  test_jr_real_data_human()
close all

addpath(strcat(pwd, '\func\'))
load('Sz1.mat')
Fs=512;
Fs_new = 512;
%EEG = resample(EEG, Fs_new,Fs);
start_ind = 1*Fs_new*60;
end_ind = length(EEG)-1*Fs_new*60;
ch=27;
param_comb = [2 6]; %B and b
y=EEG(start_ind:end_ind,ch);
T=length(y);
n_s=6;
n_p=2;
p = 90 + 30.*randn(1,T);
H=zeros(1,n_s+n_p);
H(2)=1; H(3)=-1; 
f_nmm=@nmm_jr_ukf_gen;
h=@h_meas;
Q0=1e-10*eye(n_s+n_p);
Q0(7,7)=1e-2;
Q0(8,8)=1e-2;
R0=1;
x_hat0=zeros(n_s+n_p,1);
x_hat0(7,1) = 70;
x_hat0(8,1) = 50;
P_xx0=eye(n_s+n_p);
ukf_params.alpha=1;
ukf_params.beta=0;
ukf_params.kappa=0;%3-(n_s+n_p);
[x_hat_f,P_xx_f,x_hat_s,P_xx_s,~, ~]=uks_em_nmm_gen(ukf_params,f_nmm,h,x_hat0,P_xx0, y', H, Q0, R0, p, param_comb, n_s);

end