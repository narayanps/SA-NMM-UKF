function [Y, mi, sigma] =run_morris_analysis()
addpath(genpath(strcat(pwd, '/func/')))
addpath(genpath(strcat(pwd, '/ext/safe_R1.1/')))

M=13;
xmin = [0  0  0  1  25  6.5  0  0  0  0  2 0.5 0.3];
xmax = [10 50 200 1000 140 110 1350 1080 340 340 9 7.5 0.8];
DistrFun  = 'unif'  ;                                                          
DistrPar = cell(M,1);                                                          
for i=1:M 
    DistrPar{i} = [ xmin(i) xmax(i) ] ; 
end                             
%param_labels = {'A','B', 'mean', 'std', 'a', 'b', 'C1','C2','C3' ,'C4' } ;
r = 10000 ; % Number of Elementary Effects                                                                                                                                           
SampStrategy = 'lhs' ; 
design_type = 'radial';                                                        
X = OAT_sampling(r,M,DistrFun,DistrPar,SampStrategy,design_type);  
n_s=6; %number of states in NMM JR
T_ignore=10000;
T=100000;
T_tot=T+T_ignore;
state(:,1) = zeros(n_s, 1);
dt=0.001;
H=[0 1 -1 0 0 0 0 ];
h=@h_meas;
tic
parfor itr=1:length(X)
state = zeros(n_s, 1);
Y_tot=[];
for t=2:T_tot
    state= nmm_jr(X(itr,:)', state, dt);
    Y_tot = [Y_tot; h(H(1:n_s), state)];
end
Y(:,itr)=Y_tot(T_ignore+1:end);
end
toc
design_type='radial';
for t=1:size(Y,1)
[ mi(t,:), sigma(t,:) ] = EET_indices(r,xmin,xmax,X,Y(t,:)',design_type);
end
end