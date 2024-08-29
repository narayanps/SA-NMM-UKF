function [x_hat_f,P_xx_f,x_hat_s,P_xx_s,Q_hist, loglikes]=uks_em_nmm_gen(ukf_params,f,h,x_hat0,P_xx0, y, H, Q, R, p, unknown_param_id, ns)

% ------------------------------------------------------------------------------
% Author: Narayan P. Subramaniyam
% Affiliation: MET Faculty, Tampere University
% Email: [narayan.subramaniyam@tuni.fi]
%
% Description:
% This MATLAB code is developed as part of my research. Feel free to reuse 
% or modify this code, provided that you give proper attribution by citing 
% the associated paper. 
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
N_ignore=5000;
for itr=1:10
    Q_hist(1,:,:) = Q;
   % B(:,:,1)=Q;
    R_hist(1) = R;
    loglike=0;
    x_hat=x_hat0;
    P_xx=P_xx0;
    
    n=length(x_hat); %state dimension
    m=2*n + 1; %no. of sigma-points
    k=size(y,1); %measurement dimension
    T=size(y,2); %no. of time-points
    
    alpha = ukf_params.alpha;
    beta = ukf_params.beta;
    kappa = ukf_params.kappa;
    lambda=alpha^2*(n+kappa)-n;
    
    %initialize
    %predicted mean and cov
    x_hat1 = zeros(n,T);
    P_xx1 = zeros(n,n,T);
    
    %filtered mean and cov
    x_hat_f = zeros(n,T);
    P_xx_f = zeros(n,n,T);
    
    %smoothed mean and cov
    x_hat_s = zeros(n,T);
    P_xx_s = zeros(n,n,T);
    
    
    %set weights
    [W_mean, W_cov] = set_weights(lambda,alpha,beta,n);
    
    
    %UKF
    for t=1:1:T
%           if itr>1
%               Q=B(:,:,t);
%               Q=Q.*eye(n);
%           end
        %generate sigma-points
        [X] = generate_sp(lambda, n, x_hat, P_xx);
        %X(end,:)=min(p_max, max(p_min, X(end,:)));
        E_xf=X-x_hat(:,ones(1,m));
        %time-update
        
        [x_hat, y_hat, P_xx, P_yy, P_xy, E_x1] = ukf_time_upd(f, h, W_mean, W_cov, X, n, m, k ,H, Q, R, p(t),unknown_param_id, ns);
        
        %store predicted value
        x_hat1(:,t) = x_hat;
        P_xx1(:,:,t) = P_xx;
        C(:,:,t) = E_xf*diag(W_cov)*E_x1';
        
        %measurement update
        [x_hat, P_xx] = ukf_meas_upd(y(:,t), x_hat, y_hat, P_xx, P_yy, P_xy);
        
        %store filtered value
        x_hat_f(:,t) = x_hat;
        P_xx_f(:,:,t) = P_xx;
        loglike = loglike + compute_loglike(y(:,t),y_hat,P_yy);
    end
    loglikes(itr) = loglike;
    %UKS
    x_hat_s(:,end) = x_hat_f(:,end);
    P_xx_s(:,:,end) = P_xx_f(:,:,end);
    
    for t=T-1:-1:1
        D(:,:,t) = C(:,:,t+1)/P_xx1(:,:,t+1);
        x_hat_s(:,t) = x_hat_f(:,t) + D(:,:,t)*(x_hat_s(:,t+1) - x_hat1(:,t+1));
        P_xx_s(:,:,t) = P_xx_f(:,:,t) + D(:,:,t)*(P_xx_s(:,:,t+1) - P_xx1(:,:,t+1))*D(:,:,t)';
        z_hat_s(:,t) = [x_hat_s(:,t) ; x_hat_s(:,t+1)];
        S(:,:,t) = [P_xx_s(:,:,t) D(:,:,t)*P_xx_s(:,:,t+1); P_xx_s(:,:,t+1)*D(:,:,t)' P_xx_s(:,:,t+1)];
    end
     
     for t=1+N_ignore:T-1
        
        [xi(:,:,t)] = generate_sp(lambda, 2*n, z_hat_s(:,t), squeeze(S(:,:,t)));
    end
    [W_mean_s, W_cov_s]=set_weights(lambda,alpha,beta,2*n);

B=0;
ETA=0;

for t=1+N_ignore:T-1
    P_ff=0;
    P_xf=0;
    f_hat=0;
    for i=1:1:size(xi,2)
        fs = f(xi(1:n,i,t), p(t),unknown_param_id, ns);
        f_hat = f_hat + (W_mean_s(i)*fs);
    end
    U = (x_hat_s(:,t+1)-f_hat)*(x_hat_s(:,t+1)-f_hat)';
    for i=1:1:size(xi,2)
        fs = f(xi(1:n,i,t),p(t),unknown_param_id, ns);
        P_xf = P_xf + W_cov_s(i)*(xi(n+1:end,i,t) - x_hat_s(:,t+1))*(fs-f_hat)';
        P_ff = P_ff  + W_cov_s(i)*(fs-f_hat)*(fs-f_hat)';
    end
    ETA=ETA + (y(:,t)-H*x_hat_s(:,t))*(y(:,t)-H*x_hat_s(:,t))' + H*P_xx_s(:,:,t)*H';
   % B(:,:,t+1) = P_xx_s(:,:,t+1) + P_ff - P_xf - P_xf' + U;
    B = B + P_xx_s(:,:,t+1) + P_ff - P_xf - P_xf' + U;
end

Q=B./T;
%Q(1:n_s, 1:n_s)=1e-10.*eye(n_s);
Q=Q.*eye(n);
R = ETA./(T+1) %1; %(GAMMA - H*LAMBDA' - LAMBDA*H' + H*PSI*H')./T;
%mean(x_hat_s(7,:),2)
%mean(x_hat_s(8,:),2)
%     if sum(eig(Q) < 0) > 0
%         disp('hej')
%         Q = sqrtm(Q*Q');
%     end
    
    
    if itr>1
       abs(loglikes(itr)-loglikes(itr-1))
    end
    if itr > 1 && abs(loglikes(itr)-loglikes(itr-1)) < 1e-2
        abs(loglikes(itr)-loglikes(itr-1))
        return;
    end
end
function [W_mean, W_cov]=set_weights(lambda,alpha,beta,n)

W_mean=[lambda/(n+lambda) 0.5./((n+lambda)+zeros(1,2*n))];           %weights for meansfigure('DefaultAxesFontSize',18

W_cov=[(lambda/(n+lambda) + (1-alpha.^2 + beta)) 0.5./((n+lambda)+zeros(1,2*n))]; %weights for covariance

function [X] = generate_sp(lambda, n, x_hat, P_xx)
eta=sqrt(lambda+n);
[U,S,~] = svd(P_xx);
L = eta*U*sqrtm(S);
temp_x = repmat(x_hat, [1 n]);
X = [x_hat temp_x+L temp_x-L];
clear temp_x

function [x, y, Pxx, Pyy, Pxy, Ex] = ukf_time_upd(f, h, W_mean, W_cov, X, n, m, k, H, Q, R, p, unknown_param_id, ns)
x=zeros(n,1);
y=zeros(k,1);
X1=zeros(n,m);
Eta1=zeros(k,m);
for i=1:m
    X1(:,i)= f(X(:,i), p,unknown_param_id, ns);
     Eta1(:,i)= h(H, X1(:,i));
     x=x+W_mean(i)*X1(:,i);
     y=y+W_mean(i)*Eta1(:,i);
end
Ex=X1-x(:,ones(1,m));
Ey=Eta1-y(:,ones(1,m));
Pxx=Ex*diag(W_cov)*Ex'+ Q;   
% Pyy=Ey*diag(W_cov)*Ey'+R;   
% Pxy =Ex*diag(W_cov)*Ey';
Pyy=H*Pxx*H'+R; %Ey*diag(W_cov)*Ey'+R;   
Pxy = Pxx*H'; %Ex*diag(W_cov)*Ey';

function [x, P_xx] = ukf_meas_upd(y, xp, yp, Pp_xx, Pp_yy, Pp_xy)
K = Pp_xy*inv(Pp_yy);
x = xp + K*(y-yp);
P_xx = Pp_xx - K*Pp_yy*K';


function ll = compute_loglike(y,mu, sigma)
m=size(y,1);
CONST = -m * 0.5 * log(2*pi);
logDetSigma = 2*sum(log(diag(chol(sigma))));
ll = CONST - 0.5*logDetSigma - 0.5*((y-mu)'/sigma * (y-mu));







