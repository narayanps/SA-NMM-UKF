function res=sobol_nmm(f,dim,n,varargin)
    % Find Sobol's indices (<a href="matlab:a=fileparts(which('CODES.install'));file=strcat(a,'/+doc/html/sobol.html');web(file);">HTML</a>)
    %
    % Syntax
    %   res=CODES.sensitivity.sobol(f,dim,n) compute first order,
    %   second order and total global sensitivity indices S, Sij and St
    %   respectively of a function f. The problem dimensions dim and sample
    %   size n must be provided. Note that this function will call f
    %   n(2*dim+2) times.
    %   [...]=CODES.sensitivity.sobol(...,param,value) uses a list of
    %   parameter and value, please refer to the <a
    %   href="matlab:a=fileparts(which('CODES.install'));file=strcat(a,'/+doc/html/sobol.html');web(file);">HTML</a>
    %   documentation
    %
    % Example
    %   f=@(x)1/8*prod(3*x.^2+1,2);
    %   dim=3;
    %   n=1e3;
    %   res=CODES.sensitivity.sobol(f,dim,n);
    %
    % See also
    % CODES.sampling.anti_lock, CODES.sampling.edsd, CODES.sampling.gmm
    %
    % Copyright 2013-2015 Computational Optimal Design of Engineering
    % Systems (CODES) laboratory
    
    % The original code belongs to CODES laboratory. I have modified it to
    % work with NMM ODEs.
    % ------------------------------------------------------------------------------
% Additional Author: Narayan P. Subramaniyam
% Affiliation: MET Faculty, Tampere University
% 
%
% Description:
% This MATLAB code is developed as part of my research. Feel free to reuse 
% or modify this code, provided that you give proper attribution by citing 
% the associated paper and CODES laboratory package
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

    
    H=[0 1 -1 0 0 0];
    T_tot=100000;
    ns=6;
    input=inputParser;
    input.KeepUnmatched=false;
    input.PartialMatching=false;
    input.addRequired('f',@(x)isa(x,'function_handle'));
    input.addRequired('dim',@isnumeric);                % Problem dimension
    input.addRequired('n',@isnumeric);                  % Number of samples to be used (x2)
    input.addOptional('lb',[],@isnumeric);              % Lower bound for side constraints
    input.addOptional('ub',[],@isnumeric);              % Upper bound for side constraints
    input.addOptional('IDF',[],...
        @(x)isa(x,'function_handle'));                  % Marginal inverse distribution functions
    input.addOptional('sampler','rand',@(x)...
        strcmp(x,'rand')||...
        strcmp(x,'halton')||...
        strcmp(x,'sobol')||...
        strcmp(x,'lhs')||...
        strcmp(x,'cvt'));                               % Sampler type
    input.addOptional('vectorized',false,@islogical);   % Is function f vectorized
    input.addOptional('f_parallel',false,@islogical);   % When f is not vectorized, should f be evaluated in parallel
    input.addOptional('conv_seq',[],@isnumeric);        % Convergence sequence for convergence plot
    input.addOptional('conv_leg',true,@islogical);      % Whether to add legend on convergence plot or not
    input.addOptional('bar_plot',false,@islogical);     % Whether to plot a bar plot or not
    input.addOptional('bar_leg',true,@islogical);       % Whether to add legend on bar plot or not
    input.addOptional('CI_boot',false,@islogical);      % Whether to compute bootstrap CI or not
    input.addOptional('nb_boot',200,@isnumeric);        % Number of bootstrap
    input.addOptional('alpha',0.05,@isnumeric);         % Significance level
    input.addOptional('boot_type','bca',@(x)...
        strcmp(x,'bca')||...
        strcmp(x,'norm')||...
        strcmp(x,'per')||...
        strcmp(x,'cper'));                              % Bootstrap CI type
    input.addOptional('err_plot',false,@islogical);     % Whether to plot a error plot or not (force CI_boot=true)
    input.parse(f,dim,n,varargin{:})
    in=input.Results;
    % Checks
    if ~any(strcmp(input.UsingDefaults,'lb')) && ~any(strcmp(input.UsingDefaults,'ub'))
        assert(~isfield(in,'IDF')||any(strcmp(input.UsingDefaults,'IDF')),'Bounds and IDF can''t be defined simultaneously.');
    end
    if ~any(strcmp(input.UsingDefaults,'IDF'))
        assert(~isfield(in,'lb')||any(strcmp(input.UsingDefaults,'lb')),'Bounds and IDF can''t be defined simultaneously.');
        assert(~isfield(in,'ub')||any(strcmp(input.UsingDefaults,'ub')),'Bounds and IDF can''t be defined simultaneously.');
    end
    if ~any(strcmp(input.UsingDefaults,'lb')) || ~any(strcmp(input.UsingDefaults,'ub'))
        assert(~any(strcmp(input.UsingDefaults,'lb')) && ~any(strcmp(input.UsingDefaults,'ub')),'Both bounds must be provided.');
    end
    if ~any(strcmp(input.UsingDefaults,'conv_seq'))
        assert(max(in.conv_seq)<=in.n,'max(conv_seq) must be lower than or equal to n');
    end
    if in.err_plot
        in.CI_boot=true;
    end
    % Create 2 DOE
    if nargin==1 || ischar(n)
        x_val=in.res.function_values.x_val;
        y_val=in.res.function_values.y_val;
        x_pert_val=in.res.function_values.x_pert_val;
        y_pert_val=in.res.function_values.y_pert_val;
    else
        f=@in.f;
        switch in.sampler
            case 'rand'
                X=rand(in.n,in.dim);
                Y=rand(in.n,in.dim);
            case 'lhs'
                X=lhsdesign(in.n,in.dim);
                Y=lhsdesign(in.n,in.dim);
            case 'halton'
                qr_set=scramble(haltonset(in.dim));
                qr_set1=scramble(haltonset(in.dim,'skip',in.n));
                X=net(qr_set,in.n);
                X=X(randperm(in.n),:);
                Y=net(qr_set1,in.n);
                Y=Y(randperm(in.n),:);
            case 'sobol'
                qr_set=scramble(sobolset(in.dim));
                qr_set1=scramble(sobolset(in.dim,'skip',in.n));
                X=net(qr_set,in.n);
                X=X(randperm(in.n),:);
                Y=net(qr_set1,in.n);
                Y=Y(randperm(in.n),:);
            case 'cvt'
                X=CODES.sampling.cvt(in.n,in.dim);
                Y=CODES.sampling.cvt(in.n,in.dim,'force_new',true,'display',false);
        end
        % Transform samples
        if ~any(strcmp(input.UsingDefaults,'IDF'))          % Transform samples using IDF
            X=in.IDF(X);
            Y=in.IDF(Y);
        elseif ~any(strcmp(input.UsingDefaults,'lb'))       % Use bounds
            range=in.ub-in.lb;
            X=bsxfun(@plus,bsxfun(@times,X,range),in.lb);
            Y=bsxfun(@plus,bsxfun(@times,Y,range),in.lb);
        end
        % Evaluate them
        if in.vectorized
            
            t=in.f(X(1,:));
            assert(size(t,1)==1,'Function f must return multiple outputs in columns')
            x_val=[t;in.f(X(2:end,:))];
            y_val=in.f(Y);
            in.dim_Y=size(x_val,2);
        else
            state=zeros(ns,1);
            x_val=[];
            for ts=1:T_tot
            [state]=in.f(X(1,:), state);
            t=H*state;
            in.dim_Y=1;%size(t,2);
            x_val=[x_val t];
            end
            state=zeros(ns,1);
            y_val=[];
            for ts=1:T_tot
            [state]=in.f(Y(1,:), state);
            t=H*state;
            assert(size(t,1)==1,'Function f must return multiple outputs in columns')
            in.dim_Y=1;%size(t,2);
            y_val=[y_val t];
            end
            in.f_parallel=1;
            if in.f_parallel
                parfor i=2:in.n
                    state=zeros(ns,1);
                    tmp=[];
                    for ts=1:T_tot
                    state=f(X(i,:), state);
                    tmp=[tmp; H*state];
                    end
                    x_val(i,:) = tmp;
                    tmp=[];
                    state=zeros(ns,1);
                    for ts=1:T_tot
                    state=f(Y(i,:), state);
                    tmp=[tmp; H*state];
                    end
                    y_val(i,:) = tmp;
                end
            else
                for i=2:in.n
                    state=zeros(ns,1);
                    for ts=1:T_tot
                        state=f(X(i,:), state);
                        x_val(i,ts)= H*state;
                    end
                    state=zeros(ns,1);
                    for ts=1:T_tot
                        state=f(Y(i,:), state);
                        y_val(i,ts)=H*state;
                    end
                end
            end
        end
        % Get 3D matrix diagonal indices
        idx=bsxfun(@plus,1:(in.dim+1)*in.n:in.n*in.dim^2,(0:in.n-1)');
        % Create the map of perturbated points
        x_pert=repmat(X,[1 1 in.dim]);
        x_pert(idx)=Y;
        y_pert=repmat(Y,[1,1 in.dim]);
        y_pert(idx)=X;
        % Unfold perturbated points
        x_pert=reshape(permute(x_pert,[1 3 2]),[in.n*in.dim in.dim]);
        y_pert=reshape(permute(y_pert,[1 3 2]),[in.n*in.dim in.dim]);
        % Evaluate them
        if in.vectorized
            x_pert_val=in.f(x_pert);
            y_pert_val=in.f(y_pert);
        else
            x_pert_val=zeros(in.n*in.dim,T_tot);
            y_pert_val=zeros(in.n*in.dim,T_tot);
            in.f_parallel=1;
            if in.f_parallel
                parfor i=1:in.n*in.dim
                    statex=zeros(ns,1);
                    statey=zeros(ns,1);
                    tmpx=[];
                    tmpy=[];
                    for ts=1:T_tot
                    statex=f(x_pert(i,:),statex);
                    statey=f(y_pert(i,:),statey);
                    tmpx = [tmpx ; H*statex];
                    tmpy = [tmpy ; H*statey];
                    end
                    y_pert_val(i,:) = tmpy;
                    x_pert_val(i,:) = tmpx;
                end
            else
                for i=1:in.n*in.dim
                    statex=zeros(ns,1);
                    statey=zeros(ns,1);
                    tmpx=[];
                    tmpy=[];
                    for ts=1:T_tot
                    statex=f(x_pert(i,:),statex);
                    statey=f(y_pert(i,:),statey);
                    tmpx = [tmpx ; H*statex];
                    tmpy = [tmpy ; H*statey];
                    end
                    y_pert_val(i,:) = tmpy;
                    x_pert_val(i,:) = tmpx;
                end
            end
        end
    end
    for t=1:T_tot
    x_val_t=permute(x_val(:,t),[1 3 2]);
    y_val_t=permute(y_val(:,t),[1 3 2]);
    x_pert_val_t=permute(x_pert_val(:,t),[1 3 2]);
    y_pert_val_t=permute(y_pert_val(:,t),[1 3 2]);
    [res_raw(t).S1,res_raw(t).S2,res_raw(t).St]=sobols(x_val_t,y_val_t,x_pert_val_t,y_pert_val_t,in,in.n);
    end
    save('res_sobol.mat', 'res_raw', '-v7.3');
    % Compute bootstrap Ci if requested
    if in.CI_boot
        [res_raw.S1_CI_boot,res_raw.S2_CI_boot,res_raw.St_CI_boot]=compute_CI(...
            x_val,y_val,x_pert_val,y_pert_val,in,'alpha',in.alpha,'type',in.boot_type);
    else
        res_raw.S1_CI_boot=[];
        res_raw.St_CI_boot=[];
        res_raw.S2_CI_boot=[];
    end
    % Construct ouput
    res=CODES.build.sobol_out(...
        {res_raw.S1,res_raw.S2,res_raw.St},...
        {res_raw.S1_CI_boot,res_raw.S2_CI_boot,res_raw.St_CI_boot},...
        @()plot_all_bars(res_raw.S1,res_raw.St,res_raw.S2,0,in),...
        @()err_plot(in,res_raw),...
        in.CI_boot,...
        @(seq)conv_plot(x_val,y_val,x_pert_val,y_pert_val,in,seq),...
        @(varargin)compute_CI(x_val,y_val,x_pert_val,y_pert_val,in,varargin));
    % Bar plot
    if in.bar_plot
        plot_all_bars(res_raw.S1,res_raw.St,res_raw.S2,0,in);
    end
    % Error plot
    if in.err_plot
        err_plot(in,res_raw)
    end
    % Convergence plot, if requested
    if ~any(strcmp(input.UsingDefaults,'conv_seq'))
        conv_plot(x_val,y_val,x_pert_val,y_pert_val,in,in.conv_seq)
    end
    % Nested functions
    function [Dx,Dy]=Di(x,y,x_pert,y_pert,mean_f,in,subset)
        % Compute the Di quantities on a subset of the DOE
        x_pert=reshape(x_pert,[size(x_pert,1)/in.dim in.dim in.dim_Y]);
        y_pert=reshape(y_pert,[size(y_pert,1)/in.dim in.dim in.dim_Y]);
        Dx=bsxfun(@minus,...
            mean(bsxfun(@times,...
                        x(1:subset,:,:),...
                        y_pert(1:subset,:,:)),...
            1),mean_f.^2);
        Dy=bsxfun(@minus,...
            mean(bsxfun(@times,...
                       y(1:subset,:,:),...
                       x_pert(1:subset,:,:)),...
            1),mean_f.^2);
    end
    function [Dx,Dy]=Dmi(x,y,x_pert,y_pert,mean_f,in,subset)
        % Compute the D-i quantities on a subset of the DOE
        x_pert=reshape(x_pert,[size(x_pert,1)/in.dim in.dim in.dim_Y]);
        y_pert=reshape(y_pert,[size(y_pert,1)/in.dim in.dim in.dim_Y]);
        Dx=bsxfun(@minus,...
            mean(bsxfun(@times,...
                       x(1:subset,:,:),...
                       x_pert(1:subset,:,:)),...
            1),mean_f.^2);
        Dy=bsxfun(@minus,...
            mean(bsxfun(@times,...
                       y(1:subset,:,:),...
                       y_pert(1:subset,:,:)),...
            1),mean_f.^2);
    end
    function Dxyc=Dij(Dxi,Dyi,x_pert,y_pert,mean_f,in,subset)
        % Compute close Dij on a subset of the DOE
        x_pert=reshape(x_pert,[size(x_pert,1)/in.dim in.dim in.dim_Y]);
        y_pert=reshape(y_pert,[size(y_pert,1)/in.dim in.dim in.dim_Y]);
        Dxyc=bsxfun(@minus,...
             mean(bsxfun(@times,...
                         permute(x_pert(1:subset,:,:),[4 2 3 1]),...
                         permute(y_pert(1:subset,:,:),[2 4 3 1])),...
             4),mean_f.^2);
        Dxyc=bsxfun(@minus,Dxyc,permute(Dxi,[2 1 3]));
        Dxyc=bsxfun(@minus,Dxyc,Dyi);
    end
    function [S1,S2,St]=sobols(x_val,y_val,x_pert_val,y_pert_val,in,subset)
        % Compute Sobol indices on a subset of the DOE
        power_mean=sqrt(mean(x_val(1:subset,:,:).*y_val(1:subset,:,:)));
        power_variance=var([x_val(1:subset,:,:);y_val(1:subset,:,:)]);
        [Dxi,Dyi]=Di(x_val,y_val,x_pert_val,y_pert_val,power_mean,in,subset);
        [Dxmi,Dymi]=Dmi(x_val,y_val,x_pert_val,y_pert_val,power_mean,in,subset);
        Dxyc=Dij(Dxi,Dyi,x_pert_val,y_pert_val,power_mean,in,subset);
        S1=abs(bsxfun(@rdivide,[Dxi;Dyi],power_variance));
        S2=abs(bsxfun(@rdivide,Dxyc,power_variance));
        for ii=1:in.dim_Y
            for jj=1:in.dim
                S2(jj,jj,ii)=nan;
            end
        end
        St=abs(1-bsxfun(@rdivide,[Dxmi;Dymi],power_variance));
    end
    function vals=for_boot(x_val,y_val,x_pert_val,y_pert_val,in,boot_idx)
        % Compute Sobol indices during a bootstrap loop
        vals=zeros(4+in.dim,in.dim,in.dim_Y);
        nb_boot=length(boot_idx);
        % Keep only bs values
        x_val=x_val(boot_idx,:,:);
        y_val=y_val(boot_idx,:,:);
        % Get indices for perturbated values
        boot_idx=reshape(bsxfun(@plus,...
                                boot_idx,...
                                (0:in.dim-1)*in.n),...
                         [nb_boot*in.dim 1]);
        % Keep only bs values
        x_pert_val=x_pert_val(boot_idx,:,:);
        y_pert_val=y_pert_val(boot_idx,:,:);
        % Conmpute indices
        power_mean=sqrt(mean(x_val.*y_val));
        power_variance=var([x_val;y_val]);
        [Dxi,Dyi]=Di(x_val,y_val,x_pert_val,y_pert_val,power_mean,in,nb_boot);
        [Dxmi,Dymi]=Dmi(x_val,y_val,x_pert_val,y_pert_val,power_mean,in,nb_boot);
        Dxyc=Dij(Dxi,Dyi,x_pert_val,y_pert_val,power_mean,in,nb_boot);
        vals(1:2,:,:)=abs(bsxfun(@rdivide,[Dxi;Dyi],power_variance));
        vals(5:end,:,:)=abs(bsxfun(@rdivide,Dxyc,power_variance));
        for ii=1:in.dim
            for jj=1:in.dim_Y
                vals(4+ii,ii,jj)=55;
            end
        end
        vals(3:4,:,:)=abs(1-bsxfun(@rdivide,[Dxmi;Dymi],power_variance));
    end
    function X=st(X)
        X=permute(X,[4 3 2 1]);
    end
    function slider_callback(hObject,~,S,St,Sij,in,ind)
        hold off
        [~,fig_h]=gcbo;
        plot_bars(ind,S,St,Sij,hObject.Value,in,fig_h);
    end
    function plot_all_bars(S,St,Sij,th,in)
        for ii=1:size(S,3)
            plot_bars(ii,S,St,Sij,th,in)
        end
    end
    function plot_bars(ind,S,St,Sij,th,in,handle)
        if nargin==6
            figure('Position',[200 200 500 500])
        else
            figure(handle);
        end
        d=size(S,2);
        Colors=hsv(d+sum(1:d-1));
        leg_all=cell(d+sum(1:d-1),1);
        plot([1 2 3],[0 0 0],'k.','HandleVisibility','off','Markersize',1)
        hold on
        height=zeros(2,3);
        S_displayed=ones(2,d+sum(1:d-1));
        kk=d+1;
        for ii=1:d
            if any(S(:,ii,ind)>th)
                fill([0.75 1 1 0.75],height(1,1)+[0 0 S(1,ii,ind) S(1,ii,ind)],Colors(ii,:))
                fill([1 1.25 1.25 1],height(2,1)+[0 0 S(2,ii,ind) S(2,ii,ind)],Colors(ii,:),'HandleVisibility','off')
                height(:,1)=height(:,1)+S(:,ii,ind);
            else
                S_displayed(1,ii)=0;
            end
            if any(St(:,ii,ind)>th)
                if S_displayed(1,ii)==1
                    fill([1.75 2 2 1.75],height(1,2)+[0 0 St(1,ii,ind) St(1,ii,ind)],Colors(ii,:),'HandleVisibility','off')
                else
                    fill([1.75 2 2 1.75],height(1,2)+[0 0 St(1,ii,ind) St(1,ii,ind)],Colors(ii,:))
                end
                fill([2 2.25 2.25 2],height(2,2)+[0 0 St(2,ii,ind) St(2,ii,ind)],Colors(ii,:),'HandleVisibility','off')
                height(:,2)=height(:,2)+St(:,ii,ind);
            else
                S_displayed(2,ii)=0;
            end
            leg_all{ii}=['$S_{' num2str(ii) '}$,$S_{' num2str(ii) '}^T$'];
        end
        for ii=1:d
            for jj=ii+1:d
                if Sij(ii,jj)>th || Sij(jj,ii)>th
                    fill([2.75 3 3 2.75],height(1,3)+[0 0 Sij(ii,jj,ind) Sij(ii,jj,ind)],Colors(kk,:))
                    fill([3 3.25 3.25 3],height(2,3)+[0 0 Sij(jj,ii,ind) Sij(jj,ii,ind)],Colors(kk,:),'HandleVisibility','off')
                    height(:,3)=height(:,3)+[Sij(ii,jj,ind);Sij(jj,ii,ind)];
                else
                    S_displayed(:,kk)=0;
                end
                leg_all{kk}=['$S_{' num2str(ii) ',' num2str(jj) '}$'];
                kk=kk+1;
            end
        end
        if in.bar_leg
            hh=legend(leg_all{any(S_displayed,1)},'location','best');
            set(hh,'interpreter','latex')
        end
        plot([0.5 3.5],[1 1],'k--')
        plot([0.5 3.5],[th th],'r--')
        set(gca,'xtick',[1 2 3],'xticklabel',{'$S_i$','$S_i^T$','$S_{ij}$'},'ticklabelinterpreter','latex')
        uicontrol('Style','text','Position',[5 10 90 20],...
        'String','Threshhold');
        uicontrol('Style','slider','Min',0,'Max',max(max(St(:,:,ind)))-0.001,...
            'Value',th,'Position', [100 10 350 20],'Callback',@(a,b)slider_callback(a,b,S,St,Sij,in,ind));
        uicontrol('Style','text','Position',[460 10 40 20],...
        'String',num2str(th,3));
        title(['Sobol indices for $Y_' num2str(ind) '$'],'interpreter','latex')
    end
    function err_plot(in,res_raw)
        for k=1:size(res_raw.S1,3)
            % First order sobol
            figure('Position',[200 200 500 500])
            hold on
            leg=cell(in.dim,1);
            for ii=1:in.dim
                set(gca,'ColorOrderIndex',ii)
                h=plot(ii,res_raw.S1(1,ii,k),'x');
                plot(ii,res_raw.S1(2,ii,k),'x','HandleVisibility','off','Color',get(h,'Color'));
                plot(ii+[-0.25 0],res_raw.S1_CI_boot(1,ii,k)*[1 1],'-','HandleVisibility','off','Color',get(h,'Color'));
                plot(ii+[-0.25 0],res_raw.S1_CI_boot(2,ii,k)*[1 1],'-','HandleVisibility','off','Color',get(h,'Color'));
                plot(ii+[0.25 0],res_raw.S1_CI_boot(3,ii,k)*[1 1],'-','HandleVisibility','off','Color',get(h,'Color'));
                plot(ii+[0.25 0],res_raw.S1_CI_boot(4,ii,k)*[1 1],'-','HandleVisibility','off','Color',get(h,'Color'));
                plot(ii+[0 0],[res_raw.S1_CI_boot(1,ii,k) res_raw.S1_CI_boot(2,ii,k)],'-','HandleVisibility','off','Color',get(h,'Color'));
                plot(ii+[0 0],[res_raw.S1_CI_boot(3,ii,k) res_raw.S1_CI_boot(4,ii,k)],'-','HandleVisibility','off','Color',get(h,'Color'));
                leg{ii}=['$S_{' num2str(ii) '}$'];
            end
            set(gca,'xtick',1:in.dim,'xticklabel',leg,'ticklabelinterpreter','latex')
            title(['First order indices of $Y_' num2str(k) '$'],'interpreter','latex')
            % Total indices
            figure('Position',[200 200 500 500])
            hold on
            leg=cell(in.dim,1);
            for ii=1:in.dim
                set(gca,'ColorOrderIndex',ii)
                h=plot(ii,res_raw.St(1,ii,k),'x');
                plot(ii,res_raw.St(2,ii,k),'x','HandleVisibility','off','Color',get(h,'Color'));
                plot(ii+[-0.25 0],res_raw.St_CI_boot(1,ii,k)*[1 1],'-','HandleVisibility','off','Color',get(h,'Color'));
                plot(ii+[-0.25 0],res_raw.St_CI_boot(2,ii,k)*[1 1],'-','HandleVisibility','off','Color',get(h,'Color'));
                plot(ii+[0.25 0],res_raw.St_CI_boot(3,ii,k)*[1 1],'-','HandleVisibility','off','Color',get(h,'Color'));
                plot(ii+[0.25 0],res_raw.St_CI_boot(4,ii,k)*[1 1],'-','HandleVisibility','off','Color',get(h,'Color'));
                plot(ii+[0 0],[res_raw.St_CI_boot(1,ii,k) res_raw.St_CI_boot(2,ii,k)],'-','HandleVisibility','off','Color',get(h,'Color'));
                plot(ii+[0 0],[res_raw.St_CI_boot(3,ii,k) res_raw.St_CI_boot(4,ii,k)],'-','HandleVisibility','off','Color',get(h,'Color'));
                leg{ii}=['$S_{' num2str(ii) '}^T$'];
            end
            set(gca,'xtick',1:in.dim,'xticklabel',leg,'ticklabelinterpreter','latex')
            title(['Total indices of $Y_' num2str(k) '$'],'interpreter','latex')
            % Total indices
            figure('Position',[200 200 500 500])
            hold on
            leg=cell(sum(1:in.dim-1),1);
            kk=1;
            for ii=1:in.dim
                for jj=ii+1:in.dim
                    set(gca,'ColorOrderIndex',kk)
                    h=plot(kk,res_raw.S2(ii,jj,k),'x');
                    plot(kk,res_raw.S2(jj,ii,k),'x','HandleVisibility','off','Color',get(h,'Color'));
                    plot(kk+[-0.25 0],res_raw.S2_CI_boot(ii,jj,k)*[1 1],'-','HandleVisibility','off','Color',get(h,'Color'));
                    plot(kk+[-0.25 0],res_raw.S2_CI_boot(ii+in.dim,jj,k)*[1 1],'-','HandleVisibility','off','Color',get(h,'Color'));
                    plot(kk+[0.25 0],res_raw.S2_CI_boot(jj,ii,k)*[1 1],'-','HandleVisibility','off','Color',get(h,'Color'));
                    plot(kk+[0.25 0],res_raw.S2_CI_boot(jj+in.dim,ii,k)*[1 1],'-','HandleVisibility','off','Color',get(h,'Color'));
                    plot(kk+[0 0],[res_raw.S2_CI_boot(ii,jj,k) res_raw.S2_CI_boot(ii+in.dim,jj,k)],'-','HandleVisibility','off','Color',get(h,'Color'));
                    plot(kk+[0 0],[res_raw.S2_CI_boot(jj,ii,k) res_raw.S2_CI_boot(jj+in.dim,ii,k)],'-','HandleVisibility','off','Color',get(h,'Color'));
                    leg{kk}=['$S_{' num2str(ii) ',' num2str(jj) '}$'];
                    kk=kk+1;
                end
            end
            set(gca,'xtick',1:in.dim,'xticklabel',leg,'ticklabelinterpreter','latex')
            title(['Second order indices for $Y_' num2str(k) '$'],'interpreter','latex')
        end
    end
    function conv_plot(x_val,y_val,x_pert_val,y_pert_val,in,conv_seq)
        assert(max(conv_seq)<=in.n,'max(conv_seq) must be lower than or equal to n');
        dim_Y=size(x_val,3);
        S_table=zeros(2,in.dim,dim_Y,length(conv_seq));
        St_table=zeros(2,in.dim,dim_Y,length(conv_seq));
        Sij_table=zeros(in.dim,in.dim,dim_Y,length(conv_seq));
        conv_seq=floor(conv_seq);
        for ii=1:length(conv_seq)
            % Compute indices on a subset
            [S_table(:,:,:,ii),Sij_table(:,:,:,ii),St_table(:,:,:,ii)]=sobols(x_val,y_val,x_pert_val,y_pert_val,in,conv_seq(ii));
        end
        for k=dim_Y
            % Create a figure
            % First order
            figure('Position',[200 200 500 500])
            if in.conv_leg
                leg=cell(in.dim,1);
            end
            hold on
            for ii=1:in.dim
                set(gca,'ColorOrderIndex',ii)
                h=plot(conv_seq,st(S_table(1,ii,k,:)),'o-');
                plot(conv_seq,st(S_table(2,ii,k,:)),'x-','HandleVisibility','off','Color',get(h,'Color'));
                if in.conv_leg
                    leg{ii}=['$S_{' num2str(ii) '}$'];
                end
            end
            if in.conv_leg
                h=legend(leg,'location','best');
                set(h,'interpreter','latex')
            end
            set(gca,'xscale','log')
            title(['First order indices for $Y_' num2str(k) '$'],'interpreter','latex')
            grid on
            % Total indices
            figure('Position',[200 200 500 500])
            if in.conv_leg
                leg=cell(in.dim,1);
            end
            hold on
            for ii=1:in.dim
                set(gca,'ColorOrderIndex',ii)
                h=plot(conv_seq,st(St_table(1,ii,k,:)),'o-');
                plot(conv_seq,st(St_table(2,ii,k,:)),'x-','HandleVisibility','off','Color',get(h,'Color'));
                if in.conv_leg
                    leg{ii}=['$S_{' num2str(ii) '}^T$'];
                end
            end
            if in.conv_leg
                h=legend(leg,'location','best');
                set(h,'interpreter','latex')
            end
            set(gca,'xscale','log')
            title(['Total indices for $Y_' num2str(k) '$'],'interpreter','latex')
            grid on
            % Second order
            figure('Position',[200 200 500 500])
            if in.conv_leg
                leg=cell(sum(1:in.dim-1),1);
            end
            hold on
            kk=1;
            for ii=1:in.dim
                for jj=ii+1:in.dim
                    set(gca,'ColorOrderIndex',kk)
                    h=plot(conv_seq,st(Sij_table(ii,jj,k,:)),'o-');
                    plot(conv_seq,st(Sij_table(jj,ii,k,:)),'x-','HandleVisibility','off','Color',get(h,'Color'));
                    if in.conv_leg
                        leg{kk}=['$S_{' num2str(ii) ',' num2str(jj) '}$'];
                        kk=kk+1;
                    end
                end
            end
            if in.conv_leg
                h=legend(leg,'location','best');
                set(h,'interpreter','latex')
            end
            set(gca,'xscale','log')
            title(['Second order indices for $Y_' num2str(k) '$'],'interpreter','latex')
            grid on
        end
    end
    % Compute CI_boot
    function [S1_CI,S2_CI,St_CI]=compute_CI(x_val,y_val,x_pert_val,y_pert_val,in,varargin)
        inp=inputParser;
        inp.KeepUnmatched=false;
        inp.PartialMatching=false;
        inp.addOptional('alpha',0.05,@isnumeric);         % Significance level
        inp.addOptional('type','bca',@(x)......
            strcmp(x,'bca')||...
            strcmp(x,'norm')||...
            strcmp(x,'per')||...
            strcmp(x,'cper'));                             % Bootstrap CI type
        while iscell(varargin{1})
            varargin=varargin{1};
        end
        inp.parse(varargin{:})
        inn=inp.Results;
        CIs=bootci(in.nb_boot,{@(x)for_boot(x_val,y_val,x_pert_val,y_pert_val,in,x),(1:in.n)'},'alpha',inn.alpha,'type',inn.type);
        for ii=1:in.dim
            for jj=1:in.dim_Y
                CIs(1:2,4+ii,ii,jj)=nan;
            end
        end
        S1_CI=permute([CIs(1:2,1,:,:);CIs(1:2,2,:,:)],[1 3 4 2]);
        St_CI=permute([CIs(1:2,3,:,:);CIs(1:2,4,:,:)],[1 3 4 2]);
        S2_CI=[permute(CIs(1,5:end,:,:),[2 3 4 1]);permute(CIs(2,5:end,:,:),[2 3 4 1])];
    end
end
