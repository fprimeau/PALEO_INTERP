function S = run9b(fname,A3D,modeflag,rbf_type)

    tmp = load('myMTMJan2015.mat','M3d','grd','MTM');
    msk = tmp.M3d; % OCIM wet-dry mask (wet == 1, dry == 0);
    grd = tmp.grd; % OCIM mesh size  
    MTM = tmp.MTM; % OCIM tools to construct diffusion operator
    %% load the sediment core C-14 data                                %
    if (strcmp(modeflag,'unique'))
        core_data = load_unique_core_data(['unique_',fname]);
    elseif(strcmp(modeflag,'all'))
        core_data = load_unique_core_data(['all_',fname]);
    end        
    j = core_data(:,4); % indices into iwet for OCIM grid

    % load in the precomputed diffusive irfs
    tmp = load('OCIM_DIFFU_DIST1.mat');
    DE = tmp.DE;
    J = tmp.J;
    crit = tmp.crit;
    logad = tmp.logad;
    dt = tmp.dt;
    
    %% interplate the anomalies                                        %
    d = core_data(:,8);     % age anomaly (paleo - glodap)
    N = length(d);          % number of data points
    
    
    % standardize the data
    mean_d = mean(d);
    var_d = var(d);
    % standardize the data
    d = (d-mean_d)/sqrt(var_d);

    % make a function to undo the standardization
    undo_standardization = @(x) x*sqrt(var_d)+mean_d;
    % precision matrix for c14-age anomalies
    a_err = core_data(:,9); % uncertainty
    a_err = a_err/sqrt(var_d);
    ww = 1./a_err.^2; ww = ww/min(ww);
    ww = ww/max(ww);
    WW = spdiags(ww,0,N,N); 
    
    %%%
    iwet = find(msk(:));
    VT =  sum(grd.dVt(iwet));
    
    %% Check if j is in J. if not compute new columns for DE
    m = ismember(j,J);
    k = find(m==0);
    jnew = j(find(~ismember(j,J)));
    
    if (~isempty(jnew))
        % basis functions are build using a "distance" defined as the
        % square root of the diffusive time where the time is given as
        % the first time, t, when the tracer concentration at a given
        % point reaches a crit*100 percent of its final equilibrium
        % value.  crit is a number strictly between 0 and 1. 
        %
        % To construct the distance metric, inject a unit amount at
        % each grid point  where we have observations  and let these 
        % impulses diffuse through out the ocean. For tracer and for
        % each grid point on the mesh record the first time when the
        % concentration passes crit/VT where VT is the total volume of
        % the ocean.  
        %
        i = find(ismember(j,jnew)); % find the rows of core_data corresponding to jnew

        ix = core_data(i,1);
        iy = core_data(i,2);
        iz = core_data(i,3);
        jnew = core_data(i,4); % this should NOT chnage the value of jnew

        STUFF = diffu_setup_noflux(jnew,dt,logad,MTM,grd,msk,ix,iy,iz);
        Ro = STUFF.R;   % initial impulses for constructing the distance
        R = Ro;
        D = 0*Ro+2000;
        spa = 365*24*60^2;
        for i = 1:30000  % time-step loop
            R = mfactor(STUFF.FB,STUFF.AA*R); % implicit time-step
            c = find(abs(R(:))>crit/VT);
            D(c) = min(D(c),i*STUFF.dt/spa);

            % ileft is the number of grid boxes that have not yet been assigned a distance
            
            ileft = find(D(:)==2000);
            fprintf('%i boxes have not reached cirtical level yet\n',length(ileft));

            if (length(ileft)==0)
                break
            end        
        end
        DE = [DE,D.^(1/2)-sqrt(dt/spa)];
        J = [J;jnew];
        % make a backup of the old distance file
        eval(['!mv OCIM_DIFFU_DIST1.mat OCIM_DIFFU_DIST1_',date,'.mat']);
        % save the new distance file that now includes idffusive irfs for the new data locations
        save('OCIM_DIFFU_DIST1.mat','DE','J','crit','logad','dt');
    end

    
    %% optimize the level 3 parameters                                 %T

    
    
    switch (rbf_type)
      case {'gaussian'}
        rbf = @(ep,x) exp(-(ep*x).^2);
        mm = 2;
        ep = logspace(0,6,150);
      case {'inverse_quadric'}
        rbf = @(ep,x) 1./sqrt(1+(ep*x).^2); 
        mm = 1;
        ep = logspace(-0.5,2,200);
      case {'quadric'}
        rbf = @(ep,x) sqrt(1+(ep*x).^2);
        mm = 10;
        ep = logspace(10,20,100);
      case{'exp'}
        rbf = @(ep,x) exp(-ep*x); 
        mm = 2
        ep = logspace(1,6,100);
      case{'exp_linear'}
        rbf = @(ep,x) exp(-ep*x).*(1+ep*x);
        mm = 40;
        ep = logspace(45,48,100);
      case{'exp_cubic'}
        rbf = @(ep,x) exp(-ep*x).*(3+3*ep*x+(ep*x).^2); 
        mm = 80;
        ep = logspace(90,94,100);
      case{'distance'}
        rbf = @(ep,x) ep*x;
        mm = 2;
        ep = logspace(-1,1,100);
    end

    % re-load in the distance matrix (potentiall updated)
    load('OCIM_DIFFU_DIST1.mat','DE','J');
    kk = find(ismember(J,j));
    DE = DE(:,kk); % distance matrix to evaluate the interpolation at all wet OCIM points
    DI = DE(j,:);  % distance matrix to evaluate the interpolation where we have data
    
    % different values of mm can be tested to improve the interpolation    
    DE = DE.^mm;
    DI = DI.^mm;
    %
    mx = max(DE(:));
    DE = DE/mx;
    DI = DI/mx;

    %
    inflag = 0;
    [logZmx,SMX] = xmod8b(rbf,a_err,log10(ep(1)),inflag,d,msk,DI,DE,WW,1,1);
    fprintf('-------------------------------------------------------\n');
    epmx = ep(1);
    S = SMX;
    % do a grid search for the optimal scaling paramter (ep)
    for k = 1:length(ep)
        [logZ(k),S] = xmod8b(rbf,a_err,log10(ep(k)),inflag,d,msk,DI,DE,WW,S.alpha,S.beta);
        ALPHA(k) = S.alpha;
        BETA(k) = S.beta;
        if (logZ(k)>logZmx)
            logZmx = logZ(k);
            epmx = ep(k);
            SMX = S;
        end
        %fprintf('alpha = %e beta = %e logP = %e\n',S.alpha,S.beta,logZ(k));

    end
    % once we have the optimal ep value do the final interpolation
    inflag = 1;
    [LOGZ,SMX] = xmod8b(rbf,a_err,log10(epmx),inflag,d,msk,DI,DE,WW,SMX.alpha,SMX.beta);

    % save the interpolation
    S = SMX;
    S.anom = undo_standardization(S.anom);

    
    % do a diagnostic figure
    set(0,'CurrentFigure',3);
    hold on
    zdef = sprintf('Z = ( d - %4.0f yrs )/%4.0f yrs',mean_d,sqrt(var_d));
    text(-2,-4,zdef,'FontSize',16);
    title(sprintf('%s standardized anomalies',fname(1:end-4)),'Interpreter','none');
    hold off
    set(0,'CurrentFigure',1)
    clf
    subplot(3,1,1);
    plot(log10(ep),logZ,'ob-'); hold on; plot(log10(epmx),SMX.logZ,'or'); 
    ylabel('logZ');drawnow
    subplot(3,1,2);
    plot(log10(ep),ALPHA,'ob-'); hold on; plot(log10(epmx),SMX.alpha,'or');
    ylabel('alpha');
    subplot(3,1,3);
    plot(log10(ep),BETA,'ob-'); hold on; plot(log10(epmx),SMX.beta,'or');
    xlabel('log10(ep)'); ylabel('beta');
    drawnow

    
    %% quatities for which we want error bars                          %
    L = chol(S.C,'lower');

    %global average
    W = grd.dVt(iwet)/VT;
    wem = W'*S.EM;
    avg_anom = undo_standardization(wem*S.w); % global averaged of c-14 age anomaly

    tmp = ((wem)*L);
    std_avg_anom = sqrt(tmp*tmp')*sqrt(var_d); % standard deviation of global average c-14 age anomaly
                                               
    %fprintf('global average of age anomaly %f +/- %f\n',avg_anom,std_avg_anom);

    
    avg_age = W'*(S.anom(iwet)+A3D(iwet));
    std_avg_age = std_avg_anom;
    
    %fprintf('global average of age %f +/- %f\n',avg_age,std_avg_age);

    mask = msk;
    mask(:,:,3:end) = 0;
    WR = mask(iwet).*grd.dVt(iwet);
    WR = WR/sum(WR);
    res_age = WR'*(S.anom(iwet)+A3D(iwet));
    res_age_anom = WR'*(S.anom(iwet));
    wem = WR'*S.EM;
    tmp = ((wem)*L)/sqrt(var_d);
    std_res = sqrt(tmp*tmp');
    S.res_age = res_age;
    S.std_res = std_res;
    %fprintf('average surface reservoir age %f +/- %f\n',res_age,std_res);
    %fprintf('average surface reservoir age anomaly %f +/- %f\n',res_age_anom,std_res);

    S.avg_anom = avg_anom;          % 
    S.std_avg_anom = std_avg_anom;
    S.avg_age = avg_age;
    S.std_avg_age = std_avg_age;
    S.glodap_age = A3D;
    S.age = A3D+S.anom;
    S.adat = undo_standardization(d);
    S.amod = undo_standardization(S.IM*S.w);    
    S.dat = d+S.glodap_age(iwet(j)); % data  (total age)
    S.mod = undo_standardization(S.IM*S.w)+S.glodap_age(iwet(j)); % model (total age)
    S.WW = WW;
    S.iwet = iwet;
    S.fname = fname;
    S.core_data = core_data;
    set(0,'CurrentFigure',1);
    eval(['print -depsc ',fname,modeflag,'_',rbf_type,'_logZVep.eps']);
    set(0,'CurrentFigure',3);
    eval(['print -depsc ',fname,modeflag,'_',rbf_type,'_modVobs.eps']);
    
end




function [logZ,S] = xmod8b(rbf,a_err,logep,inflag,d,msk,DM_data,DM_eval,WW,varargin)
    persistent alpha beta
    if (nargin==11)
        alpha = varargin{1};
        beta = varargin{2};
    end
    
    %% Prepare the evaluation and interpolation matrices               %
    ep = 10^logep;    
    IM = rbf(ep,DM_data); % Interpolation matrix
    EM = rbf(ep,DM_eval); % Evaluation matrix

    %%                                                                 %
    
    %% log likelihood  = log(prob(d|w,beta))                           %
    N = length(d);
    B = IM'*WW*IM; [U,SV,V] = svd(B); lamb = diag(SV);
    Ed = @(w) 0.5*(d-IM*w)'*WW*(d-IM*w);
    %%                                                                 %

    %% log prior  = log(prob(w|alpha))                                 %
    nw = size(IM,2);
    C = eye(nw);
    Ew = @(w) 0.5*w'*C*w;
    %%                                                                 %

    %% neg log posterior = log(prob(w|d,alpha,beta))                   %
    A = @(alpha,beta) alpha*C+beta*B;
    M = @(w,alpha,beta) alpha*Ew(w)+beta*Ed(w);
    w_hat = @(alpha,beta) A(alpha,beta)\(beta*IM'*WW*d);
    %%                                                                 %

    
    logZ_function = @(alpha, beta) -M( w_hat(alpha,beta), alpha,beta)...
        - 0.5*logdet(A(alpha,beta)) ...
        + 0.5*N*log(beta) + 0.5*(nw-1)*log(alpha);
    
    
    %% Find alpha and beta such that maximize the evidence             %
    %  2*alpha * Ew(w_hat(alpha,beta)) = gam(alpha,beta)
    %  2*beta * Ed(w_hat(alpha,beta))  = N-gam(alpha,beta)
    %
    if (isempty(alpha)) 
        alpha = 0;
    end
    if (isempty(beta)) 
        beta = 1;
    end
    wmp = w_hat(alpha,beta); % compute the most probable weights
    e = 1; tol = 1e-6; outflag = 1;
    it = 0;

    % do 15 re-estimation iterations before going to fminsearch 
    % to ensure that fminsearch converges
    % if (1)
    %    while ( (e>tol) & (it<15) )
    %        it = it+1;
    %        logP = logZ_function(alpha,beta);
    %
    %        lama = beta*lamb;
    %        gam = sum(lama./(lama+alpha));
    %        g = [2*alpha*Ew(wmp)-gam;...
    %             2*beta*Ed(wmp)-N+gam];
    %        e = sqrt(g'*g);
    %
    %        % use re-estimation to optimize alpha and beta
    %        alpha = gam/(2*Ew(wmp));
    %        beta = (N-gam)/(2*Ed(wmp));
    %        wmp = w_hat(alpha,beta);
    %        
    %        if (it == 1000)
    %            outflag = 0;
    %        end
    %        if (isinf(Ed(wmp))|isinf(logP)|isnan(logP))
    %            outflag = 0;
    %            fprintf('WTF?\n'); % something went wrong
    %            break
    %        end
    %    end
    %end
    F = @(p) -logZ_function(p(1),p(2));
    options = optimset('MaxIter',1000);
    options = optimset('TolX',1e-6);
    options = optimset('TolFun',1e-6);
    % double check that we have the most probable alpha and beta
    [phat,fval,exitflag] = fminsearch(F,[alpha;beta],options);

    alpha = phat(1);
    beta = phat(2);
    logP = logZ_function(alpha,beta);
    wmp = w_hat(alpha,beta);
    
    % make a diagnostic figure 
    % plot the standardized model against the standardized data
    set(0,'CurrentFigure',3)
    errorbar(d,IM*wmp,a_err,'d');    
    hold on
    plot([-5 5],[-5 5],'--k');
    axis([-5 5 -5 5]); axis square;
    set(gca,'XTick',[-5:1:5]);
    set(gca,'YTick',[-5:1:5]);
    xlabel('Z_{obs}'); ylabel('Z_{mod}');
    text(-4,4.5,sprintf('alpha = %6.2e beta = %6.2e',alpha,beta),'FontSize',13);
    text(-4,3.8,['log10(ep) = ',num2str(logep)],'FontSize',13)
    text(-4,2.8,['Ed = ',num2str(Ed(wmp))],'FontSize',13);
    text(-4,1.8,['logZ = ',num2str(logP)],'FontSize',13);
    set(gca,'FontSize',16);
    grid on
    drawnow
    hold off
    

    % do the interpolation to the full OCIM grid
    iwet = find(msk(:));
    S.anom = msk+nan;
    S.anom(iwet) = EM*wmp;      % age anomaly (lgm-glodap)
    S.wmp = wmp;
    
    %%                                                                 %
    if ((outflag == 1)&(exitflag==1)) % everything is OK alpha and beta converged to their most probable value
        % uncertainty of alpha and beta
        %dalpha = sqrt(2/gam);
        %dbeta = sqrt(2/(N-gam));
        logZ = logP;

        S.logZ = logZ;
        S.alpha = alpha;          % w precision in prior
        S.beta = beta;            % data precision in likelihood
    
        if (inflag == 1) %  do the interpolation onto the OCIM grid
            iwet = find(msk(:));
            S.anom = msk+nan;
            S.anom(iwet) = EM*wmp;      % age anomaly (lgm-glodap)
            S.C = nearestSPD(inv(A(alpha,beta))); % w*w' covariance matrix
            S.w = wmp;                % interpolation weights
            S.alpha = alpha;          % w precision in prior
            S.beta = beta;            % data precision in likelihood
            S.ep = 10.^logep;         % basis function shape parameter
            S.logZ = logZ;            % log of evidence(ep)
            S.EM = EM;
            S.IM = IM;
        end
    else % alpha and beta did not converge to their most probable values
        logZ = -Inf;
        S.alpha = 0;
        S.beta = 1;
    end
end

function DM = DistanceMatrix(dsites,ctrs)
% DM = DistanceMatrix(dsites,ctrs)
% Forms the distance matrix of two sets of points in R^s,
% i.e., DM(i,j) = || datasite_i - center_j ||_2.
% Input
%   dsites: Mxs matrix representing a set of M data sites in R^s
%              (i.e., each row contains one s-dimensional point)
%   ctrs:   Nxs matrix representing a set of N centers in R^s
%              (one center per row)
% Output
%   DM:     MxN matrix whose i,j position contains the Euclidean
%              distance between the i-th data site and j-th center
    
    [M,s] = size(dsites); [N,s] = size(ctrs);
    DM = zeros(M,N);
    % Accumulate sum of squares of coordinate differences
    % The ndgrid command produces two MxN matrices:
    %   dr, consisting of N identical columns (each containing
    %       the d-th coordinate of the M data sites)
    %   cc, consisting of M identical rows (each containing
    %       the d-th coordinate of the N centers)
    for d=1:s
        [dr,cc] = ndgrid(dsites(:,d),ctrs(:,d));
        DM = DM + (dr-cc).^2;
    end
    DM = sqrt(DM);
end


