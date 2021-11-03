function S = run11(fname,A3D,modeflag,rbf_type,version)
%WTF    
    %% load the OCIM stuff                                              %
    tmp = load('myMTMJan2015.mat','M3d','grd','MTM');
    msk = tmp.M3d; % OCIM wet-dry mask (wet == 1, dry == 0);
    grd = tmp.grd; % OCIM mesh size  
    MTM = tmp.MTM; % OCIM tools to construct diffusion operator
    iwet = find(msk(:));     % indices into the wet OCIM gridboxes
    VT =  sum(grd.dVt(iwet));% volumes of the wet OCIM gridboxes


    %% load the sediment core C-14 data                                %
    if (strcmp(modeflag,'unique'))
        in = importdata(['unique_',fname]);
    elseif(strcmp(modeflag,'all'))
        in = importdata(['all_',fname]);
    end        

    j = in.data(:,7); % indices into iwet 
    a = in.data(:,9); % core c14-age minus glodap c14 age
    glodap = in.data(:,8); % glodap c14 age
    a_err = in.data(:,10); % core c14 age error estimate (assumed 1 s.d.)
    total_age = a + glodap;        % this is the total age (NOT just the anomaly)
    
    N = length(a);    

    % load in the precomputed diffusive irfs
    tmp = load(['OCIM_DIFFU_DIST1_',version,'.mat']);
    DE = tmp.DE;
    J = tmp.J;
    crit = tmp.crit;
    logad = tmp.logad;
    dt = tmp.dt;
    
    % log-transform the age and interpolate the log-transformed age to ensure a positive age interpolation
    d = log(total_age);

    % precision matrix for log c14-age
    sig = log(1+a_err./total_age);
    sig2 = sig.^2;

    ww = 1./sig2; ww = ww/min(ww);
    ww = ww/max(ww);
    WW = spdiags(ww,0,N,N); 
    
    %%%
    
    %% Check if j is in J. if not issue error message                  %
    m = ismember(j,J);
    k = find(m==0);
    jnew = j(find(~ismember(j,J)));
    if (~isempty(jnew))
        fprintf('You are missing the diffusive length from one of the core locations...\n');
        keyboard
    end
    
    %% optimize the level 3 parameters                                  %
       
    switch (rbf_type)
      case {'gaussian'}
        rbf = @(ep,x) exp(-(ep*x).^2);
        mm = 1;
        ep = logspace(-1,6,150);
      case {'inverse_quadric'}
        rbf = @(ep,x) 1./sqrt(1+(ep*x).^2); 
        mm = 1;
        ep = logspace(0,2,20); %[logspace(-0.5,0,10),fliplr(logspace(-1.0,0,10))];
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
    glodap = abs(A3D(iwet));
    prior_target = log(glodap);

    
    [logZmx,SMX] = xmod10(rbf,sig,log10(ep(1)),inflag,d,prior_target,msk,DI,DE,WW,1,1);
    fprintf('-------------------------------------------------------\n');
    epmx = ep(1);
    S = SMX;
    % do a grid search for the optimal scaling paramter (ep)
    for k = 1:length(ep)
        fprintf('k = %i\n',k);
        [logZ(k),S] = xmod10(rbf,sig,log10(ep(k)),inflag,d,prior_target,msk,DI,DE,WW,S.alpha,S.beta);
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
    [LOGZ,SMX] = xmod10(rbf,sig,log10(epmx),inflag,d,prior_target,msk,DI,DE,WW,SMX.alpha,SMX.beta);

    % save the interpolation
    S = SMX;
    S.age = msk+nan;
    S.age(iwet) = exp(S.mod3d(iwet));
    % do a diagnostic figure
    set(0,'CurrentFigure',3);
    hold on
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
    avg_age = wem*S.w;
    tmp = ((wem)*L);
    std_avg_age = sqrt(tmp*tmp'); % standard deviation of global average c-14 age
    avg_age = W'*(S.age(iwet));
    
    fprintf('global average of age %f +/- %f\n',avg_age,std_avg_age);

    mask = msk;
    mask(:,:,3:end) = 0;
    WR = mask(iwet).*grd.dVt(iwet);
    WR = WR/sum(WR);
    res_age = WR'*(S.age(iwet));
    wem = WR'*S.EM;
    %tmp = ((wem)*L)/sqrt(var_d);
    tmp = wem*L;
    std_res = sqrt(tmp*tmp');
    S.res_age = res_age;
    S.std_res = std_res;
    fprintf('average surface reservoir age %f +/- %f\n',res_age,std_res);


    S.avg_age = avg_age;
    S.std_avg_age = std_avg_age;
    S.glodap_age = A3D;
    S.dat = d;
    S.mod = S.IM*S.w;
    S.WW = WW;
    S.iwet = iwet;
    S.fname = fname;

end




function [logZ,S] = xmod10(rbf,sig,logep,inflag,d,prior_target,msk,DM_data,DM_eval,WW,varargin)
    persistent alpha beta
    if (nargin==12)
        alpha = varargin{1};
        beta = varargin{2};
    end

    %% Prepare the evaluation and interpolation matrices               %
    ep = 10^logep;    
    IM = rbf(ep,DM_data); % Interpolation matrix
    EM = rbf(ep,DM_eval); % Evaluation matrix
    
    %%                                                                 %

    
    %% log prior  = log(prob(w|alpha))                                 %
    nw = size(IM,2);
    wp = EM\prior_target; 
    C = EM'*EM;
    % reparameterize w by shifting, rotating, and rescaling
    [Vc,Dc] = svd(C);
    Dc = diag(sqrt(diag(Dc)));
    transform_back = @(dw) Dc*(Vc\dw);
    I = eye(size(Dc,1));
    Ew = @(dw) 0.5*dw'*I*dw;    
    %%                                                                 %

    %% log likelihood  = log(prob(d|w,beta))                           %
    N = length(d);
    sd = d-IM*wp; % shift d so that dw is centered at 0
    IMT = IM*Vc/Dc;
    Ed = @(dw) 0.5*(IMT*dw-sd)'*WW*(IMT*dw-sd);
    B = IMT'*WW*IMT; 
    %%                                                                 %


    %% neg log posterior = log(prob(w|d,alpha,beta))                   %
    A = @(alpha,beta) alpha*I+beta*B;
    M = @(w,alpha,beta) alpha*Ew(w)+beta*Ed(w);
    w_hat = @(alpha,beta) A(alpha,beta)\(beta*IMT'*WW*sd);

    %%                                                                 %

    
    logZ_function = @(alpha, beta) -M( w_hat(alpha,beta), alpha,beta)...
        - 0.5*logdet(A(alpha,beta)) ...
        + 0.5*N*log(beta) + 0.5*(nw-1)*log(alpha);
    

    
    %% Find alpha and beta that maximize the evidence             %
    %  2*alpha * Ew(w_hat(alpha,beta)) = gam(alpha,beta)
    %  2*beta * Ed(w_hat(alpha,beta))  = N-gam(alpha,beta)
    %
    
    if (isempty(alpha)) 
        alpha = 10;
    end
    
    if (isempty(beta)) 
        beta = 1;
    end
    wmp = w_hat(alpha,beta); % compute the most probable weights
    e = 1; tol = 1e-6; outflag = 1;
    it = 0;
    

    % do 15 re-estimation iterations before going to fminsearch 
    % to ensure that fminsearch converges

    if (1)
        while ( (e>tol) & (it<50) )
            it = it+1;
            logP = logZ_function(alpha,beta);
            
            %lama = beta*lamb;
            %gam = sum(lama./(lama+alpha));
            gam = nw - alpha*trace(inv(A(alpha,beta)));
            g = [2*alpha*Ew(wmp)-gam;...
                 2*beta*Ed(wmp)-N+gam];
            e = sqrt(g'*g);
            fprintf(' it = %i ||g|| = %e logP = %e \n',it,e,logP);
            
            % use re-estimation to optimize alpha and beta
            alpha_test = gam/(2*Ew(wmp));
            beta_test = (N-gam)/(2*Ed(wmp));
            if beta_test<0
                beta_test = beta;
            end
 
            logP_test = logZ_function(alpha_test,beta_test);
            if (~isreal(logP_test))
                fprintf('imaginary!!!????\n');
                alpha_test = alpha_test+1e-6;
                beta_test = beta_test;
                logP_test = logZ_function(alpha_test,beta_test);
                if (isreal(logP_test))
                    fprintf('OK now!:-)\n');
                    alpha = alpha_test;
                    beta = beta_test;
                else
                    fprintf('Still not OK. :-( \n');
                    break
                end
            else
                alpha = alpha_test;
                beta = beta_test;
            end
            
            wmp = w_hat(alpha,beta);
                        
            if (it == 1000)
                outflag = 0;
            end
            if ((it>2)&(isinf(Ed(wmp))|isinf(logP)|isnan(logP)))
                outflag = 0;
                fprintf('WTF?\n'); % something went wrong
                keyboard
                break
            end
        end
    end
    
    
    F = @(p) -logZ_function(exp(p(1)),exp(p(2)));
    options = optimset(@fminsearch);
    options = optimset(options,'MaxIter',1000);
    options = optimset(options,'TolX',1e-4);
    options = optimset(options,'TolFun',1e-4);
    options = optimset(options,'MaxFunEvals',1000);
    options = optimset(options,'Display','off');

    % double check that we have the most probable alpha and beta
    [phat,fval,exitflag] = fminsearch(F,[log(alpha);log(beta)],options);
    alpha = exp(phat(1));
    beta = exp(phat(2));
    logP = -fval;
    wmp = w_hat(alpha,beta);
    fprintf('alpha = %f beta = %f logP = %f exitflag = %i\n',alpha,beta,logP,exitflag);

   
    wmp = wp+(Vc/Dc)*wmp;

    % make a diagnostic figure 
    % plot the standardized model against the standardized data
    set(0,'CurrentFigure',3)
    std_d = std(d);
    mean_d = mean(d);
    errorbar((d-mean_d)/std(d),(IM*wmp-mean_d)/std(d),sig,'d');    
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
    S.mod3d = msk+nan;
    S.mod3d(iwet) = EM*wmp;      % total age (NOT the anomaly)
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


