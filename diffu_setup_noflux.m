function [STUFF] = diffu_setup_noflux(j,dt,logad,MTM,grd,M3d,ix,iy,iz);
% build a set of diffusion bassis functions
    logai = log(1e4); % scale ad relative to ai=1e3 m^2/s
                      % build the diffusion operator
    iocn = find(M3d(:));
    d0 = @(x) spdiags([x(:)],[0],length(x(:)),length(x(:)));
    iocn = find(M3d(:)); n = length(iocn);
    %
    
    %III = M3d+nan;
    %III(iocn) = 1:n;
    %IIE = III; IIE = IIE(:,[2:end,1],:); 
    %IIW = III; IIW = IIW(:,[end,1:end-1],:);
    %IIN = III; IIN = IIN([2:end,1],:,:);
    %IIS = III; IIS = IIS([end,1:end-1],:,:);
    %IIT = III; IIT = IIT(:,:,[end,1:end-1]);
    %IIB = III; IIB = IIB(:,:,[2:end,1]);
    
    %I0= speye(n);
    %IE = I0(IIE(iocn),:);
    %IW = I0(II
    
    D2y = MTM.BCY*d0(exp(logai))*MTM.TY               + MTM.BCY*d0(exp(logai))*d0(MTM.SyN)*MTM.TYZ;
    D2x = MTM.BCX*d0(exp(logai))*MTM.TX               + MTM.BCX*d0(exp(logai))*d0(MTM.SxE)*MTM.TXZ;
    D2zI = MTM.BCZ*d0(exp(logai))*d0(MTM.SxT)*MTM.TZX + MTM.BCZ*d0(exp(logai))*d0(MTM.SyT)*MTM.TZY + ...
           MTM.BCZ*d0(exp(logai))*d0(MTM.S2)*MTM.TZ;
    D2zD = MTM.BCZ*d0(exp(logad))*MTM.TZ;
    DIFI = MTM.VI*(D2x+D2y+0.1*D2zI);
    DIFD = MTM.VI*D2zD;
    D = DIFI+DIFD;
    
    % use trapezoidal rule for time integration
    % R(n+1)-R(n) = 0.5*dt*D*(R(n+1)+R(n))
    % => [I-0.5*dt*D]*R(n+1) = [I+0.5*dt*D]*R(n)
    B = [speye(n)-0.5*dt*D];
    AA = [speye(n)+0.5*dt*D];
    fprintf('Factoring the big diffusion matrix...');
    FB = mfactor(B);
    fprintf('\n')
    % make the initial condition
    jj = union(j,j); % use only unique centers
    
    %find the grid box indices of the impulses
    for k = 1:length(jj)
        ii = find(jj(k)==j);
        iix(k) = ix(ii(1));
        iiy(k) = iy(ii(1));
        iiz(k) = iz(ii(1));
    end
    
    R = spdiags(1./grd.dVt(iocn),0,length(iocn),length(iocn));
    R = full(R(:,jj));
    % make the projection operator onto the data points
    for ii = 1:length(j)
        k(ii) = find(jj == j(ii));
    end
    H = speye(n);
    H = H(:,k);
    H = H';
    Hj = speye(n);
    Hj = Hj(:,j);
    
    
    STUFF.Hj = Hj;
    STUFF.FB = FB;
    STUFF.AA = AA;
    STUFF.H = H;
    STUFF.R = R;
    STUFF.Ro = R;
    STUFF.jj = jj;
    STUFF.iix = iix;
    STUFF.iiy = iiy;
    STUFF.iiz = iiz;
    STUFF.dt = dt;
    %figure(1); figure(2)
    %tau = STUFF.dt*2;
    %is = find(STUFF.Ro(:,1)==1);
    %R = STUFF.Ro(:,10);
    %for i = 1:1000
    %    R = mfactor(STUFF.FB,STUFF.AA*R+(1/tau)*(STUFF.Ro(:,1)-R));
    %    R = R./max(R(:));
    %    G = M3d+nan;
    %    G(iocn) = R(:,1);
    %    set(0,'CurrentFigure',1);
    %    plot(squeeze(G(STUFF.iiy(1),STUFF.iix(1),:)),-grd.zt,'-o'); 
    %    title(num2str(i)); drawnow
    %    set(0,'CurrentFigure',2);
    %    pcolor(grd.xt,grd.yt,G(:,:,1)); shading flat; colorbar; ...
    %        title(num2str(i)); drawnow
    %end
end

