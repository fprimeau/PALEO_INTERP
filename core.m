function [fname,basin_names,A3D] = core(slice,modeflag)
% if modeflag == 1 average all the data in the same OCIM grid box and assign an uncertainty to the average
% if modeflag == 0 keep the multiple data points that fall in the same OCIM grid box 
    fprintf('Time-slice options:\n')
    switch(slice)
      case 1,
        fprintf('1. LGM (interpolation and offsets)\n')
        xlsname = 'LGM_R_2020.xlsx';
      case 2,
        fprintf('2. HS1 (interpolation and offsets)\n')
        xlsname = 'HS1_R_2020.xlsx';
      case 3,
        fprintf('3. BA (interpolation and offsets)\n')
        xlsname = 'BA_R_2020.xlsx';
      case 4,
        fprintf('4. YD (interpolation and offsets)\n')
        xlsname = 'YD_R_2020.xlsx'; 
      case 5,
        fprintf('5. HOL (scatter of HOL data on glodap)\n')
        xlsname = 'HOL_R_2020.xlsx';
      case 6,
        fprintf('6. EHOL (interpolation and offsets)\n')
        xlsname = 'EHOL_R_2020.xlsx';
      case 7,
        fprintf('7. EHOL (debug case)\n')
        xlsname = 'debug.xlsx';
    end
    
    %
    % Use OCIM grid to do the interpolation
    %
    
    %
    % read in the data from the xls spread sheet
    %
    
    d = xlsread(xlsname);
    lat = d(:,1); % latitutde
    lon = d(:,2); % longitude E
    z = d(:,3);   % depth
    res = d(:,13);
    res_err = d(:,14);
    bp = d(:,10); % bp
    bp_err = d(:,11);
    batm = d(:,15); % b-atm
    batm_err = d(:,16);
    
    
    % remove obviously bad or missing data
    ikp = find(~isnan(lat+lon+z+batm));
    lat = lat(ikp);
    lon = lon(ikp);
    z = z(ikp);
    res = res(ikp);
    res_err = res_err(ikp);
    
    bp = bp(ikp);
    bp_err = bp_err(ikp);
    
    batm = batm(ikp);
    batm_err = batm_err(ikp);
    
    % read in the OCIM grid and masks
    load('transport_Redi_Jan2013.mat','MSKS');
    tmp = load('myMTMJan2015.mat','M3d','grd');
    grd = tmp.grd; % structure with the mesh definitino
    msk = tmp.M3d; % the wet-dry mask (wet == 1, dry == 0)
    x = lon+(lon<0).*360;
    xt = grd.xt;
    yt = grd.yt;
    zt = grd.zt;
    [ny,nx,nz] = size(msk);
    ii = 1:ny;
    jj = 1:nx;
    kk = 1:nz;
    iy = interp1(yt,ii,lat,'nearest');
    ix = interp1(xt,jj,x,'nearest');
    iz = interp1(zt,kk,z,'nearest');
    
    iwet = find(msk);
    J = msk + nan;
    J(iwet) = 1:length(iwet);
    [iix, iiy] = meshgrid(1:nx,1:ny);
    
    fprintf('1st PASS:\n')
    for kk = 1:length(iy)
        if (msk(iy(kk),ix(kk),iz(kk))==1)
            fprintf('   %4i (%3i,%3i,%3i) = (%6.2f,%6.2f,%7.2f)\n',kk,iy(kk),ix(kk),iz(kk),lat(kk),x(kk),z(kk));
        else
            fprintf('???%4i (%3i,%3i,%3i) = (%6.2f,%6.2f,%7.2f) \n',kk,iy(kk),ix(kk),iz(kk),lat(kk),x(kk),z(kk));
        end
    end
    
    topo = sum(msk,3); % water depth in levels
    fname = sprintf('%s.txt',xlsname(1:end-5));
    fid = fopen(fname,'w');
    fprintf(fid,' lat  lon depth  iy   ix iz  iy  ix iz       j   batm  batm_err   res  res_err basin\n');
    fprintf(    ' lat  lon depth  iy   ix iz  iy  ix iz       j   batm  batm_err   res  res_err basin\n');
    for k = 1:length(ix) % loop over all the data
        fprintf(fid,'%4.0f %4.0f %5.0f %3i %4i %2i ',...
                lat(k),x(k),z(k),iy(k),ix(k),iz(k));
        fprintf('%4.0f %4.0f %4.0f %3i %4i %2i ',...
                lat(k),x(k),z(k),iy(k),ix(k),iz(k));
        iflag = 0;
        for i0 = 0:2 % depth shift
            for i1=0:8  % outer shift index loop
                for i2 = 0:8 % inner shift index loop
                    for i3 = 0:8 % inner most shift index loop
                        for i4 = 0:8; % inter most shift index loop
                            test = nb(nb(nb(nb(topo,i1),i2),i3),i4);
                            if (test(iy(k),ix(k))>=iz(k)-i0) % test water column depth
                                                             % shift the vector index
                                jj = nb(nb(nb(nb(J,i1),i2),i3),i4);      
                                j(k) = jj(iy(k),ix(k),iz(k)-i0);
                                jres(k) = jj(iy(k),ix(k),1);
                                % shift the cartesian index
                                xx = nb(nb(nb(nb(iix,i1),i2),i3),i4);  six = xx(iy(k),ix(k));
                                yy = nb(nb(nb(nb(iiy,i1),i2),i3),i4);  siy = yy(iy(k),ix(k));
                                iz(k) = iz(k)-i0;
                                iflag = 1;          
                                break
                            end
                        end 
                        if (iflag==1)
                            break
                        end
                    end
                    if (iflag ==1)
                        break
                    end
                end 
                if (iflag==1)
                    break
                end
            end 
            if (iflag == 1)
                break
            end
        end
        
        if (iflag == 0) 
            j(k) = nan;
            basin = nan;
        else
            if (MSKS.ATL(iwet(j(k))) == 1)
                if (grd.YT3d(iwet(j(k)))>0)
                    basin = 1;
                else
                    basin = 2;
                end
            elseif (MSKS.PAC(iwet(j(k))) == 1)
                if (grd.YT3d(iwet(j(k)))>0)
                    basin = 3;
                else
                    basin = 4;
                end
            elseif (MSKS.IND(iwet(j(k))) == 1)
                basin = 5;
            elseif(MSKS.ARC(iwet(j(k)))==1)
                basin = 6;
            elseif(MSKS.MED(iwet(j(k)))==1)
                basin = 7;
            end
            if (grd.YT3d(iwet(j(k))) < -36)
                basin = 8;
            end
        end
        basin_names{1} = 'NATL';
        basin_names{2} = 'SATL';
        basin_names{3} = 'NPAC';
        basin_names{4} = 'SPAC';
        basin_names{5} = ' IND';
        basin_names{6} = ' ARC';
        basin_names{7} = ' MED';
        basin_names{8} = '  SO';
        fprintf(fid,'%3i %3i %2i  %6i   %4.0f   %7.0f %5.0f  %7.0f %10i\n',...
                siy,six,iz(k),j(k),batm(k),batm_err(k),res(k),res_err(k),basin);
        fprintf('%3i %3i %2i  %6i %4.0f   %7.0f   %5.0f  %7.0f %10i\n',...
                siy,six,iz(k),j(k),batm(k),batm_err(k),res(k),res_err(k),basin);
    end 
    
    fclose(fid);
    
    fprintf('2nd PASS:\n')
    d = importdata(fname);
    ix = d.data(:,8);
    iy = d.data(:,7);
    iz = d.data(:,9);
    for kk = 1:length(iy)
        %  if (~isnan(ix(kk)+iy(kk)+iz(kk)))
        if (msk(iy(kk),ix(kk),iz(kk))==1)
            fprintf('   %4i (%3i,%3i,%3i) = (%6.2f,%6.2f,%7.2f)\n',kk,iy(kk),ix(kk),iz(kk),lat(kk),x(kk),z(kk));
        else
            fprintf('???%4i (%3i,%3i,%3i) = (%6.2f,%6.2f,%7.2f) \n',kk,iy(kk),ix(kk),iz(kk),lat(kk),x(kk),z(kk));
        end
    end 
    

    
    % read in the data file we just created
    [a,a_err,j,x,y,z,ix,iy,iz] = load_core_data(fname);
    dat = importdata(fname); basin = dat.data(:,15); clear dat; % this is for plotting purposes only
    
    % read in the glodap c14 age file
    [glodap_c14_age,A3D] = load_glodap_im(j,msk);
    % create the age anomaly (paleo minus modern)
    anom = a-glodap_c14_age;
    
    
    if (modeflag == 0)
        %
        % Second file with no averaging, all the data are kept 
        % even if there are multiple measurements in the same box
        %
        name = ['all_',fname];
        fid = fopen(name,'w');
        fprintf(fid,' lat   lon depth   iy    ix  iz      j glodap <batm> <batm_err> basin \n');
        for i = 1:length(j)
            fprintf(fid,'%4.0f  %4.0f  %4.0f  %3i  %4i  %2i ',...
                    y(i),x(i),-z(i),iy(i),ix(i),iz(i));
            fprintf(fid, '%6i   %4.0f  %4.0f       %4.0f  %4i\n',...
                    j(i),glodap_c14_age(i),anom(i),a_err(i),basin(i));
        end
        fclose(fid);
    elseif (modeflag == 1)
        %
        % Second file with ONLY ONE AGE PER OCIM GRID BOX
        %
        
        % assemble a dataset of unique locations by performing a weighted average of all the data points
        % within a given grid box.
        uj = unique(j);
        
        d = zeros(length(uj),1);
        err = zeros(length(uj),1);
        n = zeros(length(uj),1);
        bsn = zeros(length(uj),1);
        xx = zeros(length(uj),1);
        yy = zeros(length(uj),1);
        zz = zeros(length(uj),1);
        ixx = zeros(length(uj),1);
        iyy = zeros(length(uj),1);
        izz  = zeros(length(uj),1);
        
        for i = 1:length(uj)
            k = find(j == uj(i));
            fprintf('%i ',k); fprintf('\n');
            
            n(i) = length(k);
            bsn(i) = basin(k(1));
            xx(i) = x(k(1));
            yy(i) = y(k(1));
            zz(i) = z(k(1));
            ixx(i) = ix(k(1));
            iyy(i) = iy(k(1));
            izz(i) = iz(k(1));
            
            % estimate the mean and standard error using a Bayesian estimation method 
            % using a prior with a finite cut off at mu +/- RANGE years
            RANGE = 1500;
            [mu,sig] = fitit(anom(k),a_err(k),RANGE,xx(i),yy(i),zz(i));
            fprintf('%i mu = %f sig = %f\n',i,mu, sig);
            
            d(i) = mu; 
            err(i) = sig;
        end
        name = ['unique_',fname];
        fid = fopen(name,'w');
        fprintf(fid,' lat   lon depth   iy    ix  iz      j glodap <batm> <batm_err> basin \n');
        for i = 1:length(uj)
            fprintf(fid,'%4.0f  %4.0f  %4.0f  %3i  %4i  %2i ',...
                    yy(i),xx(i),-zz(i),iyy(i),ixx(i),izz(i));
            fprintf(fid, '%6i   %4.0f  %4.0f       %4.0f  %4i\n',...
                    uj(i),glodap_c14_age(i),d(i),err(i),bsn(i));
        end
        fclose(fid);
    end
end
function [r1] = nb(M,n)
% NEIGHBOUR creates a matrix containing the n'th neighbour for each
%           element of M. The identification of the neighbour is given
%           by the following nine-point stencil:
%
%               .    .    .
%               6    3    5
%            ^
%            |  .    .    .
%            |  2    0    1
%         i  |    
%            |  .    .    .
%            |  7    4    8
%            +----------->
%                  j
%
%          NOTE (1) the direction of increase of i and j in M(i,j).
%               (2) neighbours outside of M are assumed to be periodic
%
% --input parameters--
%        M is the matrix of elements for which the neighbours will be
%          found.
%
%        n is the index identifying which neighbour to be found
%          according to the stencil given above.
%
% --output parameters--
%
%        r1 is the matrix of neighbours
%
    if     (n == 0)
        r1 = M;
    elseif (n == 1)
        r1 = M(:,[2:end,1],:);
    elseif (n == 2)
        r1 = M(:,[end,1:end-1],:,:);
    elseif (n == 3)
        r1 = M([2:end,1],:,:);
    elseif (n == 4)
        r1 = M([end,1:end-1],:,:);
    elseif (n == 5)
        r1 = M([2:end,1],[2:end,1],:);
    elseif (n == 6)
        r1 = M([2:end,1],[end,1:end-1],:);
    elseif (n == 7)
        r1 = M([end,1:end-1],[end,1:end-1],:);
    elseif (n == 8)
        r1 = M([end,1:end-1],[2:end,1],:);
    end
    
end


function [mu,sig] = fitit(d,err,RANGE,lat,lon,depth)
    % Max standard deviation is RANGE/2
    N = length(d);  % number of observations at a given station
    
    W = @(psig) inv( diag( err.^2 + exp(2*psig) ) ); % the precision matrix
    
    % form the negative log-likelihood
    L = @(mu,psig) -0.5*log( det( W(psig) ) ) + ... 
        0.5 * ( mu - d )' * W(psig) * ( mu - d );

    options = optimoptions('fminunc','SpecifyObjectiveGradient',true);
    
    p0 = [mean(d),log(std(d)+0.1*mean(d))]; % try a 10% relative error as the initial guess
    [popt,Lval,exitflag,output] = fminsearch( @(p) L(p(1),p(2)),p0 );
    fprintf('exitflag = %i\n',exitflag);
    mu = popt(1);
    x = linspace(mu-RANGE,mu+RANGE,200);
    y = logspace(-1,4,200);
    [X,Y] = meshgrid(x,y);
    LL = 0*X;
    for k = 1:length(X(:))
        LL(k) = -L(X(k),log(Y(k)));
    end
    LL = exp(LL-max(LL(:)));
    
    prob  = trapz(y,LL,1);
    Z = trapz(x,prob);
    prob = prob/Z;
    
    set(0,'CurrentFigure',3);
    plot(x,prob,'-k','LineWidth',2); 
    grid on
    xlabel('paleo minus glodap c14 age anomaly (yrs)');
    ylabel('probability density (1/yrs)')
    set(gca,'FontSize',16);
    mu = trapz(x,x.*prob);
    var = trapz(x,(x-mu).^2.*prob);
    sig = sqrt(var);
    text(mu-RANGE,1.1*max(prob),sprintf('mean = %4.0f yrs, std = %4.0f yrs, n = %2i',mu,sig,N),'FontSize',15);
    axis([mu-RANGE,mu+RANGE,0,1.2*max(prob)]);
    title(sprintf('lat = %3.0f lon = %3.0f depth = %4.0f m',lat,lon,abs(depth)))
    drawnow
end


