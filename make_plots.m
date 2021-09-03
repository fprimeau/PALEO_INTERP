rbf_type = 'inverse_quadric'
if(strcmp(rbf_type,'gaussian'))
    load gaussian.mat
elseif(strcmp(rbf_type,'inverse_quadric'))
    load inverse_quadric.mat
end
   
load('transport_Redi_Jan2013.mat','M3d','MSKS','grd');
PAC = MSKS.PAC;
ATL = MSKS.ATL;


Batm = S{1}.age;
inan = find(isnan(Batm(:)));
Batm(inan) = 0;
pac = PAC;
pac(inan) = 0;
atl = ATL;
atl(inan) = 0;

% make a pacific land mask
pmsk = squeeze(sum(pac.*M3d,2))';
imsk = find(pmsk(:)==0);
pmsk = 0*pmsk;
pmsk(imsk) = 1;
% make an atlantic land mask
amsk = squeeze(sum(atl.*M3d,2))';
imsk = find(amsk(:)==0);
amsk = 0*amsk;
amsk(imsk) = 1;

clevels = [-2000:200:2000];

figure(1)
scr_sz = get(0,'ScreenSize');
set(gcf,'Position',[100,100,floor(scr_sz(3)*0.4),floor(scr_sz(4)*0.8)]);

%
ttl1 = 'PACIFIC E.HOL-YD'; ttl2 = 'ATLANTIC E.HOL-YD';
A = S{6}.age-S{4}.age; A(inan) = 0;
PA = squeeze(sum(A.*pac.*M3d.*grd.dVt,2))./squeeze(sum(pac.*M3d.*grd.dVt,2));
AA = squeeze(sum(A.*atl.*M3d.*grd.dVt,2))./squeeze(sum(atl.*M3d.*grd.dVt,2));
plotit(PA',AA',grd,amsk,pmsk,clevels,ttl1,ttl2,0.6);


%
ttl1 = 'PACIFIC YD-BA'; ttl2 = 'ATLANTIC YD-BA';
B = S{4}.age-S{3}.age; B(inan) = 0;
PB = squeeze(sum(B.*pac.*M3d.*grd.dVt,2))./squeeze(sum(pac.*M3d.*grd.dVt,2));
AB = squeeze(sum(B.*atl.*M3d.*grd.dVt,2))./squeeze(sum(atl.*M3d.*grd.dVt,2));
plotit(PB',AB',grd,amsk,pmsk,clevels,ttl1,ttl2,0.4);


%
ttl1 = 'PACIFIC BA-HS1'; ttl2 = 'ATLANTIC BA-HS1';
C = S{3}.age-S{2}.age; C(inan) = 0;
PC = squeeze(sum(C.*pac.*M3d.*grd.dVt,2))./squeeze(sum(pac.*M3d.*grd.dVt,2));
AC = squeeze(sum(C.*atl.*M3d.*grd.dVt,2))./squeeze(sum(atl.*M3d.*grd.dVt,2));
plotit(PC',AC',grd,amsk,pmsk,clevels,ttl1,ttl2,0.2);


%
ttl1 = 'PACIFIC HS1-LGM'; ttl2 = 'ATLANTIC HS1-LGM';
D = S{2}.age-S{1}.age; D(inan) = 0;
PD = squeeze(sum(D.*pac.*M3d.*grd.dVt,2))./squeeze(sum(pac.*M3d.*grd.dVt,2));
AD = squeeze(sum(D.*atl.*M3d.*grd.dVt,2))./squeeze(sum(atl.*M3d.*grd.dVt,2));
plotit(PD',AD',grd,amsk,pmsk,clevels,ttl1,ttl2,0);

if (strcmp(rbf_type,'gaussian'))
    print -dpng fig_gaussian.png
elseif (strcmp(rbf_type,'inverse_quadric'))
    print -dpng fig_inverse_quadric.png
end

function g = plotit(P,A,grd,amsk,pmsk,clevels,ttl1,ttl2,pos)
    subplot('position',[0.1 pos+0.1 0.39 0.7/4])
    contourf(grd.yt,-grd.zt,inpaint_nans(P),clevels);
    hold on
    zt = [grd.zt(:); grd.zt(end)+grd.dzt(end)/2];
    zz = [zt(sum(~pmsk)+1)-grd.dzt(end)/2; zt(end)+grd.dzt(end)/2; zt(end)+grd.dzt(end)/2];
    yy = [grd.yt(:); grd.yt(end);grd.yt(1)];
    fill(yy,-zz,[0.3 0.3 0.3])
    set(gca,'TickDir','out')
    set(gca,'YTick',[-5000:1000:0]);
    set(gca,'YTickLabel',{'5000','4000','3000','2000','1000','   0'});
    ylabel('depth (m)')

    if (pos~=0) 
        set(gca,'XTick',[]);
    else
    set(gca,'XTick',[-90:30:90])
    end
    title(ttl1);
    axis([-70 70 -5000 0]);
    set(gca,'FontSize',15);
    
    subplot('position',[0.5 pos+0.1 0.39 0.7/4])
    contourf(grd.yt,-grd.zt,inpaint_nans(A),clevels);
    hold on
    zz = [zt(sum(~amsk)+1)-grd.dzt(end)/2; zt(end)+grd.dzt(end)/2; zt(end)+grd.dzt(end)/2];
    yy = [grd.yt(:); grd.yt(end);grd.yt(1)];
    fill(yy,-zz,[0.3 0.3 0.3])
    set(gca,'TickDir','out')
    set(gca,'YTickLabel',[]);
    if (pos~=0) 
        set(gca,'XTick',[]);
    else
    set(gca,'XTick',[-90:30:90])
    end
    title(ttl2);
    axis([-70 70 -5000 0]);    
    set(gca,'FontSize',15);
    
    if(pos==0)
        subplot('position',[0.3 0.05 0.4 0.08/4]);
        contourf(clevels,[0 1],[clevels;clevels],clevels);
        set(gca,'YTick',[]);
        set(gca,'TickDir','out')
        set(gca,'XTick',sort([clevels(1:4:end),0]))
        colormap(bluewhitered);
        caxis([-2000 2000]);
        xlabel('^{14}C-age (years)');
        set(gca,'FontSize',15);
        drawnow
        
        g = gca;
    end
end