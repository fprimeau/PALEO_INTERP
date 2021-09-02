function [a,A3D] = load_glodap_im(j,M3d)
    iocn = find(M3d);
    tauC14 = 5730/log(2);
    load MYOBSdata.mat
    A3D = M3d+nan;
    A3D(iocn) = -log(MYOBS.c14star)*tauC14;
    for k = 1:24
        A3D(:,:,k) = inpaint_nans(A3D(:,:,k));
    end
    inan = find(M3d(:)==0);
    A3D(inan) = nan;
    a = A3D(iocn(j));
end
