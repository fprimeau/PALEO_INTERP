function [a,a_err,j,x,y,z,ix,iy,iz] = load_core_data(fname)
% [a,a_err,j,x,y,z,ix,iy,iz] = load_core_data('file.txt')
    d = importdata(fname);
    j = d.data(:,10);  % index into 3-d array for model grid
    batm = d.data(:,11);
    batm_err = d.data(:,12);
    res = d.data(:,13);
    res_err = d.data(:,14);
    a = batm;                   % c-14 age
                                
    % uncertainty propagation, if using b-p
    % a_err = sqrt(bp_err.^2+res_err.^2);
    a_err = batm_err;
    % to ignore the measurement error uncomment the next line
    % a_err = 1+0*sqrt(bp_err.^2+res_err.^2);
    
    z = -d.data(:,3);  % vertical coordinate
    x = d.data(:,2);   % longitude
    y = d.data(:,1);   % latitude
    iz = d.data(:,9);  % OCIM depth index
    ix = d.data(:,8);  % OCIM longitude index
    iy = d.data(:,7);  % OCIM latitude indes
    
end
