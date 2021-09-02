function [core_data] = load_unique_core_data(fname)
% [a,a_err,j] = load_unique_core_data('file.txt')
    d = importdata(fname);

    j = d.data(:,7);      % index into 3-d array for model grid
    a = d.data(:,9);      % c-14 age anomaly (yrs)
    a_err = d.data(:,10); % uncertainty (yrs)


    z = -d.data(:,3);  % vertical coordinate
    x = d.data(:,2);   % longitude
    y = d.data(:,1);   % latitude
    iz = d.data(:,6);  % OCIM depth index
    ix = d.data(:,5);  % OCIM longitude index
    iy = d.data(:,4);  % OCIM latitude indes

    core_data = [ix,iy,iz,j,x,y,z,a,a_err];
end
