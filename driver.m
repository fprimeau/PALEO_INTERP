% initiaze figure windows
figure(1); figure(3);

%
% load the OCIM grid 
%
OCIM_STUFF = load('myMTMJan2015.mat');

%
% load the diffusion-based distance computd on the OCIM grid 
%
%   J: list of OCIM wet points for which the diffusion-based distance has
%      already been computed
%
%   DE: distance from points in J to all the other OCIM grid points
%
versions = {'V1','V2','V3','V4','V5','A2','A3'};
for v = 7:7
    version = versions{v};
    load(['OCIM_DIFFU_DIST1_',version,'.mat'],'DE'); 
    load(['OCIM_DIFFU_DIST1_',version,'.mat'],'J'); 
    
    slice_name = {'LGM','HS1','BA','YD','HOL','EHOL'};
    % step through the six time slices
    for slice = 1:6
        %
        % 1. read in the xls data spread sheet
        % 2. bin the data onto the OCIM grid
        % 3. find the basin names where the points with data are located
        % 4. write a text file with the data
        %
        
        % modeflag choices:
        % 'unique': first average the data that falls into the same OCIM grid box and use the average
        % 'all': do not average keep all the data and let the interpolation decide how to treat the inconsistent measurements
        % I think the 'all' option is better
        modeflag = 'all';  

        %modeflag = 'unique';

        % choose an radial basic function
        %   rbf_type = 'gaussian';
        %   rbf_type = 'inverse_quadric';
        %   rbf_type = 'quadric';
        %   rbf_type = 'exp_linear';
        %   rbf_type = 'exp_cubic';
        %   rbf_type = 'distance';
        %   rbf_type = 'exp';
        %
        % I tested only the gaussian and the inverse quadric and both seem to work well.
        rbf_type = 'inverse_quadric';
        
        
        % 
        % extract the data from the spread-sheat
        %
        [fname,basin_names,A3D] = core(slice,modeflag);    

        %
        % do the interpolation
        %
        %S{slice} = run9b(fname,A3D,modeflag,rbf_type,version);
        S{slice} = run11(fname,A3D,modeflag,rbf_type,version);
        S{slice}.rbf_type = rbf_type;
    
        figure(1)
        eval(['print -dpng ',slice_name{slice},'_',version,'.png']);
        figure(3)
        eval(['print -dpng ',slice_name{slice},'_',version,'_modVobs.png']);
        
    end
    % save the interpolations 
    eval(['save ',rbf_type,'_',version,'.mat S']);
    
end
    