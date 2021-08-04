function [utility_3D_k_diff, utility_3D_klead_diff] = utility_diff(utility_3D)

global nkgrid kstep 
utility_3D_k_forward = nan(size(utility_3D));
utility_3D_k_backward = nan(size(utility_3D));

utility_3D_k_forward(:,1:nkgrid-1,:) = utility_3D(:,2:nkgrid,:); 
utility_3D_k_forward(:,nkgrid,:) = utility_3D_k_forward(:,nkgrid-1,:);

utility_3D_k_backward(:,2:nkgrid,:) = utility_3D(:,1:nkgrid-1,:);
utility_3D_k_backward(:,1,:) = utility_3D(:,2,:);

utility_3D_klead_forward = nan(size(utility_3D));
utility_3D_klead_backward = nan(size(utility_3D));

utility_3D_klead_forward(:,:,1:nkgrid-1) = utility_3D(:,:,2:nkgrid); 
utility_3D_klead_forward(:,:,nkgrid) = utility_3D(:,:,nkgrid-1); 

utility_3D_klead_backward(:,:,2:nkgrid) = utility_3D(:,:,1:nkgrid-1); 
utility_3D_klead_backward(:,:,1) = utility_3D(:,:,2); 

% differentiation wrt k
utility_3D_k_diff = ...
    (utility_3D_k_forward - utility_3D_k_backward)/(2*kstep); 

% differentiation wrt klead
utility_3D_klead_diff = ...
    (utility_3D_klead_forward - utility_3D_klead_backward)/(2*kstep);

end

