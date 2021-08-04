function tm_2D_final = tran_mat(tm_A,tm_z)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
global nAgrid nzgrid
tm_A_4D = repmat(tm_A,1,1,nzgrid,nzgrid);

tm_z_4D = repmat(tm_z,1,1,nAgrid,nAgrid);
tm_z_4D = permute(tm_z_4D,[3,4,1,2]);

tm_4D = permute(tm_A_4D.*tm_z_4D,[1,3,2,4]);
tm_2D = reshape(tm_4D,[nAgrid*nzgrid,nAgrid*nzgrid]);

tm_2D_final = sparse(tm_2D);
end

