function log_err_euler = ...
    EulerErr(utility_3D_k_diff, utility_3D_klead_diff,tm_2D,opt_k_idx,idx_A_s,idx_z_s)

global nAgrid nzgrid nkgrid beta v_kgrid
% evaluate at the optimal k level
utility_2D_k_diff = nan(nAgrid*nzgrid,nkgrid);
utility_2D_klead_diff = nan(nAgrid*nzgrid,nkgrid);
for state = 1:nAgrid*nzgrid
    for iter_k = 1:nkgrid
        temp = opt_k_idx(state,iter_k);
        % LHS
        utility_2D_klead_diff(state,iter_k) = utility_3D_klead_diff(state,iter_k,temp);
        % RHS
        utility_2D_k_diff(state,iter_k) = utility_3D_k_diff(state,temp,opt_k_idx(state,temp));
    end
end

err_euler = utility_2D_klead_diff+...
    beta*tm_2D*utility_2D_k_diff;
% notice that utility_2D_k_diff is negative.

% compute the logged euler error
log_err_euler = log(abs(err_euler(:,2:nkgrid-1)));
% steady state logged euler error
log_err_euler_s = reshape(log_err_euler,[nAgrid,nzgrid,nkgrid-2]);
log_err_euler_s = log_err_euler_s(idx_A_s,idx_z_s,:);

% plot figures
figure();

sgtitle('Logged Euler Error')
subplot(1,2,1);
plot(v_kgrid(2:nkgrid-1),permute(log_err_euler_s,[3,1,2]),...
    'b','LineWidth',3.0);
title('Steady State');xlabel('$k$')

subplot(1,2,2);
surf(log_err_euler);colorbar;view(0,90);
title('All State');xlabel('$k$'),ylabel('state');

end

