function NoOutput = ImpulseResponse(time_window,idx_A_s,idx_z_s,k_s,opt_k_idx_3D,v_kgrid,c1_vec,c2_vec,l_vec)

%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

global nkgrid nAgrid nzgrid

k_idx_path_IRF_A = nan(1,time_window);
k_idx_path_IRF_z = nan(1,time_window);
idx_k_s = floor((k_s-min(v_kgrid))/(max(v_kgrid)-min(v_kgrid))*(nkgrid-1)+1);
k_idx_path_IRF_A(1) = idx_k_s;
k_idx_path_IRF_z(1) = idx_k_s;

A_path_idx_base = ones(1,time_window)*idx_A_s;
z_path_idx_base = ones(1,time_window)*idx_z_s;
A_path_idx = A_path_idx_base;
z_path_idx = z_path_idx_base;
A_path_idx(2) = nAgrid;
z_path_idx(2) = nzgrid;

for i = 2:time_window
    k_idx_path_IRF_A(i) = opt_k_idx_3D(A_path_idx(i),z_path_idx_base(i),k_idx_path_IRF_A(i-1));
    k_idx_path_IRF_z(i) = opt_k_idx_3D(A_path_idx_base(i),z_path_idx(i),k_idx_path_IRF_z(i-1));
end

% make the vector of variables with same length
time_window = time_window-1;
klead_idx_path_IRF_A = k_idx_path_IRF_A(1:time_window);
klead_idx_path_IRF_z = k_idx_path_IRF_z(1:time_window);
A_path_idx = A_path_idx(1:time_window);
A_path_idx_base = A_path_idx_base(1:time_window);
z_path_idx = z_path_idx(1:time_window);
z_path_idx_base = z_path_idx_base(1:time_window);
k_idx_path_IRF_A = k_idx_path_IRF_A(1:time_window);
k_idx_path_IRF_z = k_idx_path_IRF_z(1:time_window);

% change the 4D entry location into a 1D location
idx_IRF_A = sub2ind([nAgrid nzgrid nkgrid nkgrid],A_path_idx,z_path_idx_base,k_idx_path_IRF_A,klead_idx_path_IRF_A);
idx_IRF_z = sub2ind([nAgrid nzgrid nkgrid nkgrid],A_path_idx_base,z_path_idx,klead_idx_path_IRF_A);

% convert into percentage for comparison
k_path_IRF_A = v_kgrid(k_idx_path_IRF_A);
k_path_IRF_A = (k_path_IRF_A - k_path_IRF_A(1))/k_path_IRF_A(1);
k_path_IRF_z = v_kgrid(k_idx_path_IRF_z);
k_path_IRF_z = (k_path_IRF_z - k_path_IRF_z(1))/k_path_IRF_z(1);
c1_path_IRF_A = c1_vec(idx_IRF_A);
c1_path_IRF_A = (c1_path_IRF_A-c1_path_IRF_A(1))/c1_path_IRF_A(1);
c1_path_IRF_z = c1_vec(idx_IRF_z);
c1_path_IRF_z = (c1_path_IRF_z-c1_path_IRF_z(1))/c1_path_IRF_z(1);
c2_path_IRF_A = c2_vec(idx_IRF_A);
c2_path_IRF_A = (c2_path_IRF_A-c2_path_IRF_A(1))/c2_path_IRF_A(1);
c2_path_IRF_z = c2_vec(idx_IRF_z);
c2_path_IRF_z = (c2_path_IRF_z-c2_path_IRF_z(1))/c2_path_IRF_z(1);
l_path_IRF_A = l_vec(idx_IRF_A);
l_path_IRF_A = (l_path_IRF_A-l_path_IRF_A(1))/l_path_IRF_A(1);
l_path_IRF_z = l_vec(idx_IRF_z);
l_path_IRF_z = (l_path_IRF_z-l_path_IRF_z(1))/l_path_IRF_z(1);

% Plot the impulse response function
figure(3);
sgtitle('Impulse Response Function');

subplot(4,2,1);
plot(1:time_window,k_path_IRF_A,'r','LineWidth',3.0);
title('saving')

subplot(4,2,2);
plot(1:time_window,k_path_IRF_z,'b','LineWidth',3.0);
title('saving')

subplot(4,2,3);
plot(1:time_window,c1_path_IRF_A,'r','LineWidth',3.0);
title('consumption good 1')
subplot(4,2,4);
plot(1:time_window,c1_path_IRF_z,'b','LineWidth',3.0);
title('consumption good 1')

subplot(4,2,5);
plot(1:time_window,c2_path_IRF_A,'r','LineWidth',3.0);
title('consumption good 2')
subplot(4,2,6);
plot(1:time_window,c2_path_IRF_z,'b','LineWidth',3.0);
title('consumption good 2')

subplot(4,2,7);
plot(1:time_window,l_path_IRF_A,'r','LineWidth',3.0);
xlabel('Time')
title('labor supply')
subplot(4,2,8);
plot(1:time_window,l_path_IRF_z,'b','LineWidth',3.0);
xlabel('Time')
title('labor supply')
end

