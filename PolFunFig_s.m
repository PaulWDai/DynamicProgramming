function NoOutput =...
    PolFunFig_s(c1_vec,c2_vec,l_vec,opt_k3D,opt_k_idx_3D,v_kgrid,idx_A_s,idx_z_s)

global nAgrid nzgrid nkgrid 

% policy function at steady state
v_PlcFun_ss = permute(opt_k3D(idx_A_s,idx_z_s,:),[3,1,2]);


opt_k_idx_s = permute(opt_k_idx_3D(idx_A_s,idx_z_s,:),[3,1,2]);
ind_s = sub2ind([nAgrid nzgrid nkgrid nkgrid],...
    idx_A_s*ones(1,nkgrid),...
    idx_z_s*ones(1,nkgrid),...
    [1:nkgrid],...
    opt_k_idx_s');

opt_c1_s = c1_vec(ind_s);
opt_c2_s = c2_vec(ind_s);
opt_l_s = l_vec(ind_s);


% policy function of k fixed A0 and z0
figure(2);
subplot(2,2,1);
plot(v_kgrid,v_PlcFun_ss,'LineWidth',3.0);hold on;
plot(v_kgrid,v_kgrid,'r-','LineWidth',2.0)
xlabel('$k$');ylabel('$k''$');title('saving')

subplot(2,2,2);
plot(v_kgrid,opt_c1_s,'LineWidth',3.0);
xlabel('$k$');ylabel('$c_1$');title('consumption good 1')

subplot(2,2,3);
plot(v_kgrid,opt_c2_s,'LineWidth',3.0);
xlabel('$k$');ylabel('$c_2$');title('consumption good 2')

subplot(2,2,4);
plot(v_kgrid,opt_l_s,'LineWidth',3.0);
xlabel('$k$');ylabel('$l$');title('labor supply')

end

