%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Project 1: Numerical Methods
% Author: Weifeng Dai, weifeng.dai@outlook.com or daiweifeng@pku.edu.cn
% Student ID: 1700018613
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear;
close;

%% Section 1: Parameters
global beta alpha psi delta eta ...
    nAgrid nzgrid nkgrid v_kgrid kstep ...
    options iter_max iter_err 

beta = 0.96;    % discount rate
alpha = 0.5;    % cobb douglas in consumption bundle
psi = 2;        % 1 + 1/ frisch
delta = 0.1;    % depreciation
eta = .33;      % cobb douglas: good 1 capital

% states of z and A
v_zgrid = [-.0673, -.0336, 0, .0336, .0673];
v_Agrid = [.9, 1, 1.1];

% transitional matrix

% good 1 logged TFP z
tm_z =  [...
0.9727, 0.0273, 0,      0,      0;...
0.0041, 0.9806, 0.0153, 0,      0;...
0,      0.0082, 0.9836, 0.0082, 0;...
0,      0,      0.0153, 0.9806, 0.0041;...
0,      0,      0,      0.0273, 0.9727;...
];

% good 2 TFP A
tm_A = [...
0.9,   0.1,    0;...
0.05,  0.9,    0.05;...
0,     0.1,    0.9;...
];

% setting

% optimization setting
options = optimset('Algorithm', 'levenberg-marquardt','Display','off');
%options = optimset('Algorithm', 'levenberg-marquardt');

% graph setting
set(0,'defaultTextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

% stopping rule for iteration setting
iter_err = 10^(-6);
iter_max = 400;

% grid setting
nkgrid = 50;
%% Section 2: Steady state

global z0 A0

z0 = 0;    % steady state
A0 = 1;

guess_s = [1,1];
SteadyStateSol = fsolve(@SteadyState, guess_s,options);
c1_s = exp(SteadyStateSol(1)); % to constrain c is greater than 0
k_s = SteadyStateSol(2);

%% Section 3: Value function iteration

%% 3.1 Value function iteration with a fixed grid

nAgrid = length(v_Agrid);
nzgrid = length(v_zgrid);

% construct the capital grid with equal step
kmin = 0.7*k_s;
kmax = 1.3*k_s;
v_kgrid = linspace(kmin,kmax,nkgrid);
kstep = (kmax-kmin)/(nkgrid-1);

% instantaneous utility

[c1_4D,c2_4D,l_4D, utility_4D,exitflag_4D] = ...
    instantaneous_utility(v_Agrid,v_kgrid,v_zgrid,nAgrid,nzgrid,nkgrid);

utility_4D(c1_4D<=0|c2_4D<=0|l_4D<=0) = -999;
utility_3D = reshape(utility_4D,[nAgrid*nzgrid,nkgrid,nkgrid]);

% transitional matrix
tm_2D = tran_mat(tm_A,tm_z);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Value Function Iteration with Fixed Grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[V2D,opt_k_idx] = ...
    VFI_fixed(utility_3D,tm_2D);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Value Function at Steady State
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% k for different z
idx_A_s = find(v_Agrid==A0);
idx_z_s = find(v_zgrid==z0);

% ValFunFig_s(V2D,v_kgrid,idx_A_s,idx_z_s);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Policy Functions at Steady State
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

opt_k3D = reshape(v_kgrid(opt_k_idx),[nAgrid,nzgrid,nkgrid]);
opt_k_idx_3D = reshape(opt_k_idx,[nAgrid,nzgrid,nkgrid]);

c1_vec = reshape(c1_4D,[nAgrid*nzgrid*nkgrid*nkgrid,1]);
c2_vec = reshape(c2_4D,[nAgrid*nzgrid*nkgrid*nkgrid,1]);
l_vec = reshape(l_4D,[nAgrid*nzgrid*nkgrid*nkgrid,1]);

% plot
% PolFunFig_s(c1_vec,c2_vec,l_vec,opt_k3D,opt_k_idx_3D,v_kgrid,idx_A_s,idx_z_s)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Euler Error (Unfinished)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[utility_3D_k_diff, utility_3D_klead_diff] = ...
    utility_diff(utility_3D);

log_err_euler = ...
    EulerErr(utility_3D_k_diff, utility_3D_klead_diff,tm_2D,opt_k_idx,idx_A_s,idx_z_s);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Impulse Response function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time_window = 20; % the time window for observing response
% calculate the sequence and then plot
ImpulseResponse(time_window,idx_A_s,idx_z_s,k_s,opt_k_idx_3D,v_kgrid,c1_vec,c2_vec,l_vec)

%% 3.2 Value function iteration with an endogenous grid
iter = 0;
err = 10^9;
K2D_idx = repmat(1:nkgrid,[nAgrid*nzgrid,1]);
while iter<iter_max && err>iter_err
    utility_diff_2D = nan(nAgrid*nzgrid,nkgrid);
    for state = 1:nAgrid*nzgrid
        for iter_k = 1:nkgrid
            opt_k_idx_temp = opt_k_idx(state,iter_k);
            opt_opt_k_idx_temp = opt_k_idx(state,opt_k_idx_temp);
            utility_diff_2D(state,iter_k) = ...
                utility_3D_klead_diff(state,opt_k_idx_temp,opt_opt_k_idx_temp);
        end
    end
    rhs = repmat(tm_2D*utility_diff_2D,[1,1,nkgrid]);
    temp = beta*permute(rhs,[1,3,2]) + utility_3D_k_diff;
    [euler_err,opt_k_idx] = min(temp,[],3);
    err = max(euler_err,[],'all')
    K2D_idx = opt_k_idx;
end


%% 3.3 Accelerator

[V2D_accelerator,opt_k_idx_accelerator] = ...
    VFI_accelerator(utility_3D,tm_2D);

%% 3.4 Multigrid

% initialization of multigrid
step1 = 20;     % the coarsest grid number
step2 = 5;      % the second-coarsest grid number

iter_threshold1 = 30;   % times of iteration in the coarsest grid
iter_threshold2 = 50;   % times of iteration in the second-coarsest grid

[V2D_multigrid,opt_k_idx_multigrid] = ...
    VFI_multigrid(step1,step2,iter_threshold1,iter_threshold2,utility_4D,utility_3D,tm_2D);

%% Section 4: Projection

%% Section 5: Perturbation




