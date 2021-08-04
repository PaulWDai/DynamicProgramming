% naive version of policy function iteration

clear;
clc;
close;

alpha = .7;
beta = .98;
delta = .1;
v_Agrid = [.9, 1, 1.1];
nAgrid = length(v_Agrid);

v_kgrid = 10:40:8000;
nkgrid = length(v_kgrid);
min_k = min(v_kgrid);
max_k = max(v_kgrid);
step_k = (max_k-min_k)/(nkgrid-1);

tm_2D = [...
0.9,   0.1,    0;...
0.05,  0.9,    0.05;...
0,     0.1,    0.9;...
];

% initial guess
%K2D = 1:1:length(v_kgrid);
%K2D = repmat(K2D,[nAgrid,1]);
K2D = ones(nAgrid,nkgrid);
iter = 0;
err = 10^9;

iter_max = 1000;
err_max = 1;

[m_A,m_k,m_klead]  = ndgrid(v_Agrid,v_kgrid,v_kgrid);
m_c = m_A.* m_k.^alpha + (1-delta).*m_k - m_klead;
m_u_c = log(m_c);

[m_A2,m_k2] = ndgrid(v_Agrid,v_kgrid);

m_r = alpha*m_klead.^(alpha-1);

while iter<iter_max && err>err_max
    iter = iter+1;
    m_kleadlead_idx = nan(nAgrid,nkgrid);
    for state = 1:nAgrid
        for k_iter = 1:nkgrid
            temp = K2D(state,k_iter);
            m_kleadlead_idx(state,k_iter) = K2D(state,temp);
        end
    end
    m_kleadlead = v_kgrid(m_kleadlead_idx);
    m_u2 = m_A2.*m_k2.^(alpha-1)-(1-delta)*m_k2 + m_kleadlead;
    m_u2(m_u2<=0) = -.000001;
    m_u_prime2 = 1./m_u2;
    log_m_u_prime2 = log(m_u_prime2);
    log_m_u_prime2(m_u_prime2<=0) = -999999;
    rhs= beta*tm_2D*(m_kleadlead.*log_m_u_prime2);
    klead = m_A2.*m_k2.^alpha+(1-delta).*m_k2-1./rhs;
    klead_idx = floor((klead-min_k)/step_k+1);
    klead_idx(klead_idx>nkgrid) = nkgrid;
    klead_idx(klead_idx<1) = 1;
    err_temp = abs(klead_idx-K2D);
    err = sum(err_temp,'all');
    err_vec(iter) = err;
    K2D = klead_idx;
end

K2D = v_kgrid(K2D);

row1 = K2D(1,:)-v_kgrid;
row2 = K2D(2,:)-v_kgrid;
row3 = K2D(3,:)-v_kgrid;

% plot the policy function
figure(1);
subplot(1,3,1);
plot(v_kgrid,K2D(1,:),'r','LineWidth',3.0);hold on;
plot(v_kgrid,v_kgrid,'b','LineWidth',3.0);
subplot(1,3,2);
plot(v_kgrid,K2D(2,:),'r','LineWidth',3.0);hold on;
plot(v_kgrid,v_kgrid,'b','LineWidth',3.0);
subplot(1,3,3);
plot(v_kgrid,K2D(3,:),'r','LineWidth',3.0);hold on;
plot(v_kgrid,v_kgrid,'b','LineWidth',3.0);


