function output =  SteadyState(x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Global:
%   z0: steady state level of z (logged TFP of good 1)
%   A0: steady state level of z (TFP of good 2)

%   Input:
%   [logged c1, k] % to ensure a positive consumption


%   Output:
%   If solved succesfully, a vector of zeros.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global z0 A0 beta alpha psi delta eta 

c1 = exp(x(1));
k = x(2);

% steady state r
r = 1/beta - (1-delta);

% good 1 production FOC
w = (1-eta)*(exp(z0)*((eta/r)^eta))^(1/(1-eta));

% % good 2 prodcution FOC
p2 = w/A0;

% steady state capital
klead = k;

%  eq. combining FOC l, c1 and c2
l = (alpha*w*((1-alpha)/(alpha*p2))^(1-alpha))^(1/(psi-1));

% eq. combining FOC c1 and c2
c2 = (1-alpha)*c1/(alpha*p2);

% good 2 market clearing
l2 = c2/A0;

% good 1 production relative price
l1 =  ((r/w)*((1-eta)/eta))*k;

% budget constraint
output(1) = (c1 + p2 * c2 + klead - (1-delta) * k ) - (w*l+r*k);

% labor market clearing
output(2) = l - (l1+l2);
end
