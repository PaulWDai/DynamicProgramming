function [c1_4D_final,c2_4D_final,l_4D_final,utility_4D_final,exitflag_4D] = ...
        instantaneous_utility(v_Agrid,v_kgrid,v_zgrid,nAgrid,nzgrid,nkgrid)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Global:
%   grid number of state variables: nAgrid nzgrid nkgrid
%   setting: options (shut down the fsolve display)
%   parameters: alpha psi

%   Input:
%   v_Agrid, v_kgrid, v_zgrid: grid of state variables

%   Output:
%   solution of instantaneous utility maximixation given k and klead
%   dimension: (A,z,k,klead)
%   c1_4D_final,c2_4D_final,l_4D_final,utility_4D_final
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
global options alpha psi

guess_temp = ones(1,4);

% initialization
c1_4D = nan(nAgrid, nzgrid, nkgrid, nkgrid);
c2_4D = nan(nAgrid, nzgrid, nkgrid, nkgrid);
l_4D = nan(nAgrid, nzgrid, nkgrid, nkgrid);
exitflag_4D = nan(nAgrid, nzgrid, nkgrid, nkgrid);


for iter_A = 1:nAgrid
    for iter_z = 1:nzgrid
        for iter_k = 1:nkgrid
            for iter_klead = 1:nkgrid
                %{
                global ktemp kleadtemp Atemp ztemp
                Atemp = v_Agrid(iter_A);
                ztemp = v_zgrid(iter_z);
                ktemp = v_kgrid(1,iter_k);
                kleadtemp = v_kgrid(1,iter_klead);
                %}
                
                % turn the display off to get rid of Matlab report
                % further check of success solving or not, see exitflag.
                known_value_temp = [v_Agrid(iter_A),v_zgrid(iter_z),v_kgrid(1,iter_k),v_kgrid(1,iter_klead)];
                sol_utility = @(x) utility(x,known_value_temp);
                [solution,fval,exitflag,output] = ...
                    fsolve(sol_utility, guess_temp,options);
                
                c1_4D(iter_A,iter_z,iter_k,iter_klead) = solution(1);
                c2_4D(iter_A,iter_z,iter_k,iter_klead) = solution(2);
                l_4D(iter_A,iter_z,iter_k,iter_klead) = solution(3);
                exitflag_4D(iter_A,iter_z,iter_k,iter_klead) = exitflag;
                
                guess_temp = solution;
                % update the new guess, such that the guessed value is
                % close to the next solution
                % surprisingly, I found it shorten the time a half.
            end
        end
    end
end
%{
c1_4D_final = exp(c1_4D);
c2_4D_final = exp(c2_4D);
l_4D_final = exp(l_4D);
%}
c1_4D_final = c1_4D;
c2_4D_final = c2_4D;
l_4D_final = l_4D;

utility_4D_final = ...
    c1_4D_final.^(alpha).*c2_4D_final.^(1-alpha)-(1/psi).*l_4D_final.^psi;


% report whether all the nonlinear equation is solvable
disp('------------------------------------------------------')
if min(exitflag_4D,[],'all')>=0
    disp('All nonlinear equations are successfully solved')
else
    disp('There exists some nonlinear equation are not solved')
end
disp('------------------------------------------------------')
disp('Time: Calculation of Instantaneous Utility')
toc
disp('------------------------------------------------------')
end

