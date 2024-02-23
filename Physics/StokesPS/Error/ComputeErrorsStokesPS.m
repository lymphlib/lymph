%> @file  ComputeErrorsStokesPS.m
%> @author Ilario Mazzieri
%> @date 16 February 2024
%> @brief  Compute errors for the unsteady Stokes problem
%>
%==========================================================================
%> @section classComputeErrorsStokesPS Class description
%> @brief  Compute errors for unsteady Stokes problem
%
%> @param Data       Struct with problem's data
%> @param femregion  Finite Element struct (see CreateDOF.m)
%> @param Matrices   Struct with problem's matrices 
%> @param Solutions  Struct with solution vectors 
%
%> @retval Error     Struct with computed errors
%>
%==========================================================================

function [Error] = ComputeErrorsStokesPS(Data, femregion, Matrices, Solutions)

%% Compute modal coefficient of the exact solution
[sigma0] = ComputeModalSolutionStokesPS(Data,femregion);

nf = femregion.ndof_f;

dt = Data.dt;
nts = floor(Data.T/dt);
counter = 1;

for t = dt : dt : nts*dt
    sigma_ex{counter}  = sigma0  * Data.sigma_t_ex{1}(t);
    %dot_sigma_ex  = sigma0  * Data.sigma_dt_ex{1}(t);
    
    tr_sigma_ex = [0.5*sigma_ex{counter}(1:nf)+0.5*sigma_ex{counter}(3*nf+1:end); ...
        0*sigma_ex{counter}(1:nf);
        0*sigma_ex{counter}(1:nf); ...
        0.5*sigma_ex{counter}(1:nf)+0.5*sigma_ex{counter}(3*nf+1:end)];
    dev_sigma_ex{counter} = sigma_ex{counter} - tr_sigma_ex;
    counter = counter + 1;
end

%% Switch from nodal to modal representation of the solution coefficients
error_dG = 0;
error_L2_dev_sigma = 0;
error_L2_sigma = 0;
error_Linf_sigma = 0;

for i = 1 : size(sigma_ex,2)

    sigma_ex{i}         = Matrices.Fluid.MPrjP\sigma_ex{i};
    dev_sigma_ex{i}     = Matrices.Fluid.MPrjP\dev_sigma_ex{i};
    % Error vectors
    error_s             = sigma_ex{i} - Solutions{i}.sigma_h;
    error_Linf_sigma    = max(error_L2_sigma, max(abs(sigma_ex{i} - Solutions{i}.sigma_h)));
    error_sigma         = error_s' * Matrices.Fluid.DGf * error_s;
    error_dG            = error_dG + dt*error_sigma;

    error_L2_sigma = max(error_L2_sigma, error_s' * Matrices.Fluid.MPrjP * error_s);
    error_dev_sigma = dev_sigma_ex{i} - Solutions{i}.dev_sigma_h;
    error_L2_dev_sigma = max(error_L2_dev_sigma,...
                             error_dev_sigma'*Matrices.Fluid.M_F*error_dev_sigma);
end

% L2 error
error_L2_vel = error_L2_dev_sigma;
error_L2     = error_L2_sigma;

% Energy error
error_Energy = sqrt(error_L2_vel + error_dG);
error_L2_v   = sqrt(error_L2_vel);
error_L2_d   = sqrt(error_L2);



%% Outputs

Error.nel = femregion.nel;
Error.h   = Data.h;
Error.p   = Data.degree;

Error.error_L2_v        = error_L2_v;
Error.error_L2_d        = error_L2_d;
Error.error_Energy      = error_Energy;
Error.error_L2_vel      = error_L2_vel;
Error.error_L2_dev_sigma = error_L2_dev_sigma;
Error.error_dG          = error_dG;
Error.error_L2_sigma    = error_L2_sigma;
Error.error_inf_sigma   = error_Linf_sigma;


