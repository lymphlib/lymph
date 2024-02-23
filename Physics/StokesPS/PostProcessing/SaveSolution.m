%> @file  SaveSolution.m
%> @author Ilario Mazzieri
%> @date 16 Febraury 2024
%> @brief  Save solution vector in a struct
%>
%==========================================================================
%> @section classSaveSolutionPS Class description
%> @brief  Save solution vector in a struct
%
%> @param Uh     Solution vector
%> @param femregion  Finite Element struct (see CreateDOF.m)
%
%> @retval Solutions  Struct with solution vector
%>
%==========================================================================

function  [Solutions]  = SaveSolution(Uh,femregion)

nf = femregion.ndof_f;

% fluid
sigma_h   =  Uh;
sigma_h_11 = Uh(     1 :    nf);
sigma_h_12 = Uh(  nf+1 :  2*nf);
sigma_h_21 = Uh(2*nf+1 :  3*nf);
sigma_h_22 = Uh(3*nf+1 :  4*nf);

%    sigma_ho_11 = Uold(     1 :    nf);
%    sigma_ho_12 = Uold(  nf+1 :  2*nf);
%    sigma_ho_21 = Uold(2*nf+1 :  3*nf);
%    sigma_ho_22 = Uold(3*nf+1 :  4*nf);

dev_sigma   =  [0.5*sigma_h_11-0.5*sigma_h_22;  sigma_h_12; sigma_h_21; 0.5*sigma_h_22-0.5*sigma_h_11];
%    dev_sigma_o =  [0.5*sigma_ho_11-0.5*sigma_ho_22; sigma_ho_12; sigma_ho_21; 0.5*sigma_ho_22-0.5*sigma_ho_11];

%    dot_dev_sigma_h = (dev_sigma - dev_sigma_o)/dt;

Solutions.sigma_h   = sigma_h;
Solutions.dev_sigma_h = dev_sigma;

% Solutions.dot_dev_sigma_h = dot_dev_sigma_h;

