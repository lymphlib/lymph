%> @file  GetInitalConditionsPoroAcuEla.m
%> @author Ilario Mazzieri
%> @date 5 June 2026
%> @brief Compute initial displacement and velocity
%>
%==========================================================================
%> @section classGetInitalConditionsPoroAcuEla Class description
%> @brief  GetInitalConditionsPoroAcuEla
%
%> @param Data        Struct with problem's data
%> @param mesh       mesh struct (region+neighbor)
%> @param femregion   Finite Element struct (see CreateDOF.m)
%> @param Matrices    Struct with problem's matrices 
%
%> @retval Uold       Vector with initial conditions
%> @retval Solutions  Struct with initial conditions for saving purposes
%>
%==========================================================================

function [Uold, Solutions] = GetInitalConditionsPoroAcuEla(Data, mesh, femregion, Matrices)

%% Computation of modal solution
[Solution_t0] = ComputeModalSolutionPoroAcuEla(Data, mesh, femregion, Data.t0);

%% Modal coefficients projection for the solutions
Solutions.up_h  = Matrices.Poro.MPrjP\(Solution_t0.up * Data.up_t_ex{1}(0));
Solutions.wp_h  = Matrices.Poro.MPrjP\(Solution_t0.wp * Data.wp_t_ex{1}(0));
Solutions.phi_h = Matrices.Acu.MPrjA\(Solution_t0.phi * Data.phi_t_ex{1}(0));
Solutions.ue_h  = Matrices.Ela.MPrjP\(Solution_t0.ue * Data.ue_t_ex{1}(0));

%% Modal coefficients projection for the solution time-derivatives
Solutions.dot_up_h  = Matrices.Poro.MPrjP\(Solution_t0.up * Data.dup_t_ex{1}(0));
Solutions.dot_wp_h  = Matrices.Poro.MPrjP\(Solution_t0.wp * Data.dwp_t_ex{1}(0));
Solutions.dot_phi_h = Matrices.Acu.MPrjA\(Solution_t0.phi * Data.dphi_t_ex{1}(0));
Solutions.dot_ue_h  = Matrices.Ela.MPrjP\(Solution_t0.ue * Data.due_t_ex{1}(0));

%% Solution vector
Uold = [Solutions.up_h;     Solutions.wp_h;     Solutions.phi_h;     Solutions.ue_h; ...
        Solutions.dot_up_h; Solutions.dot_wp_h; Solutions.dot_phi_h; Solutions.dot_ue_h];
