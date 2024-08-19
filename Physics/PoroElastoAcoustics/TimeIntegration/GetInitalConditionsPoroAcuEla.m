%> @file  GetInitalConditionsPoroAcuEla.m
%> @author Ilario Mazzieri
%> @date 26 June 2024
%> @brief Compute initial displacement and velocity
%>
%==========================================================================
%> @section classGetInitalConditionsPoroAcuEla Class description
%> @brief  GetInitalConditionsPoroAcuEla
%
%> @param Data        Struct with problem's data
%> @param femregion   Finite Element struct (see CreateDOF.m)
%> @param Matrices    Struct with problem's matrices 
%
%> @retval Uold       Vector with initial conditions
%> @retval Solutions  Struct with initial conditions for saving purposes
%>
%==========================================================================

function [Uold, Solutions] = GetInitalConditionsPoroAcuEla(Data,femregion, Matrices)


[up0,wp0,phi0,ue0] = ComputeModalSolutionPoroAcuEla(Data,femregion);

up0t   = up0   * Data.up_t_ex{1}(0);
wp0t   = wp0   * Data.wp_t_ex{1}(0);
phi0t  = phi0  * Data.phi_t_ex{1}(0);
ue0t   = ue0   * Data.ue_t_ex{1}(0);

dup0t   = up0   * Data.dup_t_ex{1}(0);
dwp0t   = wp0   * Data.dwp_t_ex{1}(0);
dphi0t  = phi0  * Data.dphi_t_ex{1}(0);
due0t   = ue0   * Data.due_t_ex{1}(0);

up0t   = Matrices.Poro.MPrjP\up0t;
wp0t   = Matrices.Poro.MPrjP\wp0t;
phi0t  = Matrices.Acu.MPrjA  \phi0t;
ue0t   = Matrices.Ela.MPrjP\ue0t;

dup0t   = Matrices.Poro.MPrjP\dup0t;
dwp0t   = Matrices.Poro.MPrjP\dwp0t;
dphi0t  = Matrices.Acu.MPrjA  \dphi0t;
due0t   = Matrices.Ela.MPrjP\due0t;

Uold = [up0t; wp0t; phi0t; ue0t; ...
    dup0t; dwp0t; dphi0t; due0t];


Solutions.up_h   = up0t;
Solutions.wp_h   = wp0t;
Solutions.phi_h  = phi0t;
Solutions.ue_h   = ue0t;
Solutions.dot_up_h  = dup0t;
Solutions.dot_wp_h  = dwp0t;
Solutions.dot_phi_h = dphi0t;
Solutions.dot_ue_h  = due0t;
