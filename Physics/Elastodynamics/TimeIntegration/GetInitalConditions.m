%> @file  GetInitalConditions.m
%> @author Ilario Mazzieri
%> @date 11 May 2023
%> @brief Compute initial displacement and velocity
%>
%==========================================================================
%> @section classGetInitalConditions Class description
%> @brief  GetInitalConditions
%
%> @param Data        Struct with problem's data
%> @param femregion   Finite Element struct (see CreateDOF.m)
%> @param Matrices    Struct with problem's matrices 
%
%> @retval Uold       Vector with initial conditions
%> @retval Solutions  Struct with initial conditions for saving purposes
%>
%==========================================================================

function [Uold, Solutions] = GetInitalConditions(Data,femregion, Matrices)

[ue0] = ComputeModalSolutionWave(Data,femregion);

ue0t   = ue0   * Data.ue_t_ex{1}(0);
due0t   = ue0   * Data.due_t_ex{1}(0);

ue0t   = Matrices.Ela.MPrjP \ue0t;
due0t   = Matrices.Ela.MPrjP \due0t;

Uold = [ue0t; due0t];

Solutions.ue_h  = ue0t;    
Solutions.dot_ue_h  = due0t;  
