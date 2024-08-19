%> @file  GetInitalConditions.m
%> @author Ilario Mazzieri, Stefano Bonetti
%> @date 24 July 2024
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
%>
%==========================================================================

function [Uold] = GetInitalConditions(Data,femregion, Matrices)

[ue0, ve0] = ComputeModalSolutionWave(Data,femregion,Data.t0);

ue0 = Matrices.Ela.MPrjP \ ue0;
ve0 = Matrices.Ela.MPrjP \ ve0;

Uold = [ue0; ve0];