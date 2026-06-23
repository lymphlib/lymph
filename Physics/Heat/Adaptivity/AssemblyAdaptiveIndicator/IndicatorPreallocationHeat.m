%> @file   IndicatorPreallocationHeat.m
%> @author Mattia Corti, Caterina Leimer Saglio
%> @date   6 February 2026
%> @brief Preallocation of the indicator structure for the heat equation.
%>
%==========================================================================
%> @section classIndicatorPreallocationHeat Class description
%==========================================================================
%> @brief           Preallocation of the indicator structure for heat equation.
%>
%> @param GenMatrices Struct containing the preallocated vector CellVector
%>
%> @retval Indicator  Indicator struct containing components stored using
%> cells associated with the mesh elements
%>                   
%==========================================================================

function [Indicator] = IndicatorPreallocationHeat(GenMatrices)

    Indicator.Volume.tau_E = GenMatrices.CellVector;   % Residual indicator
    
    Indicator.Faces.tau_J  = GenMatrices.CellVector;   % Jump indicator
    Indicator.Faces.tau_N  = GenMatrices.CellVector;   % Gradient normal component jump indicator
    Indicator.Faces.tau_T  = GenMatrices.CellVector;   % Gradient tangent component jump indicator
end
