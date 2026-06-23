%> @file   IndicatorPreallocationFKPP.m
%> @author Mattia Corti, Caterina Leimer Saglio
%> @date   14 April 2026
%> @brief Preallocation of the indicator structure for the FKPP equation.
%>
%==========================================================================
%> @section classIndicatorPreallocationFKPP Class description
%==========================================================================
%> @brief           Preallocation of the indicator structure for FKPP equation.
%>
%> @param GenMatrices Struct containing the preallocated vector CellVector
%>
%> @retval Indicator  Indicator struct containing components stored using
%> cells associated with the mesh elements
%>                   
%==========================================================================

function [Indicator] = IndicatorPreallocationFKPP(GenMatrices)

    Indicator.Volume.tau_E = GenMatrices.CellVector;   % Residual indicator
    
    Indicator.Faces.tau_J  = GenMatrices.CellVector;   % Jump indicator
    Indicator.Faces.tau_N  = GenMatrices.CellVector;   % Gradient normal component jump indicator
    Indicator.Faces.tau_T  = GenMatrices.CellVector;   % Gradient tangent component jump indicator
end
