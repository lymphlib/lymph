%> @file   FinalIndicatorHeat.m
%> @author Mattia Corti, Caterina Leimer Saglio
%> @date   6 February 2026
%> @brief  Final indicator for heat equation.
%>
%==========================================================================
%> @section classFinalIndicatorHeat Class description
%==========================================================================
%> @brief           Final indicator for heat equation.
%>
%> @param Indicator_loc  Indicator struct containing the local indicators
%>
%> @retval Indicator     Indicator struct containing the different indicators
%>                   
%==========================================================================

function [Indicator] = FinalIndicatorHeat(Indicator_loc)

    Indicator.Adaptivity = Indicator_loc;
    
    Indicator.tau_E = sqrt(Indicator_loc.Volume.tau_E);
    Indicator.tau_J = sqrt(Indicator_loc.Faces.tau_J);
    Indicator.tau_N = sqrt(Indicator_loc.Faces.tau_N);
    Indicator.tau_T = sqrt(Indicator_loc.Faces.tau_T);

    Indicator.tau   = sqrt(Indicator_loc.Volume.tau_E+Indicator_loc.Faces.tau_J+Indicator_loc.Faces.tau_T+Indicator_loc.Faces.tau_N);
end