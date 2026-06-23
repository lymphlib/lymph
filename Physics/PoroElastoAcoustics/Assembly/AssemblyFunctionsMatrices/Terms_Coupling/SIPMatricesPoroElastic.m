%> @file   SIPMatricesPoroElastic.m
%> @author Ilario Mazzieri, Mattia Corti
%> @date   5 June 2026
%> @brief  Global matrices structure for poroelasticity-elastodynamics coupling.
%>
%==========================================================================
%> @section classSIPMatricesPoroElastic Class description
%==========================================================================
%> @brief           Global matrices structure for poroelasticity-elastodynamics coupling.
%>
%> @param Matrices_loc  Matrices struct from assembly function
%>
%> @retval Matrices     Matrices struct containing global matrices for solver
%>                   
%==========================================================================

function [Matrices] = SIPMatricesPoroElastic(Matrices_loc)

    %% Matrices
    Matrices.C_P = [Matrices_loc.Faces.C1_P_loc  Matrices_loc.Faces.C2_P_loc;
                    Matrices_loc.Faces.C3_P_loc  Matrices_loc.Faces.C4_P_loc];

    Matrices.C_E = [Matrices_loc.Faces.C1_E_loc'  Matrices_loc.Faces.C2_E_loc'; 
                    Matrices_loc.Faces.C3_E_loc'  Matrices_loc.Faces.C4_E_loc'];
    
end
