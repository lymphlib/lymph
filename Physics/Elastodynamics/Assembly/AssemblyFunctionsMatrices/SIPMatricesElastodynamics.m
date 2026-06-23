%> @file   SIPMatricesElastodynamics.m
%> @author Ilario Mazzieri, Mattia Corti
%> @date   12 May 2026
%> @brief  Global matrices structure for elastodynamics problem with SIP implementation.
%>
%==========================================================================
%> @section classSIPMatricesElastodynamics Class description
%==========================================================================
%> @brief           Global matrices structure for elastodynamics problem with SIP implementation.
%>
%> @param Matrices_loc  Matrices struct from assembly function
%>
%> @retval Matrices     Matrices struct containing global matrices for solver
%>                   
%==========================================================================

function [Matrices] = SIPMatricesElastodynamics(Matrices_loc)

    %% Mass matrix elastic
    Matrices.M_P_rho = matlab.internal.math.blkdiag(Matrices_loc.Volume.M1_P_rho_loc,Matrices_loc.Volume.M1_P_rho_loc);

    %% Projection matrix
    Matrices.MPrjP = matlab.internal.math.blkdiag(Matrices_loc.Volume.MPrjP_1_loc,Matrices_loc.Volume.MPrjP_1_loc);

    %% Damping matrix velocity
    Matrices.Dvel = matlab.internal.math.blkdiag(Matrices_loc.Volume.D1_loc,Matrices_loc.Volume.D1_loc);
    Matrices.Ddis = matlab.internal.math.blkdiag(Matrices_loc.Volume.C1_loc,Matrices_loc.Volume.C1_loc);

    %% Absorbing matrices
    Matrices.Svel = -[Matrices_loc.Faces.ABC_S1_loc, Matrices_loc.Faces.ABC_S2_loc; 
                      Matrices_loc.Faces.ABC_S3_loc, Matrices_loc.Faces.ABC_S4_loc];

    Matrices.Rdis = -[Matrices_loc.Faces.ABC_R1_loc, Matrices_loc.Faces.ABC_R2_loc; 
                      Matrices_loc.Faces.ABC_R3_loc, Matrices_loc.Faces.ABC_R4_loc];

    %% Elastic SIP-DG matrix
    V    = [Matrices_loc.Volume.V1_loc,    Matrices_loc.Volume.V2_loc;    
            Matrices_loc.Volume.V3_loc,    Matrices_loc.Volume.V4_loc];

    IT = [Matrices_loc.Faces.IT1_loc, Matrices_loc.Faces.IT2_loc; 
          Matrices_loc.Faces.IT3_loc, Matrices_loc.Faces.IT4_loc];

    S_P  = matlab.internal.math.blkdiag(Matrices_loc.Faces.S1_P_loc,Matrices_loc.Faces.S4_P_loc);

    Matrices.A_E = V + S_P - IT - transpose(IT);

end
