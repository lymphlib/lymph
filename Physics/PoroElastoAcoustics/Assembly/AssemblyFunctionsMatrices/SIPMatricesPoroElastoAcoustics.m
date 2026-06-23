%> @file   SIPMatricesPoroElastoAcoustics.m
%> @author Ilario Mazzieri, Mattia Corti
%> @date   29 May 2026
%> @brief  Global matrices structure for poroelastoacoustics problem with SIP implementation.
%>
%==========================================================================
%> @section classSIPMatricesPoroElastoAcoustics Class description
%==========================================================================
%> @brief           Global matrices structure for poroelastoacoustics problem with SIP implementation.
%>
%> @param Matrices_loc  Matrices struct from assembly function
%>
%> @retval Matrices     Matrices struct containing global matrices for solver
%>                   
%==========================================================================

function [Matrices] = SIPMatricesPoroElastoAcoustics(Matrices_loc)

    %% SIP matrices for the elastodynamics (from Physics/Elastodynamics)
    Matrices_Ela_loc.Volume = Matrices_loc.Volume.Ela;
    Matrices_Ela_loc.Faces  = Matrices_loc.Faces.Ela;
    [Matrices.Ela] = SIPMatricesElastodynamics(Matrices_Ela_loc);
    
    %% SIP matrices for the poroelasticity
    Matrices_Poro_loc.Volume = Matrices_loc.Volume.Poro;
    Matrices_Poro_loc.Faces  = Matrices_loc.Faces.Poro;
    [Matrices.Poro] = SIPMatricesPoroelasticity(Matrices_Poro_loc);

    %% SIP matrices for the acoustics
    Matrices_Acu_loc.Volume = Matrices_loc.Volume.Acu;
    Matrices_Acu_loc.Faces  = Matrices_loc.Faces.Acu;
    [Matrices.Acu] = SIPMatricesAcoustics(Matrices_Acu_loc);

    %% SIP matrices for the poroelasticity-elasticity coupling
    Matrices_PoroEla_loc.Faces  = Matrices_loc.Faces.PoroEla;
    [Matrices.PoroEla] = SIPMatricesPoroElastic(Matrices_PoroEla_loc);

    %% SIP matrices for the poroelasticity-acoustic coupling
    Matrices_PoroAcu_loc.Faces  = Matrices_loc.Faces.PoroAcu;
    [Matrices.PoroAcu] = SIPMatricesPoroAcoustics(Matrices_PoroAcu_loc); 

    %% SIP matrices for the elasticity-acoustic coupling
    Matrices_ElaAcu_loc.Faces  = Matrices_loc.Faces.ElaAcu;
    [Matrices.ElaAcu] = SIPMatricesElastoAcoustics(Matrices_ElaAcu_loc); 

end
