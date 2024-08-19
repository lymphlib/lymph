%> @file  MatPoroAcuEla.m
%> @author Ilario Mazzieri
%> @date 17 June 2024
%> @brief  Assembly of the matrices for the poro-elasto-acoustic problem
%>
%==========================================================================
%> @section classMatPoroAcuEla Class description
%> @brief  Assembly of the matrices for the poro-elasto-acoustic problem
%
%> @param Data       Struct with problem's data
%> @param neighbor   Neighbor struct (see MakeNeighbor.m)
%> @param femregion  Finite Element struct (see CreateDOF.m)
%
%> @retval Matrices.  Matrices.Poro    = Matrices for poroelastic domain
%> Matrices.Ela     = Matrices for elastic domain
%> Matrices.Acu     = Matrices for acoustic domain
%> Matrices.PoroAcu = Matrices for poro-acoustic coupling
%> Matrices.PoroEla = Matrices for poro-elastic coupling
%> Matrices.ElaAcu  = Matrices for elasto-acoustic
%coupling
%>                   
%==========================================================================

function [Matrices] = MatPoroAcuEla(Data, neighbor, femregion)

%% BULK MATRICES
% Porous media
disp('Matrices for Porous media ... ');
 
switch Data.quadrature
    case "QF"
        [Matrices.Poro] = MatPoroQF(Data, neighbor, femregion);
    case "ST"
        [Matrices.Poro] = MatPoroST(Data, neighbor, femregion);
end

disp('Done')
disp('------------------------------------------------------')

% Acoustic media
disp('Matrices for Acoustic media ... ');
switch Data.quadrature
    case "QF"
        [Matrices.Acu] = MatAcuQF(Data, neighbor, femregion);
    case "ST"
        [Matrices.Acu] = MatAcuST(Data, neighbor, femregion);
end
disp('Done')
disp('------------------------------------------------------')

% Elastic media
disp('Matrices for Elastic media ... ');
switch Data.quadrature
    case "QF"
        [Matrices.Ela] = MatElaQF(Data, neighbor, femregion);
    case "ST"
        [Matrices.Ela] = MatElaST(Data, neighbor, femregion);
end
disp('Done')
disp('------------------------------------------------------')

%% COUPLING MATRICES
% Porous-Acoustic media
disp('Matrices for Poro-Acoustic coupling ... ');
[Matrices.PoroAcu] = MatPoroAcu(Data, neighbor, femregion);
disp('Done')
disp('------------------------------------------------------')

% Porous-Elastic media
disp('Matrices for Poro-Elastic coupling ... ');
[Matrices.PoroEla] = MatPoroEla(Data, neighbor, femregion);
disp('Done')
disp('------------------------------------------------------')

% Elastic-Acoustic media
disp('Matrices for Elasto-Acoustic coupling ... ');
[Matrices] = MatElaAcu(Data, neighbor, femregion, Matrices);
disp('Done')
disp('------------------------------------------------------')

