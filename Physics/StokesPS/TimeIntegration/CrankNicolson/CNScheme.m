%> @file  NewmarkScheme.m
%> @author Ilario Mazzieri
%> @date 16 Febraury 2024
%> @brief  Crank-Nicolson time integration scheme
%>
%==========================================================================
%> @section CNScheme Class description
%> @brief  Crank-Nicolson time integration scheme
%
%> @param Setup       Struct with problem's setup
%> @param Data        Struct with problem's data
%> @param femregion   Finite Element struct (see CreateDOF.m)
%> @param mesh        Struct with mesh geometry for output
%> @param M,K,F       Matrices and vector for MX'+ KX = F
%> @param Uold        Initial vector X(0)
%
%> @retval Solutions  Struct with solution vector
%>
%==========================================================================

function [Solutions] = CNScheme(Setup, Data, mesh, femregion, M, K, F, Uold)

% Create output folder for vtk files is it doesn't exist
if ~exist(Setup.OutFolderVTK,'dir'); mkdir(Setup.OutFolderVTK); end

% Time discretization parameters
dt  = Data.dt;
nts = floor(Data.T/dt);

% teta parameter for Crank-Nicolson
theta  = 0.5;

ML_Fluid = M + K*theta*dt;
MR_Fluid = M - K*(1-theta)*dt;

clear M K

Fold = F.f_f * Data.source_sigma_t{1}(0) + F.g_f * Data.source_sigma_d_t{1}(0);

t = 0;
disp(['Starting time: ', num2str(t)]);
disp('------------------------------------------------------')

tic;
counter = 1;

for t = dt : dt : nts*dt

    disp(['time: ', num2str(t)]);

    Fnew = F.f_f * Data.source_sigma_t{1}(t) + F.g_f * Data.source_sigma_d_t{1}(t);
    F_Fluid = MR_Fluid*Uold + dt*theta*Fnew + dt*(1-theta)*Fold;

    Uh = ML_Fluid\(F_Fluid);
    Solutions{counter} = SaveSolution(Uh,femregion);

    if (mod(counter,Data.timesave)==0)

        %% Save .mat file
        if Setup.isSaveSolution
            filename = fullfile(Setup.OutFolder,[Data.name,'_',num2str(counter),'.mat']);
            save(filename,'Data','dt','femregion','t','Solutions');
        end

        %% Save solutions paraview
        if (Setup.isSaveCSV || Setup.isSaveVTK || Setup.isPlotSolution)
            [Disp, Pressure] = GetSolutionQuadPointsPS(Data, femregion, Solutions, t);
        end

        if Setup.isSaveCSV
            csvwrite(fullfile(Setup.OutFolder,['Disp_',num2str(counter),'.csv']),Disp.Sigma);
            csvwrite(fullfile(Setup.OutFolder,['Pre_',num2str(counter),'.csv']),Pressure.Pf);
        end
        if Setup.isSaveVTK

            % pseudo stress
            tri = delaunay(Disp.Sigma(:,1),Disp.Sigma(:,2));
            fname = fullfile(Setup.OutFolderVTK, [Data.name,'_sigma11_', num2str(counter),'.vtk']);
            WriteVtk(fname, Disp.Sigma(:,1), Disp.Sigma(:,2), 0*Disp.Sigma(:,2), Disp.Sigma(:,3), tri, 's_11');
            fname = fullfile(Setup.OutFolderVTK, [Data.name,'_sigma12_', num2str(counter),'.vtk']);
            WriteVtk(fname, Disp.Sigma(:,1), Disp.Sigma(:,2), 0*Disp.Sigma(:,2), Disp.Sigma(:,4), tri, 's_12');
            fname = fullfile(Setup.OutFolderVTK, [Data.name,'_sigma21_', num2str(counter),'.vtk']);
            WriteVtk(fname, Disp.Sigma(:,1), Disp.Sigma(:,2), 0*Disp.Sigma(:,2), Disp.Sigma(:,5), tri, 's_21');
            fname = fullfile(Setup.OutFolderVTK, [Data.name,'_sigma22_', num2str(counter),'.vtk']);
            WriteVtk(fname, Disp.Sigma(:,1), Disp.Sigma(:,2), 0*Disp.Sigma(:,2), Disp.Sigma(:,6), tri, 's_22');
            fname = fullfile(Setup.OutFolderVTK, [Data.name,'_div_sigma_', num2str(counter),'.vtk']);
            WriteVtk(fname, Disp.divSigma(:,1), Disp.divSigma(:,2), 0*Disp.divSigma(:,2), [Disp.divSigma(:,3) Disp.divSigma(:,4)], tri, 'div_s');

        end

        %% Plot solution matlab
        if Setup.isPlotSolution
            ScatterPlot(Disp.Sigma(:,[1 2 3 7]),  mesh.region, Data, [Disp.SigmaTag,'11'], 0);
            ScatterPlot(Disp.Sigma(:,[1 2 4 8]),  mesh.region, Data, [Disp.SigmaTag,'12'], 0);
            ScatterPlot(Disp.Sigma(:,[1 2 5 9]),  mesh.region, Data, [Disp.SigmaTag,'21'], 0);
            ScatterPlot(Disp.Sigma(:,[1 2 6 10]), mesh.region, Data, [Disp.SigmaTag,'22'], 0);
            ScatterPlot(Pressure.Pf(:,[1 2 3 4]), mesh.region, Data, 'p', 0);

        end



    end

    Uold    = Uh;
    Fold    = Fnew;
    counter = counter + 1;

end

Solutions{counter} = SaveSolution(Uh,femregion);


toc







