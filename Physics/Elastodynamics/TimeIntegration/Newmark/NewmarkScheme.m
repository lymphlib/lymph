%> @file  NewmarkScheme.m
%> @author Ilario Mazzieri
%> @date 11 May 2023
%> @brief  Newmark time integration scheme
%>
%==========================================================================
%> @section classNewmarkScheme Class description
%> @brief  Newmark time integration scheme
%
%> @param Setup       Struct with problem's setup
%> @param Data        Struct with problem's data
%> @param femregion   Finite Element struct (see CreateDOF.m)
%> @param mesh        Struct with mesh geometry for output
%> @param A,B,C,F     Matrices and vector for AX''+ BX'+ CX = F
%> @param Uold        Initial vector [X(0), X'(0)]
%
%> @retval Solutions  Struct with solution vector
%>
%==========================================================================

function [Solutions] = NewmarkScheme(Setup, Data, femregion, mesh, A, B, C, F, Uold)

% Time discretization parameters
dt  = Data.dt;
nts = floor(Data.T/dt);

% Newmark parameter
beta_nm  = Data.BetaNM; 
gamma_nm = Data.GammaNM;

A1 = [A + dt^2*beta_nm*C, dt^2*beta_nm*B; ...
    gamma_nm*dt*C,      A + gamma_nm*dt*B];

A2 = [A-dt^2*(1/2-beta_nm)*C, dt*A - dt^2*(1/2-beta_nm)*B; ...
    -dt*(1-gamma_nm)*C,     A - dt*(1-gamma_nm)*B];

clear A B C


Fold = F.f_e * Data.source_ue_t{1}(0) + F.g_e * Data.source_ued_t{1}(0);


t = 0;
disp(['Starting time: ', num2str(t)]);
disp('------------------------------------------------------')

tic;
counter = 1;

for t = dt : dt : nts*dt
    
    disp(['time: ', num2str(t)]);
    
    Fnew = F.f_e * Data.source_ue_t{1}(t) + F.g_e * Data.source_ued_t{1}(t);
    
    rhs = [dt^2*beta_nm*Fnew + dt^2*(1/2-beta_nm)*Fold; ...
        gamma_nm*dt*Fnew  + dt*(1-gamma_nm)*Fold];
    
    Uh = A1\(A2*Uold + rhs);
    
    if (mod(counter,Data.timesave)==0)
        
        [Solutions] = SaveSolution(Uh,femregion);
        
        %% Save .mat file        
        if Setup.isSaveSolution
            filename = fullfile(Setup.OutFolder,[Data.name,'_',num2str(counter),'.mat']);
            save(filename,'Data','dt','femregion','t','Solutions');
        end
        
        %% Save solutions paraview
        if (Setup.isSaveCSV || Setup.isSaveVTK || Setup.isPlotSolution)
            [Displacement, Pressure, Velocity] = GetSolutionQuadPoints(Data, femregion, Solutions, t);
        end
        
        if Setup.isSaveCSV
            csvwrite(fullfile(Setup.OutFolder,['Disp_',num2str(counter),'.csv']),Displacement.Ue);
            csvwrite(fullfile(Setup.OutFolder,['Pre_',num2str(counter),'.csv']),Pressure.Pe);
            csvwrite(fullfile(Setup.OutFolder,['Vel_',num2str(counter),'.csv']),Velocity.Ue);
        end
        if Setup.isSaveVTK
            tri = delaunay(Displacement.Ue(:,1),Displacement.Ue(:,2));
            propertyName = 'ue';
            fname = fullfile(Setup.OutFolderVTK, [Data.name,'_',propertyName,'_', num2str(counter),'.vtk']);
            WriteVtk(fname, Displacement.Ue(:,1), Displacement.Ue(:,2), 0*Displacement.Ue(:,2), [Displacement.Ue(:,3) Displacement.Ue(:,4)], tri, propertyName);
            propertyName = 've';
            fname = fullfile(Setup.OutFolderVTK, [Data.name,'_',propertyName,'_', num2str(counter),'.vtk']);
            WriteVtk(fname, Velocity.Ue(:,1), Velocity.Ue(:,2), 0*Velocity.Ue(:,2), [Velocity.Ue(:,3) Velocity.Ue(:,4)], tri, propertyName);
            propertyName = 'pe';
            fname = fullfile(Setup.OutFolderVTK, [Data.name,'_',propertyName,'_', num2str(counter),'.vtk']);
            WriteVtk(fname, Pressure.Ue(:,1), Pressure.Ue(:,2), 0*Pressure.Ue(:,2), Pressure.Ue(:,3), tri, propertyName);
        end
        
        %% Plot solution matlab
        if Setup.isPlotSolution
           ScatterPlot(Displacement.Ue, mesh, Data, Displacement.UeTag, 1);
           ScatterPlot(Displacement.Ue, mesh, Data, Displacement.UeTag, 2);
           ScatterPlot(Velocity.Ue, mesh, Data, Velocity.UeTag, 1);
           ScatterPlot(Velocity.Ue, mesh, Data, Velocity.UeTag, 2);
        end
        
        
           
        
        
    end
    
    Uold    = Uh;
    Fold    = Fnew;
    counter = counter + 1;
    
end

[Solutions] = SaveSolution(Uh,femregion);

