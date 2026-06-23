%> @file  SimulationSetup.m
%> @author The Lymph Team
%> @date 13 May 2026
%> @brief Initialization of the simulation 
%>
%==========================================================================
%> @section classSimulationSetup Class description
%==========================================================================
%> @brief            Initialization of the simulation 
%>
%> @param ~
%>
%> @retval ~
%>
%==========================================================================

%% Import lymph and add path related to this physics.
run("ImportLymphPaths.m")
addpath(genpath(MyPhysicsPath));

%% Print the header
Header;

%% Create the setup structure
run("RunSetup.m")

%% Initiate the parallel
Setup.p = gcp('nocreate');
if Setup.isParallel && isempty(Setup.p)
    parpool(Setup.numCores);
elseif not(Setup.isParallel) && ~isempty(Setup.p)
    delete(Setup.p)
end

ps = parallel.Settings;
ps.Pool.AutoCreate = false;
