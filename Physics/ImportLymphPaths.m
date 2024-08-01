restoredefaultpath;

%% Import the necessary lymph paths for the problem resolution and print lymph header

% LymphFolder is the parent of the one containing this file
% (LymphFolder/Physics/ImportLymphPaths.m)
[LymphFolder, ~, ~] = fileparts(fileparts(mfilename('fullpath')));
PathCore    = fullfile(LymphFolder,'Core');

% Check if Polymesher is available and possibly download it.
PathPolymesher = fullfile(PathCore,'/MeshGeneration/PolyMesh/Polymesher');
PathOriginalPolyMesher = fullfile(PathPolymesher, 'SMO_12_PolyMesher_matlabcode_v1.1');
if ~exist(PathOriginalPolyMesher, 'dir')
	currentPath = pwd;
	cd(PathPolymesher);
	
	% Download and extract the original PolyMesher code - version 1.1
	websave(strcat(PathOriginalPolyMesher,'.zip'),'http://paulino.princeton.edu/journal_papers/2012/SMO_12_PolyMesher_matlabcode_v1.1.zip');
	unzip(strcat(PathOriginalPolyMesher,'.zip'), PathOriginalPolyMesher);
	delete(strcat(PathOriginalPolyMesher,'.zip'));
	
	% Delete original version of PolyMesher.m: lymph has its own modified implementation.
	delete(fullfile(PathOriginalPolyMesher, 'PolyMesher', 'PolyMesher.m'));

	cd(currentPath);
end

addpath(genpath(PathCore));

Header;