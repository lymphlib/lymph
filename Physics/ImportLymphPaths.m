restoredefaultpath;

%% Import the necessary lymph paths for the problem resolution and print lymph header

% LymphFolder is the parent of the one containing this file
% (LymphFolder/Physics/ImportLymphPaths.m)
[LymphFolder, ~, ~] = fileparts(fileparts(mfilename('fullpath')));
PathCore    = fullfile(LymphFolder,'Core');

addpath(genpath(PathCore));

Header;