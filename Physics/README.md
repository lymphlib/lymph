# Solvers based on lymph #

This directory contains solvers for different physics.

It also provides the following files:
* `ImportLymphPaths.m`
	Set up the environment paths to use `lymph`: call `ImportLymphPaths` in your `Run*.m` file (see below).
* `RunSetup.m`
	Define common parameters for post-processing and output: to use these parameters, call `RunSetup` in your `Run*.m` file (see below).

When contributing to lymph with a new physics solver, **add a new directory here** with your code.
The typical structure of such new directory (named `MYPHYSICS` below) is as follows (see the `Heat` folder for an example):

* `Run*MYPHYSICS.m`
	Callable script(s) to initialize `lymph`, set up the I/O parameters of the simulation, and launch the execution of a simulation (see also next point).
	Multiple files with such name can be present, e.g. to run single simulations, convergence analysis, etc.
	This/these file(s) **must include the following lines**:
	```
	%% Import lymph and add path related to this physics.
	run("../ImportLymphPaths.m")
	MyPhysicsPath = pwd;
	addpath(genpath(fullfile(ProblemPath,'Assembly')));
	addpath(genpath(fullfile(ProblemPath,'InputData')));
	... and same for all other directories related to MYPHYSICS ...

	%% Simulation - Setup
	run("../RunSetup.m")
	```
* `MainFunctions`
	Folder containing the actual main function of a simulation, launched by `Run*MYPHYSICS.m`.
* `Assembly`
	Folder containing the functions to assemble the system of interest.
* `Error`
	Folder containing the functions to compute the error w.r.t. an exact solution.
* `InputData`
	Folder containing functions to generate a `Data` struct, used throughout the code, encoding the physical and numerical parameters of the problem to solve and possibly the exact solution for convergence analysis.
* `PostProcessing`
	Folder containing physics-specific functions for post-processing: see also `../Core/PostProcessing` for physics-independent functions (e.g. conversion to VTK format).
* Other folders containing utilities and functions to address specific issues of a physics solver.
* The folders related to `lymph` tutorials distributed with the current release contain also a `Doxygen` folder with their documentation.
