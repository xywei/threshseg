%% Clear workspace & Import libs

clear;
clc;
close all;
addpath('../lib/');

%% Set parameters

% The most important ones
dt          = 0.03;           % \delta t
lambda      = .005;             % \lambda
n_phases    = 4;              % Number of phases
isGrayScale = false;          % Image will be converted to grayscale first if set to true
input_file  = 'flowers.jpg';                                 % Input files are under ./input
initial_contour_filename = 'flowers_init_4phase_plgs.txt';   % Used only if read from file.

% Additional functionalities
verbose = true;       % Whether to draw intermediate states.
update_energy = true; % Compute and store energy at each time step.

resize_M = 0;
resize_N = 0;

% Adding Gaussian noise
noise_level = 0.0;

% Timestep with damping
n_damping_steps = 0;       % Apply damping to the first several steps, set to <=0 to turn it off.
damping_factor  = 0.9;

normalization_method = Normalize.Max;
initialization_method = Init.GrabPolygonsFromGUI;

% Max number of iterations
MAXITER = 1e6;

%% Now run the segmenter
run_segmenter;

