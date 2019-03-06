%% Clear workspace & Import libs

clear;
clc;
close all;
addpath('../lib/');

%% Set parameters

% The most important ones
dt          = 0.01;          % \delta t
lambda      = 0.003;         % \lambda
performance_test = false;     % Set to true to turn off all unecessary computations.
                             % This option overwrites all others below.

n_phases    = 2;             % Number of phases
isGrayScale = true;          % Image will be converted to grayscale first if set to true
input_file = 'ring.png';                                % Input files are under ./input
initial_contour_filename = 'ring_init_2phase_recs.txt'; % Used only if read from file.

% Additional functionalities
verbose = true;       % Whether to draw intermediate states.
update_energy = true; % Compute and store energy at each time step.
energy_hint_steps = [7 25 70 115 155]; % Intermedia steps to be plotted.
energy_hint_fill = false; % set to true to use contourf, false to use contour.
% If not set, the default is contourf.

resize_M = 0;
resize_N = 0;

% Adding Gaussian noise
noise_level = 0.0;

% Timestep with damping
n_damping_steps = 0;       % Apply damping to the first several steps, set to <=0 to turn it off.
damping_factor  = 0.9;

normalization_method = Normalize.Max;
initialization_method = Init.ReadRecsFromFile;

% Max number of iterations
MAXITER = 1e6;

%% Now run the segmenter
run_segmenter;

