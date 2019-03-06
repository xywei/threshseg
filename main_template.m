%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ThreshSeg: An MBO-type Approach
%
% Reference: Dong Wang, Haohan Li, Xiaoyu Wei, and Xiao-Ping Wang (2016):
%            An Efficient Threshold Dynamics Method for Image Segmentation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note: This is NOT an excutable script by itself!
% To start, make a copy of this file and make sure input files are under 
% ./input/ relative to your script.

%% Clear workspace & Import libs

clear;
clc;
close all;

% Path to the library
addpath('./lib/');

%% Set parameters

% The most important ones
dt          = 0.01;           % \delta t
lambda      = .1;             % \lambda
n_phases    = 4;              % Number of phases
isGrayScale = false;          % Image will be converted to grayscale first if set to true
input_file  = 'flowers.jpg';                                 % Input files are under ./input
initial_contour_filename = 'flowers_init_4phase_recs.txt';   % Used only if read from file.

verbose = true;       % Whether to draw intermediate states as computing.
% Setting this to be true will display an animation of contour evolution
% when iteration.

update_energy = true; % Compute and store energy at each time step.

energy_hint_steps = [7 25 70 115 155]; % Intermedia steps to be plotted.
energy_hint_fill  = true;
% NOTE: If you are not sure which steps to plot (probably you are tuning parameters
% like dt or lambda), then comment this line and simply not to declare this variable.
% Also, this only works when update_energy is set to true.

% Resizing the image
% Set to 0 to use original size
% Note that M = height, N = width w.r.t. an image
resize_M = 0;
resize_N = 0;

% Adding Gaussian noise
noise_level = 0.0;

% Timestep with damping
n_damping_steps = 0;       % Apply damping to the first several steps, set to <=0 to turn it off.
damping_factor  = 0.9;

% Normalizing input images:
%  - Normalize.None:    Do nothing.
%  - Normalize.Max:     Divide all values with the max.
%  - Normalize.MinMax:  Linear map s.t. max=1 and min=0 after transform.
%  - Normalize.MeanVar: Linear map s.t. mean=0 and var=1 after transform.
normalization_method = Normalize.Max;

% Initial contours:
%  - Init.ReadRecsFromFile
%  - Init.ReadPolygonsFromFile
%  - Init.GrabRecsFromGUI
%  - Init.GrabPolygonsFromGUI
initialization_method = Init.GrabRecsFromGUI;

% Max number of iterations
MAXITER = 1e6;

%% Now run the segmenter
run_segmenter;

