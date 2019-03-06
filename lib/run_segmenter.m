%% The main procedure of segmentation
% This script must be called in such a way that the
% input images and intial conditions are saved under ./input/.
% The path above is RELATIVE to the path where this is called.

if exist('performance_test','var')
    if performance_test
        verbose = false;
        update_energy = false;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

display('*--------------------------------------------------------------------------');
display('* THRESH-SEG: An MBO-type Approach');
display('* Reference: Dong Wang, Haohan Li, Xiaoyu Wei, and Xiao-Ping Wang (2016):');
display('*            An efficient iterative thresholding method for image segmentation.');
display('*--------------------------------------------------------------------------');

I = imread(['./input/' input_file]);

% Resize
if resize_M > 0
    I = imresize(I,[resize_M,resize_N]);
end;
g = im2double(I);
if isGrayScale
    g = rgb2gray(g);
end
[M,N,n_channels]=size(g);

% Add noise
g = g + noise_level .* randn(size(g));

% Normalize
g = Normalize.apply(g, normalization_method);

% Initial Contours
% In this implementation, u[:,:,1] ~ u[:,:,n_phase-1] are foreground shapes
u = Init.getInitialValues(I, n_phases, M, N, initialization_method, initial_contour_filename);

if verbose
    figure();
    
    % Define color order
    ColOrd   = get(gca,'ColorOrder');
    n_colors = size(ColOrd, 1);
    
    subplot(2,2,1); 
    imshow(I); 
    title('Input Image');
    
    subplot(2,2,2); 
    imshow(I);
    hold on;
    for phase = 1:n_phases-1
        Col = ColOrd(mod(phase,n_colors)+1,:);
        contour(u(:,:,phase), [0.5 0.5], 'Color', Col, 'LineWidth',1.3); 
    end
    title('Initial Contour');
    subplot(2,2,3); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main Iteration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist('energy_hint_steps','var')
    hint = 0;
    u_cache = zeros([M,N,n_phases,length(energy_hint_steps)]);
end
k=0;

fprintf('Iterating ..');
is_moving   = true;
tic;
while is_moving && k < MAXITER
    if verbose
        imshow(g);
        hold on;
        for phase = 1:n_phases
            contour(u(:,:,phase), [0.5 0.5], 'r','LineWidth',1.3); 
        end
        hold off;
        title('Segmentation');
        drawnow;
    end

    k = k+1;
    % display (['Step ' int2str(k)]);
    if mod(k,50) == 1
        fprintf('\n');
    end
    fprintf('.');
    
    % Cache the results for later use in energy plot
    if update_energy
        if exist('energy_hint_steps','var')
            if any(abs(k-energy_hint_steps)<1e-10)
                hint = hint + 1;
                u_cache(:,:,:,hint) = u;
            end
        end
    end
    
    time_iter_start = toc;
    
    % Update data term
    f = compute_data_term (g, u);
    
    % Apply heat kernel convolution
    uh = apply_heat_convolution (dt, u);
    
    % Thresholding: for each pixel, assign it to the phase with lowest
    % score. If there is a tie, the phase with smaller index wins.
    % TODO: is there a faster way to implement this part? (perhaps mex C?)
    is_moving = false;
    ttmp = toc;
    for j=1:N
        for i=1:M
            % Get scores for each phase
            min_score = Inf;
            for phase = 1:n_phases
                if f(i,j,phase) - 2 * lambda / sqrt(dt) * sqrt(pi) * uh(i,j,phase) < min_score
                    min_score = f(i,j,phase) - 2 * lambda / sqrt(dt) * sqrt(pi) * uh(i,j,phase);
                    min_index = phase;
                end;
            end
            if u(i,j,min_index)==0
                is_moving = true;
            end
            u(i,j,:) = 0;
            u(i,j,min_index) = 1;
        end
    end
    
    % Damping
    if k < n_damping_steps
        dt = dt * damping_factor;
    end
    
    % Update energy
    % This may be slow because the energy array is not pre-allocated.
    % As in the proof, u_i is at step k+1,
    %  data term and G*u_i is at step k.
    if update_energy
        energy(k) = compute_energy (u, uh, f, dt, lambda);
    end
    
end

stop_time = toc;

if verbose
    subplot(2,2,4)
    imshow(I);
    hold on;
    cartoon = zeros(M,N);
    for phase = 1:n_phases
        cartoon = cartoon + phase * u(:,:,phase);
    end
    contourf(cartoon);
    hold off;
    title 'The Cartoonization';
end

fprintf('\n');
display (['Convergence achieved in ' int2str(k-1) ' steps, runtime = ' num2str(stop_time) '.']);

% Produce energy plot
if update_energy
    display ('Energy history is stored in variable "energy".');
    % Size of energy plot
    p1 = 0; p2 = 0; width = 12; height = 6;
    % Size of embeded segmentations
    sw = 1.3; sh = 1.3; offset1 = 0.05; offset2 = 0.05;
    figure('Units','Inch','Position',[p1 p2 width height]);
    plot(energy, 'LineWidth', 2);
    set(gca,'FontSize',18);
    hold on;
    if exist('energy_hint_steps','var')
        plot(energy_hint_steps, energy(energy_hint_steps), 'd', ...
            'MarkerEdgeColor','k', 'MarkerFaceColor', 'r', 'MarkerSize', 10);
    end
    xlabel 'Time Step Number';
    ylabel 'Energy';
    xl = xlim; 
    yl = ylim;
    AxesHandle=findobj(gcf,'Type','axes');
    geometry = get(AxesHandle,'Position');
    if exist('energy_hint_steps','var')
        hint = 1;
        while hint <= length(energy_hint_steps)
            subax1 = width * (geometry(1) + geometry(3) * (energy_hint_steps(hint) - xl(1)) / (xl(2)-xl(1)));
            subax2 = height * (geometry(2) + geometry(4) * (energy(energy_hint_steps(hint)) - yl(1)) / (yl(2) - yl(1)));
            axes('Units','Inch','Position', [subax1+offset1 subax2+offset2 sw sh]);
            cartoon = zeros(M,N);
            for phase = 1:n_phases
                cartoon = cartoon + phase * u_cache(:,:,phase,hint);
            end
            if exist('energy_hint_fill','var')
                if energy_hint_fill
                    imshow(I);
                    hold on;
                    contourf(cartoon);
                else
                    [~,hd] = imcontour(cartoon, 1:n_phases-1, 'k');
                    set(hd,'LineWidth',2)
                    % No box
                    axis off;
                    % With box
                    %set(gca,'YTick',[]);
                    %set(gca,'XTick',[]);              
                end
            else
                imshow(I);
                hold on;
                contourf(cartoon);
            end
            if isGrayScale
                 set(gcf, 'colormap', gray); 
            end
            h = title(int2str(energy_hint_steps(hint)));
            P = get(h, 'Position');
            set(h,'Position',[P(1) P(2)+35 P(3)]);
            set(gca,'FontSize',18);
            hint = hint + 1;
        end
    end
end