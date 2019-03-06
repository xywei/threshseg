classdef Init
    % Class for setting initial contours
    
    enumeration
        ReadRecsFromFile, ReadPolygonsFromFile, GrabRecsFromGUI, GrabPolygonsFromGUI
    end
    
    methods(Static)
        
        % Main interface
        function u = getInitialValues (I, n_phases, M, N, initialization_method, initial_contour_filename)
            if initialization_method==Init.ReadRecsFromFile
                u = Init.initReadRecsFromFile (n_phases, M, N, initial_contour_filename);
            elseif initialization_method==Init.ReadPolygonsFromFile
                u = Init.initReadPolygonsFromFile (n_phases, M, N, initial_contour_filename);
            elseif initialization_method==Init.GrabRecsFromGUI
                u = Init.initGrabRecsFromGUI (I, n_phases, M, N);
            elseif initialization_method==Init.GrabPolygonsFromGUI
                u = Init.initGrabPolygonsFromGUI (I, n_phases, M, N);
            else
                error ('Please specify correct initialization method!');
            end
        end
        
        
        % Read rectangles from a file representing initial regions.
        % Each of the subsequent n_phases-1 lines contains four real
        % numbers a, b, c, d in [0,1], representing a rectangular region of
        % [a,b]X[c,d]. (Note there is actually a transpose in the image.)
        %
        % If there is intersection, the intersection part is assigned to
        % the phase that is not specified (the background).
        % In the end the last phase fills what is left.
        function u = initReadRecsFromFile (n_phases, M, N, filename)
            fileID = fopen(['./input/' filename], 'r');
            display(['Using iniitals from file ./input/' filename '.']);
            formatSpec =  '%f %f %f %f';
            sizeA = [4 Inf];
            A = fscanf(fileID,formatSpec,sizeA);
            fclose(fileID);
            
            % Handle abnormal inputs
            if n_phases ~= size(A,2)+1
                error ('Invalid input file! (Number of rectangles must equal to n_phases-1!)');
            end
            if min(min(A)) < 0
                A = A + min(min(A));
            end
            if max(max(A)) > 1
                A = A ./ max(max(A));
            end
            
            u = zeros(M,N,n_phases);
            
            for phase = 1:n_phases-1
                lmt = floor(A(:,phase) .* [M M N N]');
                u(lmt(1):lmt(2), lmt(3):lmt(4), phase) = 1;
            end
            
            % Handle intersection
            intersection_mask = ones(M,N) - (sum(u,3) > 1);
            for phase = 1:n_phases-1
                u(:,:,phase) = min(u(:,:,phase), intersection_mask);
            end;
            
            % Background phase
            u(:,:,n_phases) = 1 - sum(u,3);
        end
        
        
        % Grab rectangles from GUI using a mouse
        % The hand picked contours can be saved and can be used
        % for future runs.
        function u = initGrabRecsFromGUI (I, n_phases, M, N)
            imshow(I);
            title (['Draw ' int2str(n_phases-1) ' rectangles with mouse ..'])
            hold on;
            
            u = zeros(M,N,n_phases);
            display('The initials are rectangles obtained through GUI.');
            display('You can copy and paste them into a file and reload for future reuse.');
            display('(Set initialization_method = Init.ReadRecsFromFile and');
            display(' initial_contour_filename correspondingly.');
            display(' ');
            display('>>> Start of initials <<<');
            
            for phase = 1:n_phases-1
                rect = getrect; % rect = [xmin, ymin, width, height]
                % lmt = [ymin, ymax, xmin, xmax] due to the transpose of
                % image show.
                rect = floor(rect);
                lmt = [rect(2), rect(2) + rect(4), rect(1), rect(1) + rect(3)];
                lmt = max(lmt, 1);
                % Handle out-of-image selections
                lmt(1:2) = min(lmt(1:2), M);
                lmt(3:4) = min(lmt(3:4), N);
                u(lmt(1):lmt(2), lmt(3):lmt(4), phase) = 1;
                contour(u(:,:,phase), [0.5 0.5], 'r', 'LineWidth',1.3);
                drawnow;
                ratios = lmt ./ [M M N N];
                display(num2str(ratios));
            end
            
            display('>>> End of initials <<<');
            display(' ');
            
            % Handle intersection
            intersection_mask = ones(M,N) - (sum(u,3) > 1);
            for phase = 1:n_phases-1
                u(:,:,phase) = min(u(:,:,phase), intersection_mask);
            end;
            
            % Background phase
            u(:,:,n_phases) = 1 - sum(u,3);
            
            close;
        end
        
        
        % Read polygons from a file representing initial regions.
        % Each of the subsequent n_phases-1 lines contains four real
        % numbers a, b, c, d in [0,1], representing a rectangular region of
        % [a,b]X[c,d]. (Note there is actually a transpose in the image.)
        %
        % If there is intersection, the intersection part is assigned to
        % the phase that is not specified (the background).
        % In the end the last phase fills what is left.
        function u = initReadPolygonsFromFile (n_phases, M, N, filename)
            A = dlmread(['./input/' filename]);
            display(['Using iniitals from file ./input/' filename '.']);
            
            % Handle abnormal inputs
            if n_phases ~= size(A,1)/2 + 1
                error ('Invalid input file! (Number of rectangles must equal to n_phases-1!)');
            end
            A = max(A, 0);
            A = min(A, 1);
            
            u = zeros(M,N,n_phases);
            x = (1:N) ./ N;
            y = (1:M) ./ M;
            [xx,yy] = meshgrid(x, y);
            
            for phase = 1:n_phases-1
                n_vertices = nnz( A(2*phase-1,:) + A(2*phase,:) );
                u(:,:,phase) = inpolygon(xx,yy, ...
                               A(2*phase-1,1:n_vertices),A(2*phase,1:n_vertices));
            end
            
            % Handle intersection
            intersection_mask = ones(M,N) - (sum(u,3) > 1);
            for phase = 1:n_phases-1
                u(:,:,phase) = min(u(:,:,phase), intersection_mask);
            end;
            
            % Background phase
            u(:,:,n_phases) = 1 - sum(u,3);
        end
        

        % Grab polygons from GUI using a mouse
        % The hand picked contours can be saved and can be used
        % for future runs.
        function u = initGrabPolygonsFromGUI (I, n_phases, M, N)
            imshow(I);
            title (['Draw ' int2str(n_phases-1) ' polygons with mouse ..'])
            hold on;
            
            u = zeros(M,N,n_phases);
            display('The initials are polygons obtained through GUI.');
            display('You can copy and paste them into a file and reload for future reuse.');
            display('(Set initialization_method = Init.ReadPolygonsFromFile and');
            display(' initial_contour_filename correspondingly.');
            display(' ');
            display('>>> Start of initials <<<');
            
            [xx,yy] = meshgrid(1:N, 1:M);
            for phase = 1:n_phases-1
                plg = getline('closed'); 
                plg = floor(plg);
                % Handle out-of-image selections
                plg = max(plg, 1);
                plg(:,1) = min(plg(:,1), N);
                plg(:,2) = min(plg(:,2), M);
            
                u(:,:,phase) = inpolygon(xx,yy,plg(:,1),plg(:,2));
                contour(u(:,:,phase), [0.5 0.5], 'r', 'LineWidth',1.3);
                drawnow;
                ratios = plg';
                ratios(1,:) = ratios(1,:) ./ N;
                ratios(2,:) = ratios(2,:) ./ M;
                display(num2str(ratios));
            end
            
            display('>>> End of initials <<<');
            display(' ');
            
            % Handle intersection
            intersection_mask = ones(M,N) - (sum(u,3) > 1);
            for phase = 1:n_phases-1
                u(:,:,phase) = min(u(:,:,phase), intersection_mask);
            end;
            
            % Background phase
            u(:,:,n_phases) = 1 - sum(u,3);
            
            close;
        end
        
    end
end