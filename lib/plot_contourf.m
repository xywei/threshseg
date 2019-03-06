function [ierr] = plot_contourf(u, isGrayScale)
% Plot the contour of segmentations
% Input:
%  - u: M-by-N-by-n_phases array, the characteritic functions.
%  - isGray: (Optional): Grayscale or color (bool).

if nargin < 2
   isGrayScale = false;
end

[M,N,n_phases] = size(u);

p1 = 0;
p2 = 0;

width = N;
height = M;

figure('Units','Pixel','Position',[p1 p2 width height]);

cartoon = zeros(M,N);
for phase = 1:n_phases
    cartoon = cartoon + phase * u(:,:,phase);
end
hold on;
contourf(cartoon);
set(gca,'position',[0 0 1 1],'units','normalized');
set(gca,'ydir','reverse');
if isGrayScale
    set(gcf, 'colormap', gray);
end
% No box
axis off;
% With box
%set(gca,'YTick',[]);
%set(gca,'XTick',[]);

ierr = 0;

end