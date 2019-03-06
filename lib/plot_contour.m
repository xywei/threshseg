function [ierr] = plot_contour(u, LineSpec, LineWidth)
% Plot the contour of segmentations
% Input:
%  - u: M-by-N-by-n_phases array, the characteritic functions.
%  - LineSpec (Optional): Line style specifier string.
%  - LineWidth (Optional): Line width (double).

if nargin < 3
   LineWidth = 2;
end
if nargin < 2
   LineSpec = 'k';
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
[~,hd] = imcontour(cartoon, 1:n_phases-1, LineSpec);
set(hd,'LineWidth',LineWidth)
% No box
axis off;
% With box
%set(gca,'YTick',[]);
%set(gca,'XTick',[]);
set(gca,'position',[0 0 1 1],'units','normalized')

ierr = 0;

end