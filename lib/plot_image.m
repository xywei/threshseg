function [ierr] = plot_image(I)
% Plot the image
% Input:
%  - I: M-by-N-by-n_channels array, the image.

[M,N,~] = size(I);

p1 = 0;
p2 = 0;

width = N;
height = M;

figure('Units','Pixel','Position',[p1 p2 width height]);
imshow(I,'Border','tight');

ierr = 0;

end