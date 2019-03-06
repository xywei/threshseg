function [uh] = apply_heat_convolution(dt, u)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fast evaluation of convelution using FFT by the basic principle
% g*h=sum(\frac{1}{dx*dy}iFFT(FFT(g).FFT(h)))
% We may do some extension to make the image to be perodic in both x, y
% direction, That depends on the profile of the image and can be modified
% if needed. In this simple code, we just assume the nonzero part are away
% from the boundary of computational domain.
% Input: dt --- artificial time step
%        u  --- characteristic function of different regions.
% Output:u_hat --- diffused value of the charecteristic function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get size info
[M,N,n_phases] = size(u);

% Precompute the weights.
% Computational domain = [-pi, pi]^2.
x1 = (0:M-1)*2*pi/M-pi;
x2 = (0:N-1)*2*pi/N-pi;
dx1 = x1(2)-x1(1);
dx2 = x2(2)-x2(1);
xx1 = repmat(x1,N,1);
xx2 = repmat(x2,M,1)';
G_dt = (1/(4*pi*dt)*exp(-(xx1.^2+xx2.^2)/(4*dt)))';
K = fft2(G_dt);

% Apply convolution via FFT
uh = zeros(size(u));
for phase = 1:n_phases
    uh(:,:,phase) = real( dx1 * dx2 * ifftshift(ifft2(fft2(u(:,:,phase)) .* K)) );
end

end