function [energy] = compute_energy (u, uh, f, dt, lambda)

% Compute the energy value.
% Input:
%  - u:  characteritic functions.
%  - uh: convolution results.
%  - f:  data term values.
%  - lambda: the model parameter.

[M,N,n_phases] = size(u);
weight = (2*pi/N) * (2*pi/M);

energy = 0;

for phase = 1:n_phases
    energy = energy + ...
             sum(sum( u(:,:,phase) .* f(:,:,phase) + ...
                      2 * lambda / sqrt(dt) * sqrt(pi) * (1-uh(:,:,phase)) .* u(:,:,phase) ...
                    )) * weight;
end

end
