function [f] = compute_data_term(g, u)
%
% Compute the data term in a vector fashion
%
% Input:
%  - g: M-by-N-by-n_channels array, the (preprocessed) image.
%  - u: M-by-N-by-n_phases array, the characteritic functions.
%
% Output:
%  - f: M-by-N-by-n_phases array, the data term f_i = f(:,:,i) described in the paper.

n_phases   = size(u, 3);
n_channels = size(g, 3);

% Mean
C = zeros(n_channels, n_phases);
for channel = 1:n_channels
    for phase = 1:n_phases
        if sum(sum(u(:,:,phase)))==0
            C(channel, phase) = 0;
        else
            C(channel, phase) = sum(sum(u(:,:,phase).*g(:,:,channel))) / sum(sum(u(:,:,phase)));
        end
    end
end

% Variance
f = zeros(size(u));
for phase = 1:n_phases
    for channel = 1:n_channels
        f(:,:,phase) = f(:,:,phase) + ( g(:,:,channel) - C(channel, phase) ).^2;
    end
end

end




