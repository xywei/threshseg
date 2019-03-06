classdef Normalize
    % Class for vector normalization.
    % If input is a matrix, it is first reshaped to a vector
    
   enumeration
      None, Max, MinMax, MeanVar
   end
   
   methods(Static)
       
       % Main interface
       function out = apply (in, normalization_method)
           out = in;
           if normalization_method==Normalize.Max
               out = Normalize.normalizeMax (in);
           elseif normalization_method==Normalize.MinMax
               out = Normalize.normalizeMinMax (in);
           elseif normalization_method==Normalize.MeanVar
               out = Normalize.normalizeMeanVar (in);
           elseif normalization_method~=Normalize.None
               error ('Please specify correct normalization method!');
           end
       end
       
       function out = normalizeMax (in)
          [M,N,L] = size(in);
          vec_tmp = reshape(in, [M*N*L,1]);
          if abs(max(vec_tmp)) < 1e-12
              vec_tmp = vec_tmp ./ 1e-12;
          else
              vec_tmp = vec_tmp ./ max(vec_tmp);
          end
          out = reshape(vec_tmp, [M,N,L]);
       end
       function out = normalizeMinMax (in)
           [M,N,L] = size(in);
           vec_tmp = reshape(in, [M*N*L,1]);
           mn = min (vec_tmp);
           mx = max (vec_tmp);
           vec_tmp = (vec_tmp - mn) ./ (mx - mn);
           out = reshape(vec_tmp, [M,N,L]);
       end
       function out = normalizeMeanVar (in)
           [M,N,L] = size(in);
           vec_tmp = reshape(in, [M*N*L,1]);
           mn = mean (vec_tmp);
           sd = std  (vec_tmp);
           vec_tmp = (vec_tmp - mn) ./ sd;
           out = reshape(vec_tmp, [M,N,L]);
       end
   end
end