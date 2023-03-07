% Function for performing batching of observations
%
%  INPUT 
%           C: matrix determining localization for obs-space  (Ny x Ny)
%      cutoff: cutoff localization coefficient for batching
%
%  OUTPUT 
%          yi: observation index after batching

function [yi] = batch_obs(C,cutoff)

% Get obs dimension
Ny = length(C(1,:));

% Initialize sampling and indexing vectors
yi = zeros(Ny,1);

for i = 1:Ny

  % Check to see if ob is already batched
  if yi(i) ~= 0, continue, end

  % Locate obs in neighborhood of current one
  ind = find(C(i,:)>=cutoff);

  % Search within obs
  for j = ind
    if yi(j) == 0
      % Store ob location
      yi(j) = i;
    end
  end   

end
