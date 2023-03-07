% This function randomly mixes a sequence of numbers

function [nums] = random_integers(nums)

% Make sure the same sample is taken each time
%RandStream.setDefaultStream (RandStream('mt19937ar','seed',1)); 

ind = []; i = 0;

sample = length(nums);

while i < sample,
  a = round(sample*rand);
  if isempty(find(ind==a)) && a > 0, 
    i = i + 1; ind(i) = a; 
  end
end

nums = nums(ind);
