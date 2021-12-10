function [quotient] = divide_first_dimension(numerator,denominator)

% Here, denomimator is a vector whose length matches the first
% dimension of numerator, which is a general array.

quotient = zeros(size(numerator));
if length(find(size(numerator)>1)) == 1
  quotient = numerator./denominator;
elseif length(size(numerator)) == 2
  for i = 1:length(denominator)
    quotient(i,:) = numerator(i,:)/denominator(i);
  end
elseif length(size(numerator)) == 3
  for i = 1:length(denominator)
    quotient(i,:,:) = numerator(i,:,:)/denominator(i);
  end
elseif length(size(numerator)) == 4
  for i = 1:length(denominator)
    quotient(i,:,:,:) = numerator(i,:,:,:)/denominator(i);
  end
else
  error('Add command to divide_first_dimension for matrix dimension > 4');
end  
