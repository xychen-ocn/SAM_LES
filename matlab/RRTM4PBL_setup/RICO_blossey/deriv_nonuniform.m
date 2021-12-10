function [deriv] = deriv_nonuniform(z,f)

%   DERIV_NONUNIFORM computes the first derivative of a function on
%   a non-uniform grid.
%
%        [DFDZ] = DERIV_NONUNIFORM(Z,F) computes a second-order
%        accurate approximation to the first derivative of the
%        function represented by the vector F which gives the values
%        of that function at the points Z.
%
%        See DIFF, INTERP1.

if length(z(:,1)) == 1
  z = z';
end

if length(f(:,1)) == 1
  f = f';
end

N = length(z);
if N ~= size(f,1)
  error(['Length of coordinate vector should match first dimension ' ...
         'of array']);
end

% Here we find second-order approximations to the first
% derivative at the points zhalf (= (z(1) + z(2))/2, etc.).
zhalf = 0.5*(z(1:end-1) + z(2:end));
df = divide_first_dimension(diff(f),diff(z));

% Interpolate these values onto the original set of points using
% linear interpolation with extrapolation allowed for the initial
% and final points.
deriv = interp1(zhalf,df,z,'linear','extrap');
% $$$ deriv = interp1(zhalf,df,z,'spline');
