function [zL, zV] = select_z_roots(z_roots, B_coef)
% select_z_roots  Robust selection of liquid and vapor Z-roots from cubic EOS.
%
%   [zL, zV] = select_z_roots(z_roots, B_coef)
%
%   Filters the three roots returned by roots() for physical values:
%   - Rejects roots with significant imaginary parts (ratio > 1e-6)
%   - Rejects roots with real part <= B_coef (required for log(Z-B) > 0)
%   Returns the smallest physical root as zL and largest as zV.
%   When only one physical root exists, zL == zV (supercritical / single phase).
%
% PARAMETERS:
%   z_roots  - [3x1] complex vector from roots(cubic_poly)
%   B_coef   - reduced co-volume B = b*P/(R*T); physical Z must exceed this
%
% RETURNS:
%   zL  - liquid compressibility factor (smallest physical root)
%   zV  - vapor compressibility factor (largest physical root)

tol = 1e-6;
real_z = real(z_roots);
% Root is "real" if |imag|/max(|real|,1) < tol
is_real = abs(imag(z_roots)) ./ max(abs(real_z), 1) < tol;
% Root is physical if real part exceeds co-volume
is_physical = is_real & (real_z > B_coef);
phys_z = sort(real_z(is_physical));

if numel(phys_z) >= 2
    zL = phys_z(1);
    zV = phys_z(end);
elseif numel(phys_z) == 1
    zL = phys_z(1);
    zV = zL;
else
    % Fallback: take root with smallest imaginary part
    [~, idx] = min(abs(imag(z_roots)));
    zL = real(z_roots(idx));
    if zL <= B_coef
        zL = B_coef + 1e-8;
    end
    zV = zL;
end
