function [fvf, d_fvf] = fvf_calcs(P, p0, cf, b0)

y = cf .* (P - p0);
fvf = b0 ./ (1 + y + (0.5 .* y .* y));
% d_fvf = (- (b0 .* 1 + cf + (cf .* (P - p0)))) ./ ((1 + y + (0.5 .* y .* y)));
d_fvf = -((b0 .* cf .* (1+y)) ./ ((1+y+(0.5 .* (y .^2))) .^ 2));

end