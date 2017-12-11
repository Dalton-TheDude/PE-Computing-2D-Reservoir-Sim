function [poro, d_poro] = poro_calcs( P, p0, phi0, cr )

poro = phi0 .* exp(cr .* (P - p0));
d_poro = cr .* poro;

end