function [accum, d_accum] = accum_calcs(P, p0, phi0, b0, cr, cf)

[phi, d_phi] = poro_calcs(P, p0, phi0, cr);
[fvf, d_fvf] = fvf_calcs(P, p0, cf, b0);

accum = fvf ./ phi;
d_accum = ((phi .* d_fvf) - (fvf .* d_phi)) ./ (phi .* phi);

end