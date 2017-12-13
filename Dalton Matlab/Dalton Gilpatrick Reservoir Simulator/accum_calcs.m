function [accum, d_accum] = accum_calcs(P, p0, phi0, b0, cr, cf)

[phi, d_phi] = poro_calcs(P, p0, phi0, cr);
[fvf, d_fvf] = fvf_calcs(P, p0, cf, b0);

accum = phi ./ fvf;
d_accum = (d_phi ./ fvf) - (phi ./ fvf ./ fvf .* d_fvf);

end