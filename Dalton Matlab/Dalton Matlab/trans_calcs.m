function [ trans, d_trans ] = trans_calcs(P, p0, cf, b0, visc)

[fvf, d_fvf] = fvf_calcs(P, p0, cf, b0);

trans = 1 ./ (visc .* fvf);
d_trans = ((- d_fvf .* visc)) ./ ((visc .* fvf) .^ 2);

end