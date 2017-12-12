function [ resid, jacob ] = discretize (P, Pold, dt, p0, phi0, b0, cr, cf, visc, kinit, dx, dy, nx, ny, conn_list, well_info)

[accum, d_accum] = accum_calcs(P, p0, phi0, b0, cr, cf);
accumold = accum_calcs(Pold, p0, phi0, b0, cr, cf);
[trans, d_trans] = trans_calcs(P, p0, cf, b0, visc);

ncells = nx .* ny;
n_vert = (nx - 1) .* ny;
nfaces = n_vert + ((ny - 1) .* nx);

kh_avg = 0.00112712 .* (dx ./ dy) .* ((2 .* (kinit(conn_list(:, 1)) .* kinit(conn_list(:, 2)))) ./ (kinit(conn_list(:, 1)) + kinit(conn_list(:, 2))));
pot_diff = P(conn_list(:, 2)) - P(conn_list(:, 1));
trans_avg = 0.5 .* ((trans(conn_list(:, 1)) + trans(conn_list(:, 2))));

d_DP_dl   = - ones(nfaces, 1);
d_DP_dr   = ones(nfaces, 1);

d_T_dl = 0.5 .* d_trans(conn_list(:,1));
d_T_dr = 0.5 .* d_trans(conn_list(:,2));

flux = kh_avg(1:end) .* pot_diff(1:end) .* trans_avg(1:end);
df_dpr = dt .* (kh_avg(1:end) .* ((d_T_dr(1:end) .* pot_diff(1:end)) + (trans_avg(1:end))));
df_dpl = dt .* (kh_avg(1:end) .* ((d_T_dl(1:end) .* pot_diff(1:end)) - (trans_avg(1:end))));

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Residual and Jacobian creation


resid = ((0.1781076 .* dx .* dy) .* (accum - accumold)) + (dt .* well_info);
jacob = zeros(ncells, ncells);

counter = 0;
for i = 1: ny
    for k = 1: nx
    counter = counter + 1;

    jacob(counter, counter) = (0.1781076 .* dx .* dy) .* d_accum(counter);

    end
end

for i = 1:nfaces
    c1 = conn_list(i ,1);
    c2 = conn_list(i, 2);
    
    resid(c1) = resid(c1) - dt .* flux(i);
    jacob(c1, c1) = jacob(c1, c1) - dt .* df_dpl(i);
    jacob(c1, c2) = jacob(c1, c2) - dt .* df_dpr(i);

    resid(c2) = resid(c2) + dt .* flux(i);
    jacob(c2, c1) = jacob(c2, c1) - dt .* df_dpl(i);
    jacob(c2, c2) = jacob(c2, c2) - dt .* df_dpr(i);
end

end