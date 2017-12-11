function [ resid, jacob ] = discretize (P, Pold, dt, p0, phi0, b0, cr, cf, visc, kinit, dx, dy, nx, ny, conn_list, well_info)

[accum, d_accum] = accum_calcs(P, p0, phi0, b0, cr, cf);
accumold = accum_calcs(Pold, p0, phi0, b0, cr, cf);
[trans, d_trans] = trans_calcs(P, p0, cf, b0, visc);

ncells = nx .* ny;
n_vert = (nx - 1) .* ny;
nfaces = n_vert + ((ny - 1) .* nx);
dt_dx(1:n_vert, 1) = dt ./ (dx .^ 2);
dt_dy(n_vert+1:nfaces, 1) = dt ./ (dy .^ 2);

kh_avg = 0.00112712 .* (dx ./ dy) .* ((2 .* (kinit(conn_list(:, 1)) .* kinit(conn_list(:, 2)))) ./ (kinit(conn_list(:, 1)) + kinit(conn_list(:, 2))));
pot_diff = P(conn_list(:, 2)) - P(conn_list(:, 1));
trans_avg = 0.5 .* ((trans(conn_list(:, 1)) + trans(conn_list(:, 2))));

d_DP_dl   = - ones(nfaces);
d_DP_dr   = ones(nfaces);

d_T_dl = 0.5 .* d_trans(conn_list(:,1));
d_T_dr = 0.5 .* d_trans(conn_list(:,2));

flux = kh_avg(1:end) .* pot_diff(1:end) .* trans_avg(1:end);
df_dpr = dt .* (kh_avg(1:end) .* ((d_trans(conn_list(:, 2)) .* pot_diff(1:end)) + (trans_avg(1:end) .* d_DP_dr)));
df_dpl = dt .* (kh_avg(1:end) .* ((d_trans(conn_list(:, 1)) .* pot_diff(1:end)) + (trans_avg(1:end) .* d_DP_dl)));

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Residual and Jacobian creation


% resid = accum - accumold - (cell_flux(1:end)) - (dt .* well_info(1:end));
resid = ((0.1781076 .* dx .* dy) .* (accum - accumold)) + (dt .* well_info);
jacob = zeros(ncells, ncells);

counter = 0;
for i = 1: ny
    for k = 1: nx
    counter = counter + 1;
%     counter_horizr = n_vert + ((k - 1) .* (ny - 1)) + i;
%     counter_horizl = n_vert + ((k - 1) .* (ny - 1)) + i - 1;
    jacob(counter, counter) = (0.1781076 .* dx .* dy) .* d_accum(counter);
%     
%     if k > 1
%         jacob(counter, counter - 1) = -df_dpl(((i-1) .* (nx - 1)) + k - 1);        
%     end
%     
%     if k < nx
%         jacob(counter, counter + 1) = -df_dpr(((i-1) .* (nx - 1)) + k); 
%     end
%     
%     if counter <= (ncells - nx)
%         jacob(counter, counter + nx) = -df_dpr(counter_horizr);
%     end
%     
%     if counter > nx
%         jacob(counter, counter - nx) = -df_dpl(counter_horizl);
%     end
%     
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