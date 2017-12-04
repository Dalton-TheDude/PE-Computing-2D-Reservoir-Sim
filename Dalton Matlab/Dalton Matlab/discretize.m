function [ resid, jacob ] = discretize (P, Pold, dt, p0, phi0, b0, cr, cf, visc, kinit, x, y, nx, ny, conn_list)

[accum, d_accum] = accum_calcs(P, p0, phi0, b0, cr, cf);
accumold = accum_calcs(Pold, p0, phi0, b0, cr, cf);
[trans, d_trans] = trans_calcs(P, p0, cf, b0, visc);

dx = x ./ nx;
dy = y ./ ny;
ncells = nx .* ny;

kh_avg = 0.00112712 .* ((2 .* (kinit(conn_list(:, 1)) .* kinit(conn_list(:, 2)))) ./ (kinit(conn_list(:, 1)) + kinit(conn_list(:, 2))));
pot_diff = P(conn_list(:, 2)) - P(conn_list(:, 1));
trans_avg = 0.5 .* ((trans(conn_list(:, 1)) + trans(conn_list(:, 2))));

f = kh_avg(1:end) .* pot_diff(1:end) .* trans_avg(1:end);

df_dpr = kh_avg(1:end) .* ((d_trans(conn_list(:, 2)) .* pot_diff(1:end)) + trans_avg(1:end));
df_dpl = kh_avg(1:end) .* ((d_trans(conn_list(:, 1)) .* pot_diff(1:end)) + trans_avg(1:end));

face_flux = f(conn_list(:, 2)) - f(conn_list(:, 1));
cell_flux = zeros(ncells, 1);

cell_flux(conn_list(:, 2)) = cell_flux(conn_list(:, 2)) - face_flux(conn_list(:, 3));
cell_flux(conn_list(:, 1)) = cell_flux(conn_list(:, 1)) + face_flux(conn_list(:, 3));

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Residual and Jacobian creation
dt_dx(1:ncells, 1) = dt ./ dx;
dt_dy(1:ncells, 1) = dt ./ dy;

n_vert = (nx - 1) .* ny;

resid = accum - accumold - ((dt_dx + dt_dy) .* (cell_flux(1:end))) ;%- cell_flux(1: end-1)));

jacob = zeros(ncells, ncells);

counter = 0;
for i = 1: ny
    for k = 1: nx
    counter = counter + 1;
    counter_horizr = n_vert + ((k - 1) .* (ny - 1)) + i;
    counter_horizl = n_vert + ((k - 1) .* (ny - 1)) + i - 1;
    jacob(counter, counter) = d_accum(counter);
    
    if k > 1
        jacob(counter, counter - 1) = df_dpl(((i-1) .* (nx - 1)) + k - 1);        
    end
    
    if k < nx
        jacob(counter, counter + 1) = df_dpr(((i-1) .* (nx - 1)) + k); 
    end
    
    
    if counter <= (ncells - nx)
        jacob(counter, counter + nx) = df_dpl(counter_horizr);
    end
    
    if counter > nx
        jacob(counter, counter - nx) = df_dpr(counter_horizl);
    end
    
    end
end


end