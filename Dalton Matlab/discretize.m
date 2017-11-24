function [ resid, jacob ] = discretize (P, Pold, dt, p0, phi0, b0, cr, cf, visc, kinit, x, y, nx, ny, nfaces)

[accum, d_accum] = accum_calcs(P, p0, phi0, b0, cr, cf);
accumold = accum_calcs(Pold, p0, phi0, b0, cr, cf);
[trans, d_trans] = trans_calcs(P, p0, cf, b0, visc);

dx = x ./ nx;
dy = y ./ ny;
face_list = zeros(nfaces, 1);
left_face = face_list(2: end);
right_face = face_list(1: end - 1);

ncells = nx .* ny;
f = zeros((nfaces) + 2, 1);
kh_avg = zeros(nfaces, 1);
pot_diff = zeros(nfaces, 1);
trans_avg = zeros(nfaces, 1);

for i = 1: ny
    for j = 2: nx
        kh_avg(((i-1) .* (ny-1)) + (j-1)) = (2 .* kinit(((i - 1) .* nx) + j) .* kinit(((i - 1) .* nx) + (j - 1))) ./ (kinit(((i - 1) .* nx) + j) + kinit(((i - 1) .* nx) + (j - 1)));
        pot_diff(((i-1) .* (ny-1)) + (j-1)) = P(((i - 1) .* nx) + j) - P(((i - 1) .* nx) + j - 1);
        trans_avg(((i-1) .* (ny-1)) + (j-1)) = 0.5 .* (trans(((i - 1) .* nx) + j) + trans(((i - 1) .* nx) + j-1));
    end
end

for i = 2: ny
    for j = 1: nx
        kh_avg(((nx - 1) .* ny) + ((nx-1) .* (j-1)) + (i-1)) = (2 .* kinit(((i - 1) .* nx) + j) .* kinit(((i - 2) .* nx) + (j))) ./ (kinit(((i - 1) .* nx) + j) + kinit(((i - 2) .* nx) + (j)));
        pot_diff(((nx - 1) .* ny) + ((nx-1) .* (j-1)) + (i-1)) = P(((i - 1) .* nx) + j) - P(((i - 2) .* nx) + (j));
        trans_avg(((nx - 1) .* ny) + ((nx-1) .* (j-1)) + (i-1)) = 0.5 .* (trans(((i - 1) .* nx) + j) + trans(((i - 2) .* nx) + (j)));
    end
end


f(2: end - 1) = kh_avg(1:end) .* pot_diff(1:end) .* trans_avg(1:end);

resid = accum - accumold - ((dt ./ dx) .* (f(2:end) - f(1: end-1)));

jaccum = spdiags(d_accum, 0, nx .* ny, nx .* ny);

end