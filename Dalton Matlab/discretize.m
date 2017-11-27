function [ resid, jacob ] = discretize (P, Pold, dt, p0, phi0, b0, cr, cf, visc, kinit, x, y, nx, ny, nfaces)

[accum, d_accum] = accum_calcs(P, p0, phi0, b0, cr, cf);
accumold = accum_calcs(Pold, p0, phi0, b0, cr, cf);
[trans, d_trans] = trans_calcs(P, p0, cf, b0, visc);

dx = x ./ nx;
dy = y ./ ny;

f = zeros((nfaces), 1);
kh_avg = zeros(nfaces, 1);
pot_diff = zeros(nfaces, 1);
trans_avg = zeros(nfaces, 1);
df_dpr = zeros(nfaces, 1);
df_dpl = zeros(nfaces, 1);

for i = 1: ny
    for j = 2: nx
        counter_facev = ((i-1) .* (ny-1)) + (j-1);
        counter_cellr = ((i - 1) .* nx) + j;
        counter_celll = ((i - 1) .* nx) + (j - 1);
        kh_avg(counter_facev) = (2 .* kinit(counter_cellr) .* kinit(counter_celll) ./ (kinit(counter_cellr) + kinit(counter_celll)));
        pot_diff(counter_facev) = P(counter_cellr) - P(counter_celll);
        trans_avg(counter_facev) = 0.5 .* (trans(counter_cellr) + trans(counter_celll));
        df_dpr(counter_facev) = kh_avg(counter_facev) .* (((d_trans(counter_cellr) .* pot_diff(counter_facev))) + trans_avg(counter_facev));
        df_dpl(counter_facev) = kh_avg(counter_facev) .* (((d_trans(counter_celll) .* pot_diff(counter_facev))) + trans_avg(counter_facev));
    end
end

for i = 2: ny
    for j = 1: nx
        counter_faceh = ((nx - 1) .* ny) + ((nx-1) .* (j-1)) + (i-1);
        counter_cellr = ((i - 1) .* nx) + j;
        counter_celll = ((i - 2) .* nx) + j;
        kh_avg(counter_faceh) = (2 .* kinit(counter_cellr) .* kinit(counter_celll) ./ (kinit(counter_cellr) + kinit(counter_celll)));
        pot_diff(counter_faceh) = P(counter_cellr) - P(counter_celll);
        trans_avg(counter_faceh) = 0.5 .* (trans(counter_cellr) + trans(counter_celll));
        df_dpr(counter_faceh) = kh_avg(counter_faceh) .* (((d_trans(counter_cellr) .* pot_diff(counter_faceh))) + trans_avg(counter_faceh));
        df_dpl(counter_faceh) = kh_avg(counter_faceh) .* (((d_trans(counter_celll) .* pot_diff(counter_faceh))) + trans_avg(counter_faceh));
    end
end

f(1:end) = kh_avg(1:end) .* pot_diff(1:end) .* trans_avg(1:end); 
        
vertical_flux = zeros(nx .* ny, 1);
horizontal_flux = zeros(nx .* ny, 1);

for i = 1: ny
    for j = 1: nx
        counter = ((i-1) .* nx) + j;
        v_counterl = ((i-1) .* (nx-1))+(j-1);
        v_counterr = ((i-1) .* (nx-1))+(j);
        if j == 1
            vertical_flux(counter) = f(v_counterr);
        else
            if j == nx
                vertical_flux(counter) = f(v_counterl)+(j-1);
            else
                vertical_flux(counter) = f(v_counterl)+f(v_counterr);
            end
        end
        
        h_counterl = ((nx-1) .* ny) + ((j-1) .* (ny-1)) + (i-1);
        h_counterr = ((nx-1) .* ny) + ((j-1) .* (ny-1)) + i;
        if i == 1
            horizontal_flux(counter) = f(h_counterr);
        else
            if i == ny
                horizontal_flux(counter) = f(h_counterl);
            else
                horizontal_flux(counter) = f(h_counterr) + f(h_counterl);          
            end
        end
        
    end
end

total_flux = vertical_flux(1:end) + horizontal_flux(1:end);

resid = accum - accumold - ((dt ./ dx) .* (total_flux(2:end) - total_flux(1: end-1)));

jacob = zeros(nx .* ny, nx .* ny + 4);

n_vert = (nx - 1) .* ny;

for i = 1: ny
    for k = 1: nx
    counter = ((i-1) .* nx) + k;
    counter_horizr = n_vert + ((k - 1) .* (ny - 1)) + i;
    counter_horizl = n_vert + ((k - 1) .* (ny - 1)) + i - 1;
    jacob(counter, counter +2) = jacob(counter, counter +2) + d_accum(counter);
    
    if k == 1
        jacob(counter, counter +2- 1) = 0;
    else
        jacob(counter, counter +2- 1) = jacob(counter, counter +2- 1) + df_dpl(((i-1) .* (nx - 1)) + k - 1);
    end
    
    if i == 1
        jacob(counter, counter +2- 2) = 0;
    else
        jacob(counter, counter +2- 2) = jacob(counter, counter +2- 2) + df_dpl(counter_horizl);
    end
    
    if k == nx
        jacob(counter, counter +2+ 1) = 0;
    else
        jacob(counter, counter +2+ 1) = jacob(counter, counter +2+ 1) + df_dpr(((i-1) .* (nx - 1)) + k);
    end
    
    if i == ny
        jacob(counter, counter +2+ 2) = 0;
    else
        jacob(counter, counter +2+ 2) = jacob(counter, counter +2+ 2) + df_dpr(counter_horizr);
    end
    
    
    end
end


end