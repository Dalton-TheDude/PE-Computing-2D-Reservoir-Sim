function connection_list = connection_list(nx, ny)

nfaces = ((nx-1) .* ny) + ((ny-1) .* nx);
connection_list = zeros(nfaces, 3); % (Left/lower cell, right/upper cell, face #)

count = 1;
for i = 1: ny
    for j = 1: nx-1
        connection_list(count, 1) = ((i-1) .* nx) + j;
        connection_list(count, 2) = ((i-1) .* nx) + j + 1;
        connection_list(count, 3) = count;
        count = count + 1;
    end
end

for i = 1: ny-1
    for j = 1: nx
        connection_list(count, 1) = (((i-1) .* nx) + j);
        connection_list(count, 2) = (((i-1) .* nx) + (j+nx));
        connection_list(count, 3) = count;
        count = count + 1;
    end
end
