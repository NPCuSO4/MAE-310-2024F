function node_pos = getBoundary(n_np, x_coor, y_coor)

node_pos = zeros(n_np, 1); % 0: interior / 1: left / 2: right / 3: lower / 4: upper

for ii = 1 : n_np
    if abs(x_coor(ii) + 1) < 1e-5
        node_pos(ii) = 1;
    elseif abs(x_coor(ii) - 1) < 1e-5
        node_pos(ii) = 2;
    elseif abs(y_coor(ii) + 1) < 1e-5
        node_pos(ii) = 3;
    elseif abs(y_coor(ii) - 1) < 1e-5
        node_pos(ii) = 4;
    end
end

end