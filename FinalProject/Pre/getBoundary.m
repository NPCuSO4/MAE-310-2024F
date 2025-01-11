function node_pos = getBoundary(n_np, x_coor, y_coor)

node_pos = zeros(n_np, 1); % 0: interior / 1: left / 2: right / 3: lower / 4: upper / 5: arc / 6: top-right corner point

for ii = 1 : n_np
    % Adjust this part for different mesh
    
    %if ((x_coor(ii)+2)^2 + (y_coor(ii)+2)^2)^0.5 < 1e-7
    %    node_pos(ii) = -1;  
    %elseif ((x_coor(ii)-2)^2 + (y_coor(ii)-2)^2)^0.5 < 1e-7
    %    node_pos(ii) = 6;
    if abs(x_coor(ii) + 1) < 1e-7
        node_pos(ii) = 1;
    elseif abs(x_coor(ii) - 1) < 1e-7
        node_pos(ii) = 2;    
    elseif abs(y_coor(ii) + 1) < 1e-7
        node_pos(ii) = 3;
    elseif abs(y_coor(ii) - 1) < 1e-7
        node_pos(ii) = 4;       
    %elseif ((x_coor(ii)+2)^2 + (y_coor(ii)+2)^2)^0.5 < 0.5+1e-7
    %    node_pos(ii) = 5;
    end

end

end