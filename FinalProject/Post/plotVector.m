function plotVector(n_el_x, n_el_y, disp)

hh_x = 1.0 / n_el_x;
hh_y = 1.0 / n_el_y;

n_np_x = n_el_x + 1;
n_np_y = n_el_y + 1;

[X, Y] = meshgrid(0 : hh_x : 1, 0 : hh_y : 1);

disp_x = reshape(disp(:, 1), n_np_x, n_np_y)';
disp_y = reshape(disp(:, 2), n_np_x, n_np_y)';

quiver(X, Y, disp_x, disp_y, 'AutoScale', 'on', 'AutoScaleFactor', 1.5);

end