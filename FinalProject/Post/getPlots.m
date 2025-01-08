function getPlots(x_coor, y_coor, disp)

figure;
quiver(x_coor, y_coor, disp(:, 1), disp(:, 2), 'AutoScale', 'on', 'AutoScaleFactor', 1.5);

figure;
scatter(x_coor, y_coor, 100, disp(:, 2), 'filled');
colorbar;

end