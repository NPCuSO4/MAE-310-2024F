function getPlots(sam_x, sam_y, sam_u, sam_strain, sam_stress)

figure;
quiver(sam_x, sam_y, sam_u(:, 1), sam_u(:, 2), 'AutoScale', 'on', 'AutoScaleFactor', 1.5);

figure;
scatter(sam_x, sam_y, 100, sam_stress(:, 1), 's', 'filled');
colorbar;

figure;
scatter(sam_x, sam_y, 100, sam_stress(:, 2), 's', 'filled');
colorbar;

end