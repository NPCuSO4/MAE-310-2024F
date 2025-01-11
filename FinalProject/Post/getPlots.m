function getPlots(sam_x, sam_y, sam_u, sam_strain, sam_stress)

% Displacement vector
figure;
quiver(sam_x, sam_y, sam_u(:, 1), sam_u(:, 2), 'AutoScale', 'on', 'AutoScaleFactor', 1.5);

% Stress vector
figure;
quiver(sam_x, sam_y, sam_stress(:, 1), sam_stress(:, 2), 'AutoScale', 'on', 'AutoScaleFactor', 1.5);

xlin = linspace(min(sam_x), max(sam_x), 100);
ylin = linspace(min(sam_y), max(sam_y), 100);
[X, Y] = meshgrid(xlin, ylin);

% Stress sigma_xx
Z = griddata(sam_x, sam_y, sam_stress(:, 1), X, Y, 'linear'); 
figure;
contourf(X, Y, Z, 20, 'LineColor', 'none');
colorbar;
xlabel('X-axis');
ylabel('Y-axis');
title('Contour Plot of \sigma_x');


end