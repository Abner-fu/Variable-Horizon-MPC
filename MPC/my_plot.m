clc;
clear;
close all;

load data_traj.mat;

obstacle_sq_1.duration  = [2 3];
obstacle_sq_1.width     = [4 3];
obstacle_sq_2.duration  = [5.3 5.8];
obstacle_sq_2.width     = [6.5 5.5];

obstacle1.center = [5 4];
obstacle1.radius = 0.5;

figure('color','w')
fig = figure(1);
hold on;
plot(xx_2(1,:), xx_2(2,:), LineWidth=1);
plot(xx_10(1,:), xx_10(2,:), LineWidth=1);
plot(xx_10_barrier(1,:), xx_10_barrier(2,:), LineWidth=1);
plot(xx_20(1,:), xx_20(2,:), LineWidth=1);
plot(xx_var(1,:), xx_var(2,:), LineWidth=1);
plot(1, 1, 'o');
plot(7, 7, 'o');
legend("$N=2$", "$N=10$", "$N=10(barrier)$", "$N=20$", "$N$ vary(ours)", "starting point", "terminal point", Location="northwest", Fontname='Times New Roman', Interpreter='latex')
xlabel("$x$(m)", Fontname='Times New Roman', Interpreter='latex')
ylabel("$y$(m)", Fontname='Times New Roman', Interpreter='latex')
rectangle('Position',[obstacle1.center(1)-obstacle1.radius, obstacle1.center(2)-obstacle1.radius, 2*obstacle1.radius, 2*obstacle1.radius], 'Curvature',[1,1]);
rectangle('Position', [min(obstacle_sq_1.duration) min(obstacle_sq_1.width) abs(obstacle_sq_1.duration(2)-obstacle_sq_1.duration(1)) abs(obstacle_sq_1.width(2)-obstacle_sq_1.width(1))], 'Curvature', 0.3)
rectangle('Position', [min(obstacle_sq_2.duration) min(obstacle_sq_2.width) abs(obstacle_sq_2.duration(2)-obstacle_sq_2.duration(1)) abs(obstacle_sq_2.width(2)-obstacle_sq_2.width(1))], 'Curvature', 0.3)

xlim([0, 8])
ylim([0, 8])
axis equal
grid on;
legend

idx_1 = find((xx_2(1,:) > 2.5) & (xx_2(1, :) < 3.5));
x_temp_1 = xx_2(1, idx_1);
y_temp_1 = xx_2(2, idx_1);

idx_10 = find((xx_10(1,:) > 2.5) & (xx_10(1, :) < 3.5));
x_temp_10 = xx_10(1, idx_10);
y_temp_10 = xx_10(2, idx_10);

idx_10_barrier = find((xx_10_barrier(1,:) > 2.5) & (xx_10_barrier(1, :) < 3.5));
x_temp_10_barrier = xx_10_barrier(1, idx_10_barrier);
y_temp_10_barrier = xx_10_barrier(2, idx_10_barrier);

idx_20 = find((xx_20(1,:) > 2.5) & (xx_20(1, :) < 3.5));
x_temp_20 = xx_20(1, idx_20);
y_temp_20 = xx_20(2, idx_20);

idx_var = find((xx_var(1,:) > 2.5) & (xx_var(1, :) < 3.5));
x_temp_var = xx_var(1, idx_var);
y_temp_var = xx_var(2, idx_var);

ax2 = axes();
ax2.Position=[0.63392857142857,0.166666666666667,0.251785714285714,0.247619047619048];
hold on;
plot(x_temp_1, y_temp_1, "Parent",ax2, 'LineWidth', 1);
plot(x_temp_10, y_temp_10, "Parent",ax2, 'LineWidth', 1);
plot(x_temp_10_barrier, y_temp_10_barrier, "Parent",ax2, 'LineWidth', 1);
plot(x_temp_20, y_temp_20, "Parent",ax2, 'LineWidth', 1);
plot(x_temp_var, y_temp_var, "Parent",ax2, 'LineWidth', 1);
grid on;

annotation('arrow', [0.4536,0.5589] ,[0.3561,0.2881]);

% export_fig 'traj' -eps