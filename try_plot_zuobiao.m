% Test program 4 - 参数化版本（r为可调阈值）
initialFun = @(X,Y)(0.3754*((sqrt(X.^2+Y.^2)-0.25)<0));

gf = @(p,X,Y)(p*0+3.1441+(0)*sin(sqrt(X.^2 + Y.^2)));
%gf = @(p,X,Y)(p*0+7.0968+(-5.9086)*sin(sqrt(X.^2 + Y.^2)));

%% ===== 用户可调参数 =====
r = 0.08;  % 定义分界阈值（可自由修改）

%% ===== 通用图形设置 =====
set(0, 'DefaultAxesFontSize', 10);
set(0, 'DefaultTextFontSize', 12);

%% ===== 函数：获取并标记交点 =====
function plot_with_intersections(X, Y, rho, t, r)
    figure('Position', [100, 100, 900, 400]);
    
    % --- 3D 密度图 ---
    subplot(1, 2, 1);
    surf(X, Y, rho, 'EdgeColor', 'none');
    title(['Tumor Density at t = ', num2str(t)], 'FontWeight', 'bold');
    xlabel('X'); ylabel('Y'); zlabel('\rho');
    shading interp;
    colormap(parula);
    colorbar;
    axis tight;
    view(-30, 30);
    grid on;
    set(gca, 'LineWidth', 1.2);
    
    % --- 2D 俯视图 ---
    subplot(1, 2, 2);
    surf(X, Y, rho, 'EdgeColor', 'none');
    view(2);
    title(['2D Interpretation at t = ', num2str(t)], 'FontWeight', 'bold');
    xlabel('X'); ylabel('Y');
    shading interp;
    colormap(parula);
    colorbar;
    axis tight;
    axis equal;
    hold on;
    
    % 获取等高线数据（使用参数r）
    C = contourc(X(1,:), Y(:,1), rho, [r, r]);
    
    % 提取交点（与y=0的交点）
    x_intersections = [];
    idx = 1;
    while idx < size(C,2)
        level = C(1,idx);
        n_points = C(2,idx);
        segment = C(:,idx+1:idx+n_points);
        
        for k = 1:n_points-1
            y1 = segment(2,k);
            y2 = segment(2,k+1);
            if y1*y2 <= 0
                t_interp = -y1/(y2-y1);
                x_inter = segment(1,k) + t_interp*(segment(1,k+1)-segment(1,k));
                x_intersections = [x_intersections, x_inter];
            end
        end
        idx = idx + n_points + 1;
    end
    
    % 绘制等高线（使用参数r）
    contour(X, Y, rho, [r, r], 'r-', 'LineWidth', 2.5);
    
    % 标记交点
    plot(x_intersections, zeros(size(x_intersections)), 'bo', 'MarkerSize', 8, 'LineWidth', 2);
    
    % 更新图例说明（包含r值）
    %legend({['\rho = ', num2str(r)], 'X-axis intersections'}, 'Location', 'northeast', 'Box', 'off');
    
    hold off;
    set(gca, 'LineWidth', 1.2);
    
    % 打印交点坐标
    %fprintf('\nFor t = %.3f (r = %.4f):\n', t, r);
    %fprintf('Intersection x-coordinates with x-axis:\n');
    disp(sort(x_intersections));
end

%% ===== 主程序 =====
% t = 1 的解（传入参数r）
[X1, Y1, rho1, ~] = solver_2D('T', 1, 'GrowthFun', gf, 'dt', 0.005/2, 'InitialFun', initialFun, 'm', 3);
plot_with_intersections(X1, Y1, rho1, 1, r);

% t = 0.875 的解（传入参数r）
[X2, Y2, rho2, ~] = solver_2D('T', 0.875, 'GrowthFun', gf, 'dt', 0.005/2);
plot_with_intersections(X2, Y2, rho2, 0.875, r);