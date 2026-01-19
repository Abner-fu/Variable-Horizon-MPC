function aver_Distance = compute_traj_dist(points)
    % 输入：
    % points - 一个二维数组，每一行是一个点，第一列是x坐标，第二列是y坐标
    % 输出：
    % totalDistance - 所有点到障碍物的距离之和

    % 初始化总距离为0
    totalDistance = 0;

    obstacle_sq_1_list = [3; 4; 4; 5];
    % 遍历每一个点
    for i = 1:size(points, 1)
        % 获取当前点的坐标
        point = points(i, :);
        
        % 计算当前点到障碍物的距离
        distance = pointToRectangleDistance(point(1), point(2), obstacle_sq_1_list(1),obstacle_sq_1_list(2), obstacle_sq_1_list(3), obstacle_sq_1_list(4));
        
        % 将距离累加到总距离中
        totalDistance = totalDistance + exp(-(8*distance-2));
    end
    aver_Distance = totalDistance;
end
