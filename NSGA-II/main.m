% 注意，这个代码要运行比较久，大约56个小时。
% Note that this code takes a long time to run, approximately 56 hours.

% 参数设置
pop_size = 50; % 种群大小
num_generations = 200; % 迭代次数
num_variables = 2; % 决策变量数量
crossover_prob = 1; % 交叉概率
mutation_prob = 0.4; % 变异概率
eta_c = 2; % 交叉分布指数（SBX参数）
eta_m = 5; % 变异分布指数（多项式变异参数）

% 变量范围设置
variable_ranges = [0.5, 3; 2, 20]; % 第一变量范围[0, 3], 第二变量范围[2, 20]

% 初始化种群
disp("initialize the population");
population = initialize_population(pop_size, num_variables, variable_ranges);

% 主循环
for gen = 1:num_generations
    % 选择、交叉和变异生成子代
    offspring_population = generate_offspring(population, crossover_prob, mutation_prob, num_variables, eta_c, eta_m, variable_ranges);

    % 合并父代和子代
    combined_population = [population; offspring_population];

    % 非支配排序和拥挤距离计算
    [sorted_population, ranks] = non_dominated_sorting(combined_population);
    sorted_population = calculate_crowding_distance(sorted_population, ranks);

    % 选择下一代种群
    population = select_next_generation(sorted_population, pop_size);

    % 显示当前代数和种群
    fprintf('-------------------------Generation %d\n------------------', gen);
    plot_population(population);
    drawnow;
end


function population = initialize_population(pop_size, num_variables, variable_ranges)
    % 初始化种群结构体数组
    empty_individual = struct('variables', [], 'objectives', [], 'rank', 0, 'distance', 0, 'domination_count', 0, 'dominated_solutions', []);
    population = repmat(empty_individual, pop_size, 1);
    for i = 1:pop_size
        variables = arrayfun(@(j) variable_ranges(j, 1) + (variable_ranges(j, 2) - variable_ranges(j, 1)) * rand, 1:num_variables);
        objectives = evaluate_objectives(variables);
        population(i).variables = variables;
        population(i).objectives = objectives;
        disp(i);
    end
end

function objectives = evaluate_objectives(variables)
    [f1, f2] = fun_name(variables(1), variables(2));
    objectives = [f1, f2];
end


function offspring_population = generate_offspring(population, crossover_prob, mutation_prob, num_variables, eta_c, eta_m, variable_ranges)
    pop_size = numel(population);
    % 初始化子代种群结构体数组
    empty_individual = struct('variables', [], 'objectives', [], 'rank', 0, 'distance', 0, 'domination_count', 0, 'dominated_solutions', []);
    offspring_population = repmat(empty_individual, pop_size, 1);
    for i = 1:pop_size/2
        % 选择父代
        parent1 = tournament_selection(population);
        parent2 = tournament_selection(population);

        % 交叉
        [child1, child2] = sbx_crossover(parent1, parent2, crossover_prob, eta_c, variable_ranges);

        % 变异
        child1 = polynomial_mutation(child1, mutation_prob, eta_m, num_variables, variable_ranges);
        child2 = polynomial_mutation(child2, mutation_prob, eta_m, num_variables, variable_ranges);

        % 评价子代目标函数
        child1.objectives = evaluate_objectives(child1.variables);
        child2.objectives = evaluate_objectives(child2.variables);

        % 添加到子代种群
        offspring_population(2*i-1) = child1;
        offspring_population(2*i) = child2;

        disp(i);
    end
end

function selected = tournament_selection(population)
    idx1 = randi(numel(population));
    idx2 = randi(numel(population));
    if population(idx1).rank < population(idx2).rank || ...
       (population(idx1).rank == population(idx2).rank && population(idx1).distance > population(idx2).distance)
        selected = population(idx1);
    else
        selected = population(idx2);
    end
end

function [child1, child2] = sbx_crossover(parent1, parent2, prob, eta_c, variable_ranges)
    if rand < prob
        child1 = struct('variables', [], 'objectives', [], 'rank', 0, 'distance', 0, 'domination_count', 0, 'dominated_solutions', []);
        child2 = struct('variables', [], 'objectives', [], 'rank', 0, 'distance', 0, 'domination_count', 0, 'dominated_solutions', []);
        for i = 1:length(parent1.variables)
            if rand <= 0.5
                if abs(parent1.variables(i) - parent2.variables(i)) > 1e-14
                    x1 = min(parent1.variables(i), parent2.variables(i));
                    x2 = max(parent1.variables(i), parent2.variables(i));
                    beta = 1.0 + (2.0 * (x1 - variable_ranges(i, 1)) / (x2 - x1));
                    alpha = 2.0 - beta^-(eta_c + 1.0);
                    rand_val = rand;
                    if rand_val <= 1.0 / alpha
                        beta_q = (rand_val * alpha)^(1.0 / (eta_c + 1.0));
                    else
                        beta_q = (1.0 / (2.0 - rand_val * alpha))^(1.0 / (eta_c + 1.0));
                    end
                    child1.variables(i) = 0.5 * ((x1 + x2) - beta_q * (x2 - x1));
                    child2.variables(i) = 0.5 * ((x1 + x2) + beta_q * (x2 - x1));
                else
                    child1.variables(i) = parent1.variables(i);
                    child2.variables(i) = parent2.variables(i);
                end
            else
                child1.variables(i) = parent1.variables(i);
                child2.variables(i) = parent2.variables(i);
            end
            % 确保变量在范围内
            child1.variables(i) = min(max(child1.variables(i), variable_ranges(i, 1)), variable_ranges(i, 2));
            child2.variables(i) = min(max(child2.variables(i), variable_ranges(i, 1)), variable_ranges(i, 2));
        end
    else
        child1 = parent1;
        child2 = parent2;
    end
end

function mutated = polynomial_mutation(individual, prob, eta_m, num_variables, variable_ranges)
    for i = 1:num_variables
        if rand < prob
            u = rand;
            if u < 0.5
                delta = (2 * u)^(1 / (eta_m + 1)) - 1;
            else
                delta = 1 - (2 * (1 - u))^(1 / (eta_m + 1));
            end
            individual.variables(i) = individual.variables(i) + delta * (variable_ranges(i, 2) - variable_ranges(i, 1));
            % 确保变量在范围内
            individual.variables(i) = min(max(individual.variables(i), variable_ranges(i, 1)), variable_ranges(i, 2));
        end
    end
    mutated = individual;
end

function [sorted_population, ranks] = non_dominated_sorting(population)
    pop_size = numel(population);
    ranks = zeros(pop_size, 1);
    fronts = cell(1, pop_size);
    num_fronts = 0;

    for i = 1:pop_size
        population(i).dominated_solutions = [];
        population(i).domination_count = 0;
        for j = 1:pop_size
            if dominates(population(i), population(j))
                population(i).dominated_solutions = [population(i).dominated_solutions, j];
            elseif dominates(population(j), population(i))
                population(i).domination_count = population(i).domination_count + 1;
            end
        end
        if population(i).domination_count == 0
            population(i).rank = 1;
            ranks(i) = 1;
            fronts{1} = [fronts{1}, i];
        end
    end

    while ~isempty(fronts{num_fronts + 1})
        num_fronts = num_fronts + 1;
        next_front = [];
        for i = fronts{num_fronts}
            for j = population(i).dominated_solutions
                population(j).domination_count = population(j).domination_count - 1;
                if population(j).domination_count == 0
                    population(j).rank = num_fronts + 1;
                    ranks(j) = num_fronts + 1;
                    next_front = [next_front, j];
                end
            end
        end
        fronts{num_fronts + 1} = next_front;
    end

    sorted_population = population;
end

function flag = dominates(individual1, individual2)
    if all(individual1.objectives <= individual2.objectives) && any(individual1.objectives < individual2.objectives)
        flag = true;
    else
        flag = false;
    end
end

function population = calculate_crowding_distance(population, ranks)
    num_objectives = numel(population(1).objectives);
    unique_ranks = unique(ranks);
    for rank = 1:length(unique_ranks)
        front = find(ranks == unique_ranks(rank));
        num_solutions = numel(front);
        if num_solutions < 3
            for i = 1:num_solutions
                population(front(i)).distance = Inf;
            end
            continue;
        end
        objectives = reshape([population(front).objectives], [], num_objectives);
        for m = 1:num_objectives
            [~, sorted_idx] = sort(objectives(:, m));
            population(front(sorted_idx(1))).distance = Inf;
            population(front(sorted_idx(end))).distance = Inf;
            for i = 2:(num_solutions - 1)
                population(front(sorted_idx(i))).distance = population(front(sorted_idx(i))).distance + ...
                    (objectives(sorted_idx(i+1), m) - objectives(sorted_idx(i-1), m));
            end
        end
    end
end

function next_population = select_next_generation(population, pop_size)
    [~, sort_idx] = sortrows([[population.rank]', -[population.distance]']);
    next_population = population(sort_idx(1:pop_size));
end

function plot_population(population)
    objectives = reshape([population.objectives], 2, []).';
    scatter(objectives(:, 1), objectives(:, 2), 'filled');
    xlabel('$realTimeFun$', Fontname='Times New Roman', Interpreter='latex');
    ylabel('$safeFun$', Fontname='Times New Roman', Interpreter='latex');
    title('NSGA-II');
end