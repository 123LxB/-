% �����Ŵ��㷨��դ�񷨻�����·���滮��MATLABʵ��
clc;
clear;
% ��������,��դ���ͼ
G=  [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
     0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
     0 0 1 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0;
     0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
     0 1 1 1 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 1 1 1 0 1 0 1 1 0 0 0 0 0 0;
     0 1 1 1 0 0 0 0 0 0 1 0 1 1 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 1 1 0 1 1 1 1 0 0 0 0 0 0;
     0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 1 1 1 1 0;
     0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0;
     0 0 1 1 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 0;
     0 0 1 1 0 0 1 1 1 0 1 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 1 1 1 0 1 0 0 0 0 0 0 1 1 0; 
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 1 0;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
p_start = 0;        % ��ʼ���
p_end = 399;        % ��ֹ���
NP = 200;           % ��Ⱥ����
max_gen = 50;       % ����������
percross = 0.8;     % �������
permutate = 0.2;    % �������
perlenth = 1;       % ·�����ȱ���
persmooth = 7;      % ·��˳���ȱ���
%init_path = [];
z = 1; 
new_pop1 = {}; % Ԫ���������飬ɶ�����Դ棬����C�����еĽṹ
new_pop2 = {};
[y, x] = size(G); % x = 20, y = 20������������
% ��������У������ұ��1.2.3...��
xstart = mod(p_start, x) + 1; 
% ��������У����ϵ��±����1.2.3...��
ystart = fix(p_start / x) + 1;
% �յ������С���
xend = mod(p_end, x) + 1; % mod����ȡ�ຯ����ǰ���Ժ�ȡ����
yend = fix(p_end / x) + 1; % fix����ȡ����������������ȡ����

% ��Ⱥ��ʼ����һ�����ؾ��ڵ�,����ʼ�������п�ʼ���ϣ���ÿ������ѡһ������դ�񣬹��ɱؾ��ڵ�
pass_num = yend - ystart + 1; % ���ڵ�ͼ�ܹ�����������ʾ
pop = zeros(NP, pass_num); % pop��ʼ�����㣬��������������*��ͼ�������ľ��󣬱�ʾ���п��ܵ�·��
for global1_i = 1 : NP
    pop(global1_i, 1) = p_start; % ѭ��NP�Σ���������pop����
    global1_j = 1;
    for global_iyk = ystart+1 : yend-1 % ��ȥ�����յ�
        global1_j = global1_j + 1;
        can = []; % ÿһ�����ϰ��ĵ����can����
        for global_ixk = 1 : x
            global_ino = (global_ixk - 1) + (global_iyk - 1) * x; % ��ÿ��դ������
            if G(global_iyk, global_ixk) == 0
                can = [can global_ino]; % �����ϰ��ĵ����can������
            end
        end
        can_num = length(can);  % can����ĳ���
        index = randi(can_num); % �������������
        pop(global1_i, global1_j) = can(index); % ÿһ�б�ʾ��һ��·���м���һ�����ϰ��ĵ�
    end
    pop(global1_i, end) = p_end; % �����յ����pop���һ��
    
    % ��Ⱥ��ʼ���ڶ������������ؾ��ڵ�������޼��·��
    % single_new_pop = generate_continuous_path(pop(i, :), G, x);
    %@@@@@@@@@@ ����������ĵ����ӳ�����·��
    generate_i = 1;
    single_new_pop = pop(global1_i,:); % ȡpop�����global1_i�ж����ǵ�generate_i���������ݣ���ʾȡ��һ��·��
    [~, single_path_num] = size(single_new_pop);
    while generate_i ~= single_path_num
        value_x_now = mod(single_new_pop(1, generate_i), x) + 1;  % ��ǰ��������
        value_y_now = fix(single_new_pop(1, generate_i) / x) + 1; % ��ǰ��������
        value_x_next = mod(single_new_pop(1, generate_i + 1), x) + 1;  % ��һ����������
        value_y_next = fix(single_new_pop(1, generate_i + 1) / x) + 1; % ��һ����������
        max_iteration = 0; % ��ʼ��ƽ��·������������
        while max(abs(value_x_next - value_x_now), abs(value_y_next - value_y_now)) ~= 1 % �жϵ�generate_i��generate_i+1�Ƿ�����������һ��,����Ͳ����м��
            x_insert = floor((value_x_next + value_x_now) / 2);
            y_insert = floor((value_y_next + value_y_now) / 2);
            if G(y_insert, x_insert) == 0 % ������դ��Ϊ����դ�񣨲����ϰ����ϣ�ʱ
                num_insert = (x_insert - 1) + (y_insert - 1) * x; % ��դ��λ��ת��Ϊ���
                single_new_pop = [single_new_pop(1, 1:generate_i), num_insert, single_new_pop(1, generate_i+1:end)]; % ��դ����Ų���single������        
            else % ����դ��Ϊ�ϰ���դ��(���ϰ����ϣ�ʱ
                % ������
                if G(y_insert, x_insert - 1) == 0 && ((x_insert - 2) + (y_insert - 1) * x ~= single_new_pop(1, generate_i)) && ((x_insert - 2) + (y_insert - 1) * x ~= single_new_pop(1, generate_i+1))
                    x_insert = x_insert - 1;
                    num_insert = (x_insert - 1) + (y_insert - 1) * x; % ��դ��λ��ת��Ϊ���
                    single_new_pop = [single_new_pop(1, 1:generate_i), num_insert, single_new_pop(1, generate_i+1:end)]; % ��դ����Ų���single������               
                % ������    
                elseif G(y_insert, x_insert + 1) == 0 && (x_insert + (y_insert - 1) * x ~= single_new_pop(1, generate_i)) && (x_insert + (y_insert - 1) * x ~= single_new_pop(1, generate_i+1))
                    x_insert = x_insert + 1;
                    num_insert = (x_insert - 1) + (y_insert - 1) * x; % ��դ��λ��ת��Ϊ���
                    single_new_pop = [single_new_pop(1, 1:generate_i), num_insert, single_new_pop(1, generate_i+1:end)]; % ��դ����Ų���single������                
                % ������
                elseif G(y_insert + 1, x_insert) == 0 && ((x_insert - 1) + y_insert * x ~= single_new_pop(1, generate_i)) && ((x_insert - 1) + y_insert * x ~= single_new_pop(1, generate_i+1))
                    y_insert = y_insert + 1;
                    num_insert = (x_insert - 1) + (y_insert - 1) * x; % ��դ��λ��ת��Ϊ���
                    single_new_pop = [single_new_pop(1, 1:generate_i), num_insert, single_new_pop(1, generate_i+1:end)]; % ��դ����Ų���single������
                % ������
                elseif  G(y_insert - 1, x_insert) == 0 && ((x_insert - 1) + (y_insert - 2) * x ~= single_new_pop(1, generate_i)) && ((x_insert - 1) + (y_insert-2) * x ~= single_new_pop(1, generate_i+1))
                    y_insert = y_insert - 1;
                    num_insert = (x_insert - 1) + (y_insert - 1) * x; % ��դ��λ��ת��Ϊ���
                    single_new_pop = [single_new_pop(1, 1:generate_i), num_insert, single_new_pop(1, generate_i+1:end)]; % ��դ����Ų���single������
                % ���������ȥ��·��
                else                
                    single_new_pop = [];
                    break
                end    
            end        
            value_x_next = x_insert;
            value_y_next = y_insert;
            max_iteration = max_iteration + 1;
            if max_iteration > 20000 % ���ִ��20000�Σ���Ϊƽ��·������������
                single_new_pop = [];
                break
            end    
        end    
        if isempty(single_new_pop)
            break
        end
        [~, single_path_num] = size(single_new_pop);
        generate_i = generate_i + 1;
    end
    if ~isempty(single_new_pop)
       new_pop1(z, 1) = {single_new_pop};
        z = z + 1;
    end
end

% �����ʼ����Ⱥ����Ӧ��
%@@@@@@@@@@ ����·������
% path_value = cal_path_value(new_pop1, x);
path_value_in_pop = new_pop1;
[value_n, ~] = size(path_value_in_pop); 
path_value_out_pop = zeros(1, value_n);
for value_i = 1 : value_n
    path_value_single_pop = path_value_in_pop{value_i, 1};
    [~, value_m] = size(path_value_single_pop);
    for value_j = 1 : value_m - 1
        % ��i�����У������ұ��1.2.3...��
        value_x_now = mod(path_value_single_pop(1, value_j), x) + 1; 
        % ��i�����У����ϵ��±����1.2.3...��
        value_y_now = fix(path_value_single_pop(1, value_j) / x) + 1;
        % ��i+1�����С���
        value_x_next = mod(path_value_single_pop(1, value_j + 1), x) + 1;
        value_y_next = fix(path_value_single_pop(1, value_j + 1) / x) + 1;
        if abs(value_x_now - value_x_next) + abs(value_y_now - value_y_next) == 1
            path_value_out_pop(1, value_i) = path_value_out_pop(1, value_i) + 1;
        else
            path_value_out_pop(1, value_i) = path_value_out_pop(1, value_i) + sqrt(2);
        end
    end
end
path_value = path_value_out_pop;

%@@@@@@@@@@ ����·��ƽ����
%path_smooth = cal_path_smooth(new_pop1, x);
path_smooth_in_pop = new_pop1;
[smooth_n, ~] = size(path_smooth_in_pop);
path_smooth_out_pop = zeros(1, smooth_n);
for smooth_i = 1 : smooth_n
    path_smooth_single_pop = path_smooth_in_pop{smooth_i, 1};
    [~, smooth_m] = size(path_smooth_single_pop);
    for smooth_j = 1 : smooth_m - 2
        % ��i�����У������ұ��1.2.3...��
        smooth_x_now = mod(path_smooth_single_pop(1, smooth_j), x) + 1; 
        % ��i�����У����ϵ��±����1.2.3...��
        smooth_y_now = fix(path_smooth_single_pop(1, smooth_j) / x) + 1;
        % ��i+1�����С���
        % smooth_x_next1 = mod(path_smooth_single_pop(1, smooth_j + 1), x) + 1;
        % smooth_y_next1 = fix(path_smooth_single_pop(1, smooth_j + 1) / x) + 1;
        % ��i+2�����С���
        smooth_x_next2 = mod(path_smooth_single_pop(1, smooth_j + 2), x) + 1;
        smooth_y_next2 = fix(path_smooth_single_pop(1, smooth_j + 2) / x) + 1;
        smooth_c2 = (smooth_x_now - smooth_x_next2)^2 + (smooth_y_now - smooth_y_next2)^2;
        if smooth_c2 < 8 && smooth_c2 > 4
            path_smooth_out_pop(1, smooth_i) = path_smooth_out_pop(1, smooth_i) + 5;
        elseif smooth_c2 <= 4 && smooth_c2 > 1
            path_smooth_out_pop(1, smooth_i) = path_smooth_out_pop(1, smooth_i) + 30;
        elseif smooth_c2 <= 1
            path_smooth_out_pop(1, smooth_i) = path_smooth_out_pop(1, smooth_i) + 5000;
        end
    end
end
path_smooth = path_smooth_out_pop;

% ��·��ƽ�����볤�ȼ���·����Ӧ�Ⱥ���
fit_value = perlenth .* path_value .^ -1 + persmooth .* path_smooth .^ -1;
mean_path_value = zeros(1, max_gen);
min_path_value = zeros(1, max_gen);

% ѭ�����������������Ŵ��㷨����Խ�ԣ�
for global2_i = 1 : max_gen
    %@@@@@@@@@@ ѡ�����
    % new_pop2 = selection(new_pop1, fit_value);
    select_in_pop = new_pop1;
    select_out_pop = cell(length(new_pop1),1); % ����һ��Ԫ�����飬ɶ�����Դ�
    [select_px, ~] = size(select_in_pop);
    total_fit = sum(fit_value);
    p_fit_value = fit_value / total_fit;
    p_fit_value = cumsum(p_fit_value);  % �����������
    % �������С��������
    select_ms = sort(rand(select_px, 1));   
    fitin = 1;
    newin = 1;
    while newin <= select_px
        if(select_ms(newin)) < p_fit_value(fitin)
            select_out_pop{newin, 1} = select_in_pop{fitin, 1};
            newin = newin+1;
        else
            fitin = fitin+1;
        end
    end
    new_pop2 = select_out_pop;

    %@@@@@@@@@@ �������
    % new_pop2 = crossover(new_pop2, percross);
    cross_in_pop = new_pop2;
    cross_out_pop = cell(length(new_pop1),1);
    [cross_px,~] = size(cross_in_pop);
    parity = mod(cross_px, 2);% �ж�·��������������ż��
    for cross_i = 1:2:cross_px-1
        singal_now_pop = cross_in_pop{cross_i, 1};
        singal_next_pop = cross_in_pop{cross_i+1, 1};
        [lia, lib] = ismember(singal_now_pop, singal_next_pop);
        [~, cross_n] = find(lia == 1);
        [~, cross_m] = size(cross_n);
        if (rand < percross) && (cross_m >= 3)
            % ����һ��2-m-1֮��������
            cross_r = round(rand(1,1)*(cross_m-3)) +2;
            crossover_index1 = cross_n(1, cross_r);
            crossover_index2 = lib(crossover_index1);
            cross_out_pop{cross_i, 1} = [singal_now_pop(1:crossover_index1), singal_next_pop(crossover_index2+1:end)];
            cross_out_pop{cross_i+1, 1} = [singal_next_pop(1:crossover_index2), singal_now_pop(crossover_index1+1:end)];       
        else
            cross_out_pop{cross_i, 1} =singal_now_pop;
            cross_out_pop{cross_i+1, 1} = singal_next_pop;
        end
        if parity == 1
            cross_out_pop{cross_px, 1} = cross_in_pop{cross_px, 1};
        end
    end
    new_pop2 = cross_out_pop;
    
    %@@@@@@@@@@ �������
    % new_pop2 = mutation(new_pop2, permutate, G, x);
    mutate_in_pop = new_pop2;
    mutate_out_pop = cell(length(new_pop1),1);
    [mutate_px, ~] = size(mutate_in_pop);
    for mutate_i = 1:mutate_px
        % ��ʼ������������
        mutate_max = 0;
        mutate_single_new_pop = mutate_in_pop{mutate_i, 1};
        [~, mutate_m] = size(mutate_single_new_pop);
        % single_new_pop_slice��ʼ��
        mutate_single_new_pop_slice = [];
        if(rand < permutate) % rand���һ��0-1����λС��
            while isempty(mutate_single_new_pop_slice)
                % ����2-��m-1�������������,������
                mutate_mpoint = sort(round(rand(1,2)*(mutate_m-3)) + [2 2]);
                mutate_single_new_pop_slice = [mutate_single_new_pop(mutate_mpoint(1, 1)-1) mutate_single_new_pop(mutate_mpoint(1, 2)+1)];
                %@@@@@@@@@@ �������ؾ��ڵ�������޼��·��
                % mutate_single_new_pop_slice = generate_continuous_path(mutate_single_new_pop_slice, G, x);
                % mutate_generate_in_pop = mutate_single_new_pop_slice;
                mutate_generate_out_pop = mutate_single_new_pop_slice;
                mutate_generate_i = 1;
                [~, mutate_single_path_num] = size(mutate_generate_out_pop);
                while mutate_generate_i ~= mutate_single_path_num
                    mutate_x_now = mod(mutate_generate_out_pop(1, mutate_generate_i), x) + 1;  % ��ǰ��������
                    mutate_y_now = fix(mutate_generate_out_pop(1, mutate_generate_i) / x) + 1; % ��ǰ��������
                    mutate_x_next = mod(mutate_generate_out_pop(1, mutate_generate_i + 1), x) + 1;  % ��һ����������
                    mutate_y_next = fix(mutate_generate_out_pop(1, mutate_generate_i + 1) / x) + 1; % ��һ����������
                    mutate_max_iteration = 0; % ��ʼ��ƽ��·������������
                    while max(abs(mutate_x_next - mutate_x_now), abs(mutate_y_next - mutate_y_now)) ~= 1 % �жϵ�generate_i��generate_i+1�Ƿ�����������һ��,����Ͳ����м��
                        mutate_x_insert = floor((mutate_x_next + mutate_x_now) / 2);
                        mutate_y_insert = floor((mutate_y_next + mutate_y_now) / 2);
                        if G(mutate_y_insert, mutate_x_insert) == 0 % ������դ��Ϊ����դ�񣨲����ϰ����ϣ�ʱ
                            mutate_num_insert = (mutate_x_insert - 1) + (mutate_y_insert - 1) * x; % ��դ��λ��ת��Ϊ���
                            mutate_generate_out_pop = [mutate_generate_out_pop(1, 1:mutate_generate_i), mutate_num_insert, mutate_generate_out_pop(1, mutate_generate_i+1:end)]; % ��դ����Ų���single������        
                        else % ����դ��Ϊ�ϰ���դ��(���ϰ����ϣ�ʱ
                            % ������
                            if G(mutate_y_insert, mutate_x_insert - 1) == 0 && ((mutate_x_insert - 2) + (mutate_y_insert - 1) * x ~= mutate_generate_out_pop(1, mutate_generate_i)) && ((mutate_x_insert - 2) + (mutate_y_insert - 1) * x ~= mutate_generate_out_pop(1, mutate_generate_i+1))
                                mutate_x_insert = mutate_x_insert - 1;
                                mutate_num_insert = (mutate_x_insert - 1) + (mutate_y_insert - 1) * x; % ��դ��λ��ת��Ϊ���
                                mutate_generate_out_pop = [mutate_generate_out_pop(1, 1:mutate_generate_i), mutate_num_insert, mutate_generate_out_pop(1, mutate_generate_i+1:end)]; % ��դ����Ų���single������               
                            % ������    
                            elseif G(mutate_y_insert, mutate_x_insert + 1) == 0 && (mutate_x_insert + (mutate_y_insert - 1) * x ~= mutate_generate_out_pop(1, mutate_generate_i)) && (mutate_x_insert + (mutate_y_insert - 1) * x ~= mutate_generate_out_pop(1, mutate_generate_i+1))
                                mutate_x_insert = mutate_x_insert + 1;
                                mutate_num_insert = (mutate_x_insert - 1) + (mutate_y_insert - 1) * x; % ��դ��λ��ת��Ϊ���
                                mutate_generate_out_pop = [mutate_generate_out_pop(1, 1:mutate_generate_i), mutate_num_insert, mutate_generate_out_pop(1, mutate_generate_i+1:end)]; % ��դ����Ų���single������                
                            % ������
                            elseif G(mutate_y_insert + 1, mutate_x_insert) == 0 && ((mutate_x_insert - 1) + mutate_y_insert * x ~= mutate_single_new_pop(1, mutate_generate_i)) && ((mutate_x_insert - 1) + mutate_y_insert * x ~= mutate_generate_out_pop(1, mutate_generate_i+1))
                                mutate_y_insert = mutate_y_insert + 1;
                                mutate_num_insert = (mutate_x_insert - 1) + (mutate_y_insert - 1) * x; % ��դ��λ��ת��Ϊ���
                                mutate_generate_out_pop = [mutate_generate_out_pop(1, 1:mutate_generate_i), mutate_num_insert, mutate_generate_out_pop(1, mutate_generate_i+1:end)]; % ��դ����Ų���single������
                            % ������
                            elseif  G(mutate_y_insert - 1, mutate_x_insert) == 0 && ((mutate_x_insert - 1) + (mutate_y_insert - 2) * x ~= mutate_generate_out_pop(1, mutate_generate_i)) && ((mutate_x_insert - 1) + (mutate_y_insert-2) * x ~= mutate_generate_out_pop(1, mutate_generate_i+1))
                                mutate_y_insert = mutate_y_insert - 1;
                                mutate_num_insert = (mutate_x_insert - 1) + (mutate_y_insert - 1) * x; % ��դ��λ��ת��Ϊ���
                                mutate_generate_out_pop = [mutate_generate_out_pop(1, 1:mutate_generate_i), mutate_num_insert, mutate_generate_out_pop(1, mutate_generate_i+1:end)]; % ��դ����Ų���single������
                            % ���������ȥ��·��
                            else                
                                mutate_generate_out_pop = [];
                                break
                            end    
                        end        
                        mutate_x_next = mutate_x_insert;
                        mutate_y_next = mutate_y_insert;
                        mutate_max_iteration = mutate_max_iteration + 1;
                        if mutate_max_iteration > 20000 % ���ִ��20000�Σ���Ϊƽ��·������������
                            mutate_generate_out_pop = [];
                            break
                        end    
                    end    
                    if isempty(mutate_generate_out_pop)
                        break
                    end
                    [~, mutate_single_path_num] = size(mutate_generate_out_pop);
                    mutate_generate_i = mutate_generate_i + 1;
                end
                mutate_single_new_pop_slice = mutate_generate_out_pop;
                % @@@@@@@@@@@@@@@@@@
                mutate_max = mutate_max + 1;
                if mutate_max >= 100000
                    break
                end
            end
            if mutate_max >= 100000
                mutate_out_pop{mutate_i, 1} = mutate_in_pop{mutate_i, 1};
            else
                mutate_out_pop{mutate_i, 1} = [mutate_single_new_pop(1, 1:mutate_mpoint(1, 1)-1), mutate_single_new_pop_slice(2:end-1), mutate_single_new_pop(1, mutate_mpoint(1, 2)+1:mutate_m)];
            end
            % single_new_pop_slice�ٴγ�ʼ��
            mutate_single_new_pop_slice = [];
        else
            mutate_out_pop{mutate_i, 1} = mutate_in_pop{mutate_i, 1};
        end
    end
    new_pop2 = mutate_out_pop;
    
    % ������Ⱥ
    new_pop1 = new_pop2;
    % ������Ӧ��ֵ
    %@@@@@@@@@@ ����·������
    % path_value = cal_path_value(new_pop1, x);
    path2_value_in_pop = new_pop1;
    [value2_n, ~] = size(path2_value_in_pop);
    path2_value_out_pop = zeros(1,value2_n);
    for value2_i = 1 : value2_n
        path2_value_single_pop = path2_value_in_pop{value2_i, 1};
        [~, value2_m] = size(path2_value_single_pop);
        for value2_j = 1 : value2_m - 1
            % ��i�����У������ұ��1.2.3...��
            value2_x_now = mod(path2_value_single_pop(1, value2_j), x) + 1; 
            % ��i�����У����ϵ��±����1.2.3...��
            value2_y_now = fix(path2_value_single_pop(1, value2_j) / x) + 1;
            % ��i+1�����С���
            value_x_next = mod(path2_value_single_pop(1, value2_j + 1), x) + 1;
            value_y_next = fix(path2_value_single_pop(1, value2_j + 1) / x) + 1;
            if abs(value2_x_now - value_x_next) + abs(value2_y_now - value_y_next) == 1
                path2_value_out_pop(1, value2_i) = path2_value_out_pop(1, value2_i) + 1;
            else
                path2_value_out_pop(1, value2_i) = path2_value_out_pop(1, value2_i) + sqrt(2);
            end
        end
    end
    path_value = path2_value_out_pop;
    
    %@@@@@@@@@@ ����·��ƽ����
    % path_smooth = cal_path_smooth(new_pop1, x);
    path2_smooth_in_pop = new_pop1;
    [smooth2_n, ~] = size(path2_smooth_in_pop);
    path2_smooth_out_pop = zeros(1,smooth2_n);
    for smooth2_i = 1 : smooth2_n
        path2_smooth_single_pop = path2_smooth_in_pop{smooth2_i, 1};
        [~, smooth2_m] = size(path2_smooth_single_pop);
        for smooth2_j = 1 : smooth2_m - 2
            % ��i�����У������ұ��1.2.3...��
            smooth2_x_now = mod(path2_smooth_single_pop(1, smooth2_j), x) + 1; 
            % ��i�����У����ϵ��±����1.2.3...��
            smooth2_y_now = fix(path2_smooth_single_pop(1, smooth2_j) / x) + 1;
            % ��i+1�����С���
            % smooth2_x_next1 = mod(path2_smooth_single_pop(1, smooth2_j + 1), x) + 1;
            % smooth2_y_next1 = fix(path2_smooth_single_pop(1, smooth2_j + 1) / x) + 1;
            % ��i+2�����С���
            smooth2_x_next2 = mod(path2_smooth_single_pop(1, smooth2_j + 2), x) + 1;
            smooth2_y_next2 = fix(path2_smooth_single_pop(1, smooth2_j + 2) / x) + 1;
            %path_smooth(1, i) = path_smooth(1, i) + abs(atan(abs(x_now - x_next1)/abs(y_now - y_next1))-atan(abs(x_next2 - x_next1)/abs(y_next2 - y_next1)));
            %a2 = (x_now - x_next1)^2 + (y_now - y_next1)^2;
            %b2 = (x_next2 - x_next1)^2 + (y_next2 - y_next1)^2;
            c2 = (smooth2_x_now - smooth2_x_next2)^2 + (smooth2_y_now - smooth2_y_next2)^2;
            %angle = (a2 + c2 - b2) / (2 * sqrt(a2) *  sqrt(c2));
            if c2 < 8 && c2 > 4
                path2_smooth_out_pop(1, smooth2_i) = path2_smooth_out_pop(1, smooth2_i) + 5;
            elseif c2 <= 4 && c2 > 1
                path2_smooth_out_pop(1, smooth2_i) = path2_smooth_out_pop(1, smooth2_i) + 30;
            elseif    c2 <= 1
                path2_smooth_out_pop(1, smooth2_i) = path2_smooth_out_pop(1, smooth2_i) + 5000;
            end
        end
    end
    path_smooth = path2_smooth_out_pop;
    
    fit_value = perlenth .* path_value .^ -1 + persmooth .* path_smooth .^ -1;
    mean_path_value(1, global2_i) = mean(path_value);
    [~, global_m] = max(fit_value);
    min_path_value(1, global2_i) = path_value(1, global_m);
end
% ��ÿ�ε���ƽ��·�����Ⱥ�����·������ͼ
figure(1)
plot(1:max_gen,  mean_path_value, 'r')
hold on;
title(['perlenth = ', num2str(perlenth)', '��persmooth = ',num2str(persmooth)','���Ż�����ͼ']); 
xlabel('��������'); 
ylabel('·������');
plot(1:max_gen, min_path_value, 'b')
legend('ƽ��·������', '����·������');
min_path_value(1, end)
% �ڵ�ͼ�ϻ�·��
[~, min_index] = max(fit_value);
min_path = new_pop1{min_index, 1};
figure(2)
hold on;
title(['perlenth = ', num2str(perlenth)', '��persmooth = ',num2str(persmooth)','�Ŵ��㷨�������˶��켣']); 
xlabel('����x'); 
ylabel('����y');
%@@@@@@@@@@ DrawMap(G);
draw_b = G;
draw_b(end+1,end+1) = 0;
colormap([1 1 1;0 0 0]);  % ������ɫ
pcolor(0.5:size(G,2) + 0.5, 0.5:size(G,1) + 0.5, draw_b); % ����դ����ɫ
set(gca, 'XTick', 1:size(G,1), 'YTick', 1:size(G,2));  % ��������
axis image xy;  % ��ÿ��������ʹ����ͬ�����ݵ�λ������һ��

[~, min_path_num] = size(min_path);
for draw_i = 1:min_path_num
    % ·���������У������ұ��1.2.3...��
    x_min_path(1, draw_i) = mod(min_path(1, draw_i), x) + 1; 
    % ·���������У����ϵ��±����1.2.3...��
    y_min_path(1, draw_i) = fix(min_path(1, draw_i) / x) + 1;
end
hold on;
plot(x_min_path, y_min_path, 'r')