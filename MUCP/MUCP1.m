classdef MUCP1 < PROBLEM
% <multi> <real/binary> <large/constrained/none>
% The multi-objective cooperative path planning of multiple USVs

    properties(Access = private)
        USV; % USV信息
        ENV; % 环境信息
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            % Parameter setting
            obj.USV.num = 2; % USV的数量
            obj.USV.vel = [10 12;9 13]; % USV的航行速度范围
            obj.USV.startPos = [0 95;5 100]; % 起始点
            obj.USV.goalPos = [100 5;95 0]; % 目标点
            obj.USV.searchRange = [100,100]; % 搜索长、宽范围
            obj.USV.density = 100; % 插值密度
            obj.USV.safeDist = 3; % 最小安全距离

            Threat = [18 51 6.5; % 威胁区域：圆心坐标a b，半径r
                65 49 8;
                50 77 7.5;
                80 35 8;
                79 69 6.5;
                29 35 7;
                50 30 7.5;
                28 75 7];
            obj.ENV.threat = Threat;

            Nodes = 10; % 导航节点数

            obj.M        = 3;
            obj.D        = 2*Nodes*obj.USV.num;
            lower = [];
            upper = [];
            encoding = [];
            for i = 1:obj.USV.num
                d = pdist([obj.USV.startPos(i,:); obj.USV.goalPos(i,:)]);
                r = d / (Nodes+1);
                scale = [1:1:Nodes/2, Nodes/2:-1:1];
                lower = [lower, zeros(1,Nodes), -r*scale];
                upper = [upper, ones(1,Nodes), r*scale];
                encoding = [encoding, 4+zeros(1,Nodes), 1+zeros(1,Nodes)];
            end
            obj.lower    = lower;
            obj.upper    = upper;
            obj.encoding = encoding; % 混合编码
        end
        %% Evaluate objective values
        function Population = Evaluation(obj,varargin)
            PopDec = varargin{1};
            Pop = PopDec;
            USV_num = obj.USV.num; % USV的数量
            Density = obj.USV.density; % 插值密度
            ThreatArea = obj.ENV.threat;

            Fitness = zeros(size(Pop, 1), obj.M);
            Constraint = zeros(size(Pop, 1), 3*USV_num + 2);

            temp_USV = obj.USV; % 为了减少通信开销的变量
            objM = obj.M;
            range = obj.USV.searchRange(1);
            USV_vel = obj.USV.vel;
            safeDist = obj.USV.safeDist;
            parfor m = 1:size(Pop, 1)
                Agent = Pop(m,:);
                multi_Path = decodePath(temp_USV, Agent); % 将搜索代理解码为航行路径
                %%% 单个USV的路径约束
                singleFitness = zeros(USV_num, objM);
                Length = zeros(1,USV_num);
                outCon = zeros(1,USV_num);
                threatCon = zeros(1,USV_num);
                smoothCon = zeros(1,USV_num);

                for k = 1:USV_num
                    Path_USV_k = multi_Path{k};
                    x1 = Path_USV_k(1,:);  y1 = Path_USV_k(2,:); % 原始坐标系中的插值节点

                    % 1. 计算路径长度
                    lengthCost = sum(sqrt(diff(x1).^2 + diff(y1).^2));

                    % 2. 是否出界
                    a1 = x1(2:end-1) < 0 | x1(2:end-1) > range;
                    a2 = y1(2:end-1) < 0 | y1(2:end-1) > range;
                    outCon(1,k) = sum(a1 | a2); % 约束条件

                    % 3. 计算威胁代价
                    ThreatCon = zeros(1,length(x1));
                    ThreatValue = zeros(1,length(x1));
                    for i = 1:length(x1)
                        Threat_A = 0;
                        Threat_B = 0;
                        for j = 1:size(ThreatArea, 1)
                            x0 = ThreatArea(j, 1:2); % 球心
                            r0 = ThreatArea(j, 3); % 半径
                            d = sqrt(sum((Path_USV_k(:,i)-x0').^2));
                            if d < r0
                                Threat_A = Threat_A + (r0-d);
                            else
                                Threat_B = Threat_B + 1/(1+(d-r0).^2); % 超参数beta = 2
                            end
                        end
                        ThreatCon(i) = Threat_A;
                        ThreatValue(i) = Threat_B;
                    end
                    threatCon(1,k) = sum(ThreatCon,2); % 约束条件
                    threatCost = sum(ThreatValue,2);

                    % 4. 路径平滑度
                    theta_max = pi/6;
                    path_direction = diff(Path_USV_k,1,2);
                    turn_angles = zeros(1, size(path_direction, 2) - 1); % 初始化转向角数组
                    % 计算相邻两列向量在 xoy 平面的夹角
                    for i = 1:size(path_direction, 2) - 1
                        vector1 = path_direction(:, i);
                        vector2 = path_direction(:, i + 1);
                        cos_theta = dot(vector1, vector2) / (norm(vector1)*norm(vector2));
                        cos_theta = min(1, max(-1, cos_theta));  % 将 cos_theta 限制在 [-1, 1] 范围内
                        turn_angles(i) = acos(cos_theta); % 计算夹角（弧度）
                    end
                    smoothCon(1,k) = sum(max(turn_angles-theta_max, 0)); % 约束条件
                    smoothCost = sum(turn_angles);

                    % 计算单个USV的路径代价，是否满足约束等
                    singleFitness(k,1) = lengthCost;
                    singleFitness(k,2) = threatCost;
                    singleFitness(k,3) = smoothCost;
                    Length(1,k) = lengthCost; % 记录每条路径的长度，用来评估协同约束
                end

                %%% 多个USV的协同约束
                % 1. 协同代价
                % 可协同时间区间
                Time_interval = zeros(USV_num,2);
                for k = 1:USV_num
                    Time_interval(k,1) = Length(1,k) ./ USV_vel(k,2); % 最短的时间
                    Time_interval(k,2) = Length(1,k) ./ USV_vel(k,1); % 最长的时间
                end
                Ideal_time = max(Time_interval(:,1)); % 理想的航行时间
                cooperativeCon = Ideal_time - min(Time_interval(:,2)); % 约束条件
                cooperativeCost = Ideal_time;

                % 2. 碰撞代价
                if cooperativeCon <= 0 % 满足此条件意味着，多架USV可以同时到达目标位置
                    flight_vel = zeros(1,USV_num);
                    PathU = zeros(Density, 2, USV_num);
                    for k = 1:USV_num
                        Path_USV_k = multi_Path{k};
                        flight_vel(1,k) = Length(1,k) ./ Ideal_time; % 计算出第k架USV的实际航行速度
                        % 计算每个位置点的实际时间
                        time_intervals = [0, cumsum(sqrt(diff(Path_USV_k(1,:)).^2 + diff(Path_USV_k(2,:)).^2))] ./ flight_vel(1,k);
                        % 用线性插值去近似第k架USV的路径
                        uniform_time_interval = linspace(0,Ideal_time,Density); % 近似线性插值节点
                        PathU(:,:,k) = interp1(time_intervals', Path_USV_k', uniform_time_interval', 'linear'); % 使用interp1函数进行插值
                    end
                    % 检查碰撞
                    distances = zeros(Density, USV_num*(USV_num-1)/2); % 多个路径点之间的欧式距离
                    for i = 1:Density
                        vectors = reshape(PathU(i,:,:), 2, USV_num)';
                        distances(i,:) = pdist(vectors); % 计算每对向量之间的欧式距离
                    end
                    collisionCon = sum(max(safeDist - min(distances,[],2), 0)); % 约束条件
                else
                    collisionCon = 1e2; % 不可协同，设置碰撞约束
                end
                multiFitness = cooperativeCost;

                %%% 计算最终的适应度
                Fitness(m, :) = sum(singleFitness,1) + multiFitness;
                Constraint(m, :) = [outCon, threatCon, smoothCon, cooperativeCon, collisionCon];
            end

            % Objective function
            PopObj = Fitness;
            % Constraints
            PopCon = Constraint;
            Population = SOLUTION(PopDec,PopObj,PopCon,varargin{2:end});
            obj.FE     = obj.FE + length(Population);
        end
        %% 定义一个公共方法，用于获取案例信息
        function [USV, ENV] = getCaseInfo(obj)
            USV = obj.USV; % USV信息
            ENV = obj.ENV; % 环境信息
        end

        %% Generate a point for hypervolume calculation
        function R = GetOptimum(obj,~)
            R = [379, 39, 21]; % 参考点先不考虑，这个点可以用来计算HV
        end
    end
end

%% 路径解码函数
function multi_Path = decodePath(USV,Agent)
% Density 插值密度
USV_num = USV.num;
Density = USV.density; % 插值密度
startPos = USV.startPos;
goalPos = USV.goalPos;
% 使用 reshape 函数将行向量重新构造为USV_num * (2*Nodes)的数组
AgentR = reshape(Agent, size(Agent,2) ./ USV_num, USV_num)';
multi_Path = cell(USV_num, 1);
for i = 1:USV_num
    AgentI = AgentR(i,:);
    AgentM = reshape(AgentI, size(AgentI,2) ./ 2, 2)';
    AgentM(1,:) = round(AgentM(1,:)); % 对第一行变量取整
    axisX = sqrt((goalPos(i,1)-startPos(i,1))^2+(goalPos(i,2)-startPos(i,2))^2); % 距离在XOY上的投影
    Xh = linspace(0, axisX, size(AgentM,2)+2); % 在齐次坐标系中的X坐标
    Path = [1,AgentM(1,:),1; Xh; 0,AgentM(2,:),0];
    Valid_path = Path(2:end,Path(1,:)==1); % 组建有效的路径点
    % 路径插值，分段三次 Hermite 插值多项式
    Xh = Valid_path(1,:);  Yh = Valid_path(2,:);
    xh = linspace(Xh(1), Xh(end), Density);
    yh = pchip(Xh,Yh,xh);
%     y1 = interp1(Xo,Yo,x1,'linear');  z1 = interp1(Xo,Zo,x1,'linear'); % 路径插值
    % 对导航节点进行旋转变换
    % HomoTrans 根据起始点和目标点位置返回齐次坐标变换矩阵
    Trans = [1 0 -startPos(i,1); 0 1 -startPos(i,2); 0 0 1];
    Theta = calculate_angle([1 0], goalPos(i,1:2)-startPos(i,1:2));
    Rot_z = [cos(Theta) sin(Theta) 0; -sin(Theta) cos(Theta) 0; 0 0 1];
    Homo = Rot_z * Trans; % 齐次坐标变换矩阵
    PathM = Homo \ [xh; yh; ones(1, Density)];
    multi_Path{i} = PathM(1:2, :);
end
end

%% 使用atan2计算方向角度
function angle = calculate_angle(v1, v2)
    theta1 = atan2(v1(2), v1(1));
    theta2 = atan2(v2(2), v2(1));
    % 计算夹角
    angle = theta2 - theta1;
    % 将夹角限制在[-pi, pi]范围内
    while angle > pi
        angle = angle - 2 * pi;
    end
    while angle < -pi
        angle = angle + 2 * pi;
    end
end