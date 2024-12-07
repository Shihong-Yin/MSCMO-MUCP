function Population = archive(Population,N)
% Select feasible and non-dominated solutions by SPEA2

    %% Select feasible solutions
    fIndex           = all(Population.cons <= 0,2);
    Population       = Population(fIndex);
    if isempty(Population)
        return
    else
        Fitness = CalFitness(Population.objs);
        Next    = Fitness < 1;
        if sum(Next) > N
            Del  = Truncation(Population(Next).objs,sum(Next)-N);
            Temp = find(Next);
            Next(Temp(Del)) = false;
        end
        Population = Population(Next);
    end
end

function Del = Truncation(PopObj,K)
% Select part of the solutions by truncation

    %% Truncation
    Distance = pdist2(PopObj,PopObj);
    Distance(logical(eye(length(Distance)))) = inf;
    Del = false(1,size(PopObj,1));
    while sum(Del) < K
        Remain   = find(~Del);
        Temp     = sort(Distance(Remain,Remain),2);
        [~,Rank] = sortrows(Temp);
        Del(Remain(Rank(1))) = true;
    end
end