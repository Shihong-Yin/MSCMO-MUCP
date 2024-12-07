function Population = EnvironmentalSelection1(Population,N,epsn)
% The environmental selection of MSCMO
   
    [~, nCon] = size(Population.cons);
    PopCon    = max(0,Population.cons);
    if sum(sum(PopCon<=epsn, 2)==nCon) > N       
        tmp        = sum(PopCon<=epsn, 2)==nCon;
        Population = Population(1:end, tmp);
        CV         = sum(max(0,Population.cons),2);
        Fitness    = CalFitness(Population.objs,CV);
        Next = Fitness < 1;
        if sum(Next) < N
            [~,Rank] = sort(Fitness);
            Next(Rank(1:N)) = true;
        elseif sum(Next) > N
            Del  = Truncation(Population(Next).objs,sum(Next)-N);
            Temp = find(Next);
            Next(Temp(Del)) = false;
        end
    else       
        CV         = sum(max(0,Population.cons),2);
        [~, rank]  = sort(CV);
        Next       = rank(1:N);       
    end   
    Population = Population(Next);
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