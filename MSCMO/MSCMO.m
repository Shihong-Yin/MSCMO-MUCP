classdef MSCMO < ALGORITHM
% <multi> <real/integer/label/binary/permutation> <constrained>
% Multi-stage constrained multi-objective evolutionary algorithm
% cp --- 5 --- Decrease trend of the dynamic constraint boundary
    
    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            cp = Algorithm.ParameterSet(5); % Set the control parameter cp, which determines the rate of reduction of the constraint boundary.
            %% Generate the random population
            Population1 = Problem.Initialization(); % Generate two random initial populations for the optimization process.
            Population2 = Problem.Initialization();            
            %% Calculate the initial dynamic constraint boundary
            [~, nCon]               = size(Population1.cons); % Get the number of constraints for the problem.
            % Calculate the initial dynamic constraint boundary for both populations.
            [initialE1, ~]          = max(max(0,Population1.cons), [], 1); % Boundary for Population1
            [initialE2, ~]          = max(max(0,Population2.cons), [], 1); % Boundary for Population2
            % If no constraints are active, set boundary values to 1 to avoid division by zero.
            initialE1(initialE1==0) = 1;
            initialE2(initialE2==0) = 1;
            % Initialize the dynamic constraint boundaries.
            epsn1                   = initialE1;
            epsn2                   = initialE2;
            % Calculate constraint violations for Population2.
            CV2                     = sum(max(0,Population2.cons),2);
            % Calculate fitness values for Population2 based on objectives and constraint violations.
            Fitness2                = CalFitness([CalSDE(Population2.objs)',CV2]);
            % Initialize an array to track the maximum constraint violation during evolution.
            MaxCV2                  = zeros(ceil(Problem.maxFE/Problem.N),1);
            MaxCV2(1)               = max(CV2);
            % Initialize a counter for adaptive boundary adjustment.
            ASC2                    = 0;
            arch                    = archive([Population1,Population2],Problem.N); % Archive the initial populations.
            %% Optimization
            while Algorithm.NotTerminated(Population1)                
                PopCon1   = max(0,Population1.cons); % Extract constraint values for both populations.
                PopCon2   = max(0,Population2.cons);
                % Check if all solutions in Population1 meet the current constraint boundary.
                if sum(sum(PopCon1<=epsn1,2)==nCon) == length(Population1)
                    % Reduce the constraint boundary dynamically for Population1.
                    epsn1 = ReduceBoundary(initialE1,ceil(Problem.FE/Problem.N),ceil(Problem.maxFE/Problem.N)-1,cp);
                end
                % Update constraint violations for Population2.
                CV2       = sum(PopCon2,2);
                % Track the maximum constraint violation for Population2.
                MaxCV2(ceil(Problem.FE/Problem.N)) = max(CV2);
                % Adjust the dynamic boundary for Population2 based on violation trends.
                if MaxCV2(ceil(Problem.FE/Problem.N))-MaxCV2(ceil((Problem.FE-Problem.N)/Problem.N))>0
                    ASC2  = ASC2 + 1;
                    epsn2 = ReduceBoundary(initialE2,ceil(Problem.FE/Problem.N)-ASC2,ceil(Problem.maxFE/Problem.N)-1,cp);
                elseif sum(sum(PopCon2<=epsn2,2)==nCon) == length(Population2)
                    ASC2  = 0;
                    epsn2 = ReduceBoundary(initialE2,ceil(Problem.FE/Problem.N),ceil(Problem.maxFE/Problem.N)-1,cp);
                end
                % Selection of mating pools using tournament selection.
                MatingPool1 = TournamentSelection(2,Problem.N,sum(max(0,Population1.cons-epsn1),2));
                MatingPool2 = TournamentSelection(2,Problem.N,Fitness2);
                % Generate offspring for both populations using genetic algorithm operators.
                Offspring1  = OperatorGAhalf(Problem,Population1(MatingPool1));
                Offspring2  = OperatorGAhalf(Problem,Population2(MatingPool2));
                % Perform environmental selection for Population1.
                Population1            = EnvironmentalSelection1([Population1,Offspring1,Offspring2],Problem.N,epsn1);
                % Perform environmental selection for Population2 and update fitness.
                [Population2,Fitness2] = EnvironmentalSelection2([Population2,Offspring1,Offspring2],Problem.N,epsn2);                
                % Archive non-dominated and feasible solutions.
                arch = [arch,Population1,Population2];
                 [~, Unduplicated] = unique(arch.objs,'rows');
                arch = arch(Unduplicated);
                arch = archive(arch,Problem.N);
                % If the maximum evaluations are reached, save Population1 as the final archive.
                if Problem.FE >= Problem.maxFE
                    Population1 = arch;
                end                        
            end
        end
    end
end