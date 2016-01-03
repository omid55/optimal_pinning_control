%Omid55
% A = adjacency matrix
% m = the number of control nodes
% withFigure = 1 or 0 for showing optimization trend
function [BestNodes , BestCost] = Optimize( A , m , weight, withFigure )

% Optimization parameters
% ------------------------------------------------------------------------------------------------
%% Nonconvex Optimization with DE
%FOR A BETTER OPTIMIZATION YOU CAN ::::::::::::
DE_config.Number_Of_Population = size(A,1);   %  INCREASE THIS SIZE
DE_config.Coef = 0.7;                                        %  INCREASE THIS PARAMETER NEAR TO 1


%Other Params
DE_config.MaxIteration = 1000;
DE_config.Beta = 0.8;
DE_config.Pr = 0.6;
DE_config.Nv = 1;
DE_config.withFigure = withFigure;
% ------------------------------------------------------------------------------------------------


%% Laplacian calculation
D = diag(sum(A));
L  = D - A;
%% DE/rand/1/binomial
[BestChromosome,BestCost] = DE_Optimizer(A, DE_config,L,m,weight);
if length(find(BestChromosome == 1)) ~= m
    disp('NOT VALID');
end
BestNodes = find(BestChromosome == 1);

end
