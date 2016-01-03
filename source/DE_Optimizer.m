%Omid55
% DE Optimizer
function [ BestChromosome,BestFitness ] = DE_Optimizer(A, DE_config,L,m,weight )

%% Parameter loading
Number_Of_Population = DE_config.Number_Of_Population;
MaxIteration = DE_config.MaxIteration;
Coef = DE_config.Coef;
Beta = DE_config.Beta;
Pr = DE_config.Pr;
Nv = DE_config.Nv;   % number of difference vectors
withFigure = DE_config.withFigure;
N = size(L,1);


%% Differential Evolution (DE)
disp('--- DE Algorithm --- ');
population = InitPopulation(A, Number_Of_Population,N,m);
fitnesses = zeros(size(population,1),1);
fitnesses = CalculateFitnesses(population,L,fitnesses,weight);
if withFigure == 1
    bests = min(fitnesses);
    means = mean(fitnesses);
end

disp('Started...');
it = 0;
while it < MaxIteration
    %tic;
    
    %lastDe = population;
    
    it = it + 1;
    if withFigure == 1
        bests = [bests; min(fitnesses)];
        means = [means; mean(fitnesses)];
        figure(1);
        plot(0:it,bests,0:it,means);
        legend('Best Cost','Mean Cost');
        xlabel('Generations');
        ylabel('Cost');
        title('Differential Evolution Algorithm');
    end
    
    [population,fitnesses] = CreateTrialVectorAndCrossOver(population,fitnesses,Beta,Pr,Nv,L,m,weight);
    
    %bb = find(fitnesses == min(fitnesses));
    %res = find(population(bb(1),:) == 1);
    %time = toc;
    % disp(['DE: iteration ' num2str(it) ' done in ' num2str(time) 's :  Best = [' num2str(res) ']']);
    
    diff = length(find(fitnesses <= mean(fitnesses)));
    if diff >= Coef * Number_Of_Population                       %if norm(population - lastDe) < Eps
        break;
    end
end

best = find(fitnesses == min(fitnesses));
bestIndex = best(1);
BestChromosome = population(bestIndex,:);
BestFitness = fitnesses(bestIndex);

end

