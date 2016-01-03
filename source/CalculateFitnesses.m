%Omid55
function [ fitnesses ] = CalculateFitnesses( population,L,fitnesses,weight )

%% Fitness Caculation
parfor i=1:size(population,1)
    fitnesses(i) = CostFunction(population(i,:),L,weight);
end

end

