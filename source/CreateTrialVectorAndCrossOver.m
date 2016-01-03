%Omid55
function [ nextGeneration,nextFitnesses ] = CreateTrialVectorAndCrossOver( population,fitnesses,Beta,Pr,Nv,L,m,weight )

%disp('Create Trial Vector And Doing CrossOver ... ');
populationCount = size(population,1);
genesNum = size(population,2);
nextGeneration = zeros(size(population));
nextFitnesses = zeros(size(fitnesses));
for i = 1 : populationCount
    
    Xi = population(i,:);
    fXi = fitnesses(i,:);
    
%     diffs = zeros(1,size(population,2));
%     list = 1:populationCount;
%     list(i) = [];
%     for j=1:Nv
%         idx = randperm(length(list));
%         selected = list(idx(1:2));
%         xs = population(selected,:);
%         diffs = diffs + xs(1,:) - xs(2,:);
%         list(list == selected(1)) = [];
%         list(list == selected(2)) = [];
%     end
% 
%     index = list(randi(length(list)));
%     x_selected = population(index,:);
%     Ui = x_selected + Beta * diffs;

    % ----------- for Nv = 1 then -----------------
    idx = randperm(populationCount);
    xs = population(idx(1:3),:);
    Ui = xs(1,:) + Beta * (xs(2,:) - xs(3,:));
    % ------------------------------------------------
    
    % Binomial Crossover
    jStar = randi(genesNum,1);
    selections = find(rand(genesNum,1) < Pr);
    J = union(selections,jStar);
        
    child = Xi;
    child(J) = Ui(J);    
    
    % make child valid
    child = round(child);
    child(child<0) = 0;
    child(child>1) = 1;
    list = find(child == 1);
    if length(list) > m
        ind = randperm(length(list));
        child(list(ind(1:length(list)-m))) = 0;
    else
        if length(list) < m
            zl = find(child == 0);
            ind = randperm(length(zl));
            child(zl(ind(1:m-length(list)))) = 1;
        end
    end
    
    childFitness = CostFunction(child,L,weight);

    if childFitness <= fXi
        nextGeneration(i,:) = child;
        nextFitnesses(i,:) = childFitness;
    else
        nextGeneration(i,:) = Xi;
        nextFitnesses(i,:) = fXi;
    end

end


end
