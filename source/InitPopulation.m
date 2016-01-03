%Omid55
function [ population ] = InitPopulation(A, populationSize,N,m ) 

%disp('InitPopulation ... ');
population = zeros(populationSize,N);
for i = 6 : populationSize
    indices = randperm(N);
    population(i,indices(1:m)) = 1;
end
D = sum(A); [D Di] = sort(D,'descend');
[EBC B] = edge_betweenness_bin(A); [B Bi] = sort(B,'descend');
C = Closeness(A); [C Ci] = sort(C,'descend'); 
E = eccentricity(A); [E Ei] = sort(E,'descend'); 
L = clustering_coef_bu(A); [L Li] = sort(L,'descend'); 
population(1,Di(1:m)) = 1;
population(2,Bi(1:m)) = 1;
population(3,Ci(1:m)) = 1;
population(4,Ei(1:m)) = 1;
population(5,Li(1:m)) = 1;

end


