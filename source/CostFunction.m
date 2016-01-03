%Omid55
function [ R ] = CostFunction( delta,L,weight )

S = L + weight*diag(delta);
lambda = eig(S);
R = lambda(end) / lambda(1);

end

