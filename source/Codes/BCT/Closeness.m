function C = Closeness(G)
% Vulnerability vector for an undirected and unweighted graph G.
% The dimension of V is the number of edges in the network G.
%
% Mahdi Jalili, STU, 2010.
D = distance_bin(G); D(isfinite(D) == 0) = 0;
C = (size(G,1)-1)./sum(D,2);