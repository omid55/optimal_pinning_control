function [kden,N,K] = density_und(CIJ)

% input:  
%           CIJ  = connection/adjacency matrix
% output: 
%           kden = connection density, number of connections present out of all possible (N^2-N)
%           N    = number of vertices
%           K    = number of edges for the entire graph

% Note: Assumes that CIJ is undirected and that there are no self-connections.
% Note: Function always returns average binary density, regardless of
% weights.
%
% Olaf Sporns, Indiana University, 2002/2007/2008

% Modification history:
% 2009-10: K fixed to sum over one half of CIJ [Tony Herdman, SFU]

N = size(CIJ,1);
K = nnz(triu(CIJ));
kden = K/((N^2-N)/2);

