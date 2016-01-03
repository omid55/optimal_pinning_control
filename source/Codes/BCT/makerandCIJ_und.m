function [CIJ] = makerandCIJ_und(N,K)

% inputs:
%           N = number of vertices
%           K = number of edges
% output:
%           CIJ = random connection matrix, nondirected (symmetrical)
%
% This function generates a random binary CIJ matrix, with size (N,K) and
% no connections on the main diagonal
%
% Olaf Sporns, Indiana University, 2007/2008

ind = triu(~eye(N));
i = find(ind);
rp = randperm(length(i));
irp = i(rp);

CIJ = zeros(N);
CIJ(irp(1:K)) = 1;
CIJ = CIJ+CIJ';         % symmetrize
