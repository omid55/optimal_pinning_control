function [id,od,deg] = degrees_dir(CIJ)

% input:  
%         CIJ  = connection/adjacency matrix
% output: 
%         id   = indegree for all vertices
%         od   = outdegree for all vertices
%         deg  = degree for all vertices
%
% Computes the indegree, outdegree, and degree (indegree + outdegree) for a
% directed binary matrix.  Weights are discarded.
%
% Note: Inputs of CIJ are assumed to be on the columns of the matrix.
%
% Olaf Sporns, Indiana University, 2002/2006/2008
% =========================================================================

% ensure CIJ is binary...
CIJ = double(CIJ~=0);

% compute degrees
id = sum(CIJ,1);    % indegree = column sum of CIJ
od = sum(CIJ,2)';   % outdegree = row sum of CIJ
deg = id+od;        % degree = indegree+outdegree


