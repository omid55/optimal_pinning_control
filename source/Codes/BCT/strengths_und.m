function [str] = strengths_und(CIJ)

% input:  
%         CIJ  = connection/adjacency matrix
% output: 
%         str  = strength for all vertices
%
% Computes the strength for a nondirected weighted matrix.
%
% Olaf Sporns, Indiana University, 2007/2008
% =========================================================================

% compute strengths
str = sum(CIJ);        % strength


