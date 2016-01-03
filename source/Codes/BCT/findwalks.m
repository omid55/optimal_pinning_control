function [Wq,twalk,wlq] = findwalks(CIJ)

% input:  
%           CIJ       connection/adjacency matrix
% outputs: 
%           Wq        3D matrix, with Wq(i,j,q) = number of walks from 
%                     'i' to 'j' of length 'q'.
%           twalk     total number of walks found
%           wlq       walk length distribution as function of 'q'
%
% Uses the powers of the adjacency matrix to produce numbers of walks
% Note that Wq grows very quickly for larger N,K,q.
% Note: Weights are discarded.
%
% Written by Olaf Sporns, Indiana University, 2002/2007/2008

% ensure CIJ is binary...
CIJ = double(CIJ~=0);

N = size(CIJ,1);
Wq = zeros(N,N,N);
CIJpwr = CIJ;
Wq(:,:,1) = CIJ;
for q=2:N
   CIJpwr = CIJpwr*CIJ;
   Wq(:,:,q) = CIJpwr;
end;

% total number of walks
twalk = sum(sum(sum(Wq)));

% walk length distribution
wlq = reshape(sum(sum(Wq)),1,N);

