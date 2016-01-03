function  [CIJ,K] = makefractalCIJ(mx_lvl,E,sz_cl)

% inputs:
%           mx_lvl   number of hierarchical levels, N = 2^mx_lvl
%           E        connection density fall-off per level
%           sz_cl    size of clusters (power of 2)
% outputs:
%           CIJ      connection matrix
%           K        number of connections present in the output CIJ
%
% NOTE: 
% Clusters have by default a connection density of 1
% Connection density decays as 1/(E^n), with n = index of hierarchical level
%
% Olaf Sporns, Indiana University, 2005/2007

% make a stupid little template
t = ones(2).*2;

% compute N and cluster size
N = 2^mx_lvl;
sz_cl = sz_cl-1;

n = [0 0 0:mx_lvl-3];

for lvl=1:mx_lvl-1
   CIJ = ones(2^(lvl+1),2^(lvl+1));
   group1 = [1:size(CIJ,1)/2];
   group2 = [size(CIJ,1)/2+1:size(CIJ,1)];
   CIJ(group1,group1) = t;
   CIJ(group2,group2) = t;
   CIJ = CIJ+ones(size(CIJ,1),size(CIJ,1));
   t = CIJ;
end;
s = size(CIJ,1);
CIJ = CIJ-ones(s,s)-mx_lvl.*eye(s);

% assign connection probablities
ee = mx_lvl-CIJ-sz_cl;
ee = (ee>0).*ee;
prob = (1./(E.^ee)).*(ones(s,s)-eye(s));
CIJ = (prob>rand(N));

% count connections
K = sum(sum(CIJ));

