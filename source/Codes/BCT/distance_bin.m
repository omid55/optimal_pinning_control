function D=distance_bin(G)
%D=distance_bin(G); distance matrix for binary undirected graph G
%Mean distance (excluding the main diagonal) equals the characteristic path length
%
%Algebraic shortest path algorithm.
%
%Mika Rubinov, UNSW, 2007 (last modified September 2008).

D=eye(length(G));
n=1;
nPATH=G;                        %n-path matrix
L=(nPATH~=0);                   %shortest n-path matrix

while find(L,1);
    D=D+n.*L;
    n=n+1;
    nPATH=nPATH*G;
    L=(nPATH~=0).*(D==0);
end

D(~D)=inf;                      %disconnected nodes are assigned d=inf;
D=D-eye(length(G));