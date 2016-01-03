function C=clustering_coef_bu(G)
%C=clustering_coef_bu(G); clustering coefficient C, for binary undirected graph G
%
%Reference: Watts and Strogatz, 1998, Nature 393:440-442
%
%Mika Rubinov, UNSW, 2007 (last modified September 2008)

n=length(G);
C=zeros(n,1);

for u=1:n
    V=find(G(u,:));
    k=length(V);
    if k>=2;                %degree must be at least 2
        S=G(V,V);
        C(u)=sum(S(:))/(k^2-k);
    end
end