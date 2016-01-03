function D=distance_wei(G)
%D=distance_wei(G); distance matrix for a weighted directed graph -
%the mean distance is the characteristic path length.
%
%The input matrix must be a mapping from weight to distance (eg. higher
%correlations may be interpreted as short distances via an inverse mapping).
%
%Dijkstra's Algorithm.
%
%Mika Rubinov, UNSW

%Modification history
%2007: original
%2009-08-04: min() function vectorized

n=length(G);
D=zeros(n); D(~eye(n))=inf;                 %distance matrix

for u=1:n
    S=true(1,n);                            %distance permanence (true is temporary)
    G1=G;
    V=u;
    while 1
        S(V)=0;                             %distance u->V is now permanent
        G1(:,V)=0;                          %no in-edges as already shortest
        for v=V
            W=find(G1(v,:));                %neighbours of shortest nodes
            D(u,W)=min([D(u,W);D(u,v)+G1(v,W)]); %smallest of old/new path lengths
        end

        minD=min(D(u,S));
        if isempty(minD)||isinf(minD),      %isempty: all nodes reached;
            break,                          %isinf: some nodes cannot be reached
        end;

        V=find(D(u,:)==minD);
    end
end