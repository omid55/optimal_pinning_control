function BC=betweenness_wei(G)
%BC=betweenness_wei(G); betweenness centrality BC for weighted directed graph
%
%The input matrix must be a mapping from weight to distance (eg. higher
%correlations may be interpreted as short distances - hence an inverse
%mapping is appropriate in that case).
%
%Betweenness may be normalised to [0,1] via BC/[(N-1)(N-2)]
%
%Brandes's modified Dijkstra's algorithm; J Math Sociol (2001) 25:163-177.
%
%Mika Rubinov, UNSW, 2007 (last modified July 2008)

n=length(G);
% E=find(G); G(E)=1./G(E);        %invert weights
BC=zeros(n,1);                  %vertex betweenness

for u=1:n
    D=inf(1,n); D(u)=0;         %distance from u
    NP=zeros(1,n); NP(u)=1;     %number of paths from u
    S=true(1,n);                %distance permanence (true is temporary)
    P=false(n);                 %predecessors
    Q=zeros(1,n); q=n;          %order of non-increasing distance

    G1=G;
    V=u;
    while 1
        S(V)=0;                 %distance u->V is now permanent
        G1(:,V)=0;              %no in-edges as already shortest
        for v=V
            Q(q)=v; q=q-1;
            W=find(G1(v,:));                %neighbours of v
            for w=W
                Duw=D(v)+G1(v,w);           %path length to be tested
                if Duw<D(w)                 %if new u->w shorter than old
                    D(w)=Duw;
                    NP(w)=NP(v);            %NP(u->w) = NP of new path
                    P(w,:)=0;
                    P(w,v)=1;               %v is the only predecessor
                elseif Duw==D(w)            %if new u->w equal to old
                    NP(w)=NP(w)+NP(v);      %NP(u->w) sum of old and new
                    P(w,v)=1;               %v is also a predecessor
                end
            end
        end

        minD=min(D(S));
        if isempty(minD), break             %all nodes reached, or
        elseif isinf(minD),                 %...some cannot be reached:
            Q(1:q)=find(isinf(D)); break	%...these are first-in-line
        end
        V=find(D==minD);
    end

    DP=zeros(n,1);                          %dependency
    for w=Q(1:n-1)
        BC(w)=BC(w)+DP(w);
        for v=find(P(w,:))
            DP(v)=DP(v)+(1+DP(w)).*NP(v)./NP(w);
        end
    end
end