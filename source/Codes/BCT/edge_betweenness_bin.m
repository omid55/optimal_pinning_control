function [EBC BC]=edge_betweenness_bin(G)
%EBC=edge_betweenness_bin(G); edge betweenness centrality EBC, for a binary graph G
%[EBC BC]=edge_betweenness_bin(G), also outputs vertex betweenness centrality BC
%
%Betweenness may be normalised to [0,1] via EBC/[(N-1)(N-2)]
%
%Brandes's modified breadth-first search; J Math Sociol (2001) 25:163-177.
%
%Mika Rubinov, UNSW, 2007 (last modified July 2008).

n=length(G);
BC=zeros(n,1);                  %vertex betweenness
EBC=zeros(n);                   %edge betweenness

for u=1:n
    D=false(1,n); D(u)=1;      	%distance from u
    NP=zeros(1,n); NP(u)=1;     %number of paths from u
    P=false(n);                 %predecessors
    Q=zeros(1,n); q=n;          %order of non-increasing distance

    Gu=G;
    V=u;
    while V
        Gu(:,V)=0;              %remove remaining in-edges
        for v=V
            Q(q)=v; q=q-1;
            W=find(Gu(v,:));                %neighbours of v
            for w=W
                if D(w)
                    NP(w)=NP(w)+NP(v);      %NP(u->w) sum of old and new
                    P(w,v)=1;               %v is a predecessor
                else
                    D(w)=1;
                    NP(w)=NP(v);            %NP(u->w) = NP of new path
                    P(w,v)=1;               %v is a predecessor
                end
            end
        end
        V=find(any(Gu(V,:),1));
    end
    if ~all(D)                              %if some vertices unreachable,
        Q(1:q)=find(~D);                    %...these are first-in-line
    end

    DP=zeros(n,1);                          %dependency
    for w=Q(1:n-1)
        BC(w)=BC(w)+DP(w);
        for v=find(P(w,:))
            DPvw=(1+DP(w)).*NP(v)./NP(w);
            DP(v)=DP(v)+DPvw;
            EBC(v,w)=EBC(v,w)+DPvw;
        end
    end
end