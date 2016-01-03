function R=randmio_und(R, ITER)
%R=randmio_und(G,ITER); randomized graph R, with equivalent degree
%sequence to the original weighted undirected graph G.
%
%Each edge is rewired (on average) ITER times. The strength distributions 
%are not preserved for weighted graphs.
%
%Rewiring algorithm: Maslov and Sneppen (2002) Science 296:910
%
%Mika Rubinov, UNSW
%
%Modification History:
%Jun 2007: Original
%Apr 2008: Edge c-d is flipped with 50% probability, allowing to explore
%          all potential rewirings (Jonathan Power, WUSTL)


[i j]=find(tril(R));
K=length(i);
ITER=K*ITER;

for iter=1:ITER
    while 1                                     %while not rewired
        while 1
            e1=ceil(K*rand);
            e2=ceil(K*rand);
            while (e2==e1),
                e2=ceil(K*rand);
            end
            a=i(e1); b=j(e1);
            c=i(e2); d=j(e2);

            if all(a~=[c d]) && all(b~=[c d]);
                break           %all four vertices must be different
            end
        end

        if rand>0.5
            i(e2)=d; j(e2)=c; 	%flip edge c-d with 50% probability
            c=i(e2); d=j(e2); 	%to explore all potential rewirings
        end
        
        %rewiring condition
        if ~(R(a,d) || R(c,b))
            R(a,d)=R(a,b); R(a,b)=0;
            R(d,a)=R(b,a); R(b,a)=0;
            R(c,b)=R(c,d); R(c,d)=0;
            R(b,c)=R(d,c); R(d,c)=0;

            j(e1) = d;          %reassign edge indices
            j(e2) = b;
            break;
        end %rewiring condition
    end %while not rewired
end %iterations