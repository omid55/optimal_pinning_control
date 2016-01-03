function R=randmio_dir(R, ITER)
%R=randmio_dir(G,ITER); randomized graph R, with equivalent degree
%sequence to the original weighted directed graph G.
%
%Each edge is rewired (on average) ITER times. The out-strength (but not
%the in-strength) distribution is preserved for weighted graphs.
%
%Rewiring algorithm: Maslov and Sneppen (2002) Science 296:910
%
%Mika Rubinov, UNSW, 2007 (last modified July 2008).

[i j]=find(R);
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

        %rewiring condition
        if ~(R(a,d) || R(c,b))
            R(a,d)=R(a,b); R(a,b)=0;
            R(c,b)=R(c,d); R(c,d)=0;

            j(e1) = d;          %reassign edge indices
            j(e2) = b;
            break;
        end %rewiring condition
    end %while not rewired
end %iterations