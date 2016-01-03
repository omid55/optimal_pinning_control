function [I Q F]=motif3funct_wei(W)
%[I Q F]=motif3funct_wei(W); weighted functional motif metrics.
%
%Input: weighted graph W (all weights [0,1]).
%Output by node: total intensity I, total coherence Q, motif count F.
%Average intensity and coherence may be obtained as I./F and Q./F.
%
%References: Onnela et al. 2005, Phys Rev E 71:065103;
%Sporns and Kotter, PLoS Biology 2004, 2:e369
%
%Mika Rubinov, UNSW, 2007 (last modified July 2008)

persistent M3 ID3 N3
if isempty(N3)
    load motif34lib M3 ID3 N3             	%load motif data
end

n=length(W);                                %number of vertices in W
I=zeros(13,n);                              %intensity
Q=zeros(13,n);                              %coherence
F=zeros(13,n);                          	%frequency

A=1*(W~=0);                                 %adjacency matrix
As=A|A.';                                   %symmetrized adjacency

for u=1:n-2                               	%loop u 1:n-2
    V1=[false(1,u) As(u,u+1:n)];         	%v1: neibs of u (>u)
    for v1=find(V1)
        V2=[false(1,u) As(v1,u+1:n)];       %v2: all neibs of v1 (>u)
        V2(V1)=0;                           %not already in V1
        V2=([false(1,v1) As(u,v1+1:n)])|V2; %and all neibs of u (>v1)
        for v2=find(V2)
            w=[W(v1,u) W(v2,u) W(u,v1) W(v2,v1) W(u,v2) W(v1,v2)];
            a=[A(v1,u);A(v2,u);A(u,v1);A(v2,v1);A(u,v2);A(v1,v2)];
            ind=(M3*a)==N3;                 %find all contained isomorphs
            m=sum(ind);                     %number of isomorphs

            M=M3(ind,:).*repmat(w,m,1);
            id=ID3(ind);
            l=N3(ind);

            x=sum(M,2)./l;                  %arithmetic mean
            M(M==0)=1;                      %enable geometric mean
            i=prod(M,2).^(1./l);            %intensity
            q=i./x;                         %coherence

            [idu j]=unique(id);             %unique motif occurences
            j=[0;j];
            mu=length(idu);                 %number of unique motifs
            i2=zeros(mu,1);
            q2=i2; f2=i2;

            for h=1:mu                      %for each unique motif
                i2(h)=sum(i(j(h)+1:j(h+1)));    %sum all intensities,
                q2(h)=sum(q(j(h)+1:j(h+1)));    %coherences
                f2(h)=j(h+1)-j(h);              %and frequencies
            end

            %then add to cumulative count
            I(idu,[u v1 v2])=I(idu,[u v1 v2])+[i2 i2 i2];
            Q(idu,[u v1 v2])=Q(idu,[u v1 v2])+[q2 q2 q2];
            F(idu,[u v1 v2])=F(idu,[u v1 v2])+[f2 f2 f2];
        end
    end
end