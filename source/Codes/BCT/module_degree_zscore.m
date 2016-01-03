function Z=module_degree_zscore(A,Ci)
%Z=module_degree_zscore(A,Ci); computes "within module degree z-score"
%
%Input: binary adjacency matrix A, community structure vector Ci.
%Output: z-score, Z.
%Output for directed graphs: "out-neighbor" z-score.
%
%Reference: Guimera R, Amaral L. Nature (2005) 433:895-900.
%
%Mika Rubinov, UNSW, 2008

n=length(A);                        %number of vertices
Z=zeros(n,1);
for i=1:max(Ci)
    Koi=sum(A(Ci==i,Ci==i),2);
    Z(Ci==i)=(Koi-mean(Koi))./std(Koi);
end

Z(isnan(Z))=0;