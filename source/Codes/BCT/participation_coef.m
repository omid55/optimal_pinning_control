function P=participation_coef(A,Ci)
%P=participation_coef(A,Ci); computes nodal "participation coefficient".
%
%Input: (binary) adjacency matrix A, community structure vector Ci.
%Output: participation coef P.
%Output for directed graphs: "out-neighbor" participation coef.
%
%Reference: Guimera R, Amaral L. Nature (2005) 433:895-900.
%
%Mika Rubinov, UNSW, 2008

n=length(A);                        %number of vertices
Ko=sum(A,2);                        %(out)degree
Ko(~Ko)=inf;                        %p_ind=0 if no (out)neighbors
Gc=A*diag(Ci);                      %neighbor community affiliation
Kc2=zeros(n,1);                     %community-specific neighbors

for i=1:max(Ci);
   Kc2=Kc2+(sum(Gc==i,2).^2);
end

P=ones(n,1)-Kc2./(Ko.^2);