function [f F]=motif4struct_bin(A)
%[f F]=motif4struct_bin(G); counts structural motif occurence
%
%Input: binary directed graph G
%Output: binary motif count f; binary motif count per vertex F
%
%Reference: Milo et al., 2002, Science
%
%Mika Rubinov, UNSW, 2007 (last modified July 2008)

persistent M4n ID4
if isempty(ID4)
    load motif34lib M4n ID4                     %load motif data
end

n=length(A);                                    %number of vertices in A
F=zeros(199,n);                                 %motif count of each vertex
f=zeros(199,1);                                 %motif count for whole graph
As=A|A.';                                       %symmetric adjacency matrix

for u=1:n-3                                     %loop u 1:n-2
    V1=[false(1,u) As(u,u+1:n)];                %v1: neibs of u (>u)
    for v1=find(V1)
        V2=[false(1,u) As(v1,u+1:n)];           %v2: all neibs of v1 (>u)
        V2(V1)=0;                               %not already in V1
        V2=V2|([false(1,v1) As(u,v1+1:n)]);     %and all neibs of u (>v1)
        for v2=find(V2)
            vz=max(v1,v2);                      %vz: largest rank node
            V3=([false(1,u) As(v2,u+1:n)]);     %v3: all neibs of v2 (>u)
            V3(V2)=0;                           %not already in V1&V2
            V3=V3|([false(1,v2) As(v1,v2+1:n)]);%and all neibs of v1 (>v2)
            V3(V1)=0;                           %not already in V1
            V3=V3|([false(1,vz) As(u,vz+1:n)]); %and all neibs of u (>vz)
            for v3=find(V3)

                s=uint64(sum(10.^(11:-1:0).*[A(v1,u) A(v2,u) A(v3,u)...
                    A(u,v1) A(v2,v1) A(v3,v1) A(u,v2) A(v1,v2)...
                    A(v3,v2) A(u,v3) A(v1,v3) A(v2,v3)]));
                ind=ID4(s==M4n);
                if nargout==2; F(ind,[u v1 v2 v3])=F(ind,[u v1 v2 v3])+1; end
                f(ind)=f(ind)+1;
            end
        end
    end
end