function BC=betweenness_bin(G)
%BC=betweenness_bin(G); betweenness centrality BC, for a binary directed graph G
%
%Betweenness may be normalised to [0,1] via BC/[(N-1)(N-2)]
%
%Algorithm of Kintali, generalised to directed and disconnected graphs
%http://www.cc.gatech.edu/~kintali/papers/bc.pdf
%
%Mika Rubinov, UNSW, 2007 (last modified July 2008)

n=length(G);                %number of nodes
I=eye(n)~=0;                %logical identity matrix
d=1;                     	%path length
NPd=G;                      %number of paths of length |d|
NSPd=NPd;                  	%number of shortest paths of length |d|
NSP=NSPd; NSP(I)=1;        	%number of shortest paths of any length
L=NSPd; L(I)=1;           	%length of shortest paths

%calculate NSP and L
while find(NSPd,1);
    d=d+1;
    NPd=NPd*G;
    NSPd=NPd.*(L==0);
    NSP=NSP+NSPd;
    L=L+d.*(NSPd~=0);
end
L(~L)=inf; L(I)=0;          %L for disconnected vertices is inf
NSP(~NSP)=1;                %NSP for disconnected vertices is 1

Gt=G.';
DP=zeros(n);            	%vertex on vertex dependency
diam=d-1;                  	%graph diameter

%calculate DP
for d=diam:-1:2
    DPd1=(((L==d).*(1+DP)./NSP)*Gt).*((L==(d-1)).*NSP);
    DP=DP + DPd1;       %DPd1: dependencies on vertices |d-1| from source
end

BC=sum(DP,1);               %compute betweenness