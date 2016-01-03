function G=SW_SmallWorld_fixed(n, k, p)

% SmallWorld generates a Strogats-Watts version of small-world random graphs.
% the disconnected nodes are randomely connected to other nodes. 
% Inputs:
% n: number of nodes
% k: k-nearest neighbors' connection (regular graph)
% p: probability of rewiring
%
% Output:
% G: A structure implemented as data structure in this as well as other
%    graph theory algorithms.
%    G.Adj   - is the adjacency matrix (1 for connected nodes, 0 otherwise)
%    G.x and G.y -   are row vectors of size nv wiht the (x,y) coordinates
%                    of each node of G
%    G.n    - number of vertices in G
%    G.e    - number of edges in G
% Created by Mahdi Jalili on 06/12/2006
% mahdi.jalili@epfl.ch
A = zeros(n,n);
for i = 1:k
    A = sparse(A+diag(ones(1,length(diag(A,i))),i)+diag(ones(1,length(diag(A,n-i))),n-i));
end
A=A+A';
B = -diag(ones(1, n), 0);
A = A+B;

for q = 1:n-k
    [i1,j1] = find( A(q, :) == 1 ); [i2,j2] = find( A(q, :) == 0 );
    for i = 1:length(i1)
        jj = ceil(length(j2)*rand); 
        % making the random distant connections
        prob = rand;
        if (prob < p) && (A(q, j2(jj)) == 0) 
            A(q, j1(i)) = 0; A(j1(i), q) = 0;
            A(q, j2(jj)) = 1; A(j2(jj), q) = 1; 
        end
    end
end

A=A-B;
% randomely conencted the disconencted nodes to other nodes
D=sum(A);
Num=randint(1, n, n-1); 
for i=1:n-1; if Num(i)==i; Num(i)=Num(i+1); end; end;
if Num(n)==n; Num(n)=Num(n-1); end;
for i=1:n
    if D(i)==0
        A(i,Num(i))=1; A(Num(i),i)=1;
    end
end

center=[0,0];
theta=linspace(0,2*pi,n+1);
rho=ones(1,n+1);%fit radius and nv
[X,Y] = pol2cart(theta',rho');
X=X+center(1);
Y=Y+center(2);
x=X(1:end-1)*10;
y=Y(1:end-1)*10;
G=A;
G=struct('Adj',G,'x',x','y',y','n',n,'e',nnz(G));