function A = Erdos_Renyi_directed(N, P)

A = zeros(N,N); PP = rand(N,N);
A(PP < P) = 1;
for i = 1:N; A(i,i) = 0; end;