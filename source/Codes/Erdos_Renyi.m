function A = Erdos_Renyi(N, P)

A = zeros(N,N); PP = rand(N,N);
A(PP < P) = 1;
A = triu(A); A = A+A';
for i = 1:N; A(i,i) = 0; end;