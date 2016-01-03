function  [CIJ] = maketoeprandCIJ(N,K,s)

% inputs:
%           N        number of vertices
%           K        number of edges
%           s        standard deviation of toeplitz
% outputs:
%           CIJ      connection matrix
%
% makes one CIJ matrix, with size = N,K, that has K connections arranged in
% a toeplitz form.
% NO RING
% no connections on main diagonal
%
% Olaf Sporns, Indiana University, 2005/2007

profile = normpdf([1:N-1],0.5,s);
template = toeplitz([0 profile],[0 profile]);
template = template.*(K./sum(sum(template)));
CIJ = zeros(N);

while ((sum(sum(CIJ)) ~= K))
   CIJ = (rand(N)<template);
end;
