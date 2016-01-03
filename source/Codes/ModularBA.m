function A = ModularBA(m,n,k,B,P_inter)

% ModularSW generates the adjacenecy matrix of modular scale-free networks 
% P_inter being the ptobability of inter-modular connections. 

A = zeros(m*n,m*n);
for i = 1:m;
    AA = BAnet(n, k, B); 
    A((i-1)*n+1:i*n,(i-1)*n+1:i*n) = AA;
end;
for i = 1:m;
    for p = (i-1)*n+1:i*n-1;
        for q = p+1:i*n;
            R = rand; if R < P_inter & A(p,q) == 1; A(p,q) = 0; A(q,p) = 0;
                B = [p,q]; B = B(randint(1,1,[1,length(B)]));
                BB = randint(1,1,[1,m*n]); 
                if B ~= BB; A(B,BB) = 1; A(B,BB) = 1; end;
            end;
        end;
    end;
end;