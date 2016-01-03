function [Min,Mout,Mall] = matching_ind(CIJ)

% input:
%           CIJ  = connection/adjacency matrix
% output:
%           Min  = matching index for incoming connections
%           Mout = matching index for outgoing connections
%           Mall = matching index for all connections

% Does not use self- or cross connections for comparison.
% Does not use connections that are not present in BOTH i and j.
% All output matrices are calculated for upper triangular only (symmetrical).
%
% Olaf Sporns, Indiana University, 2002/2007/2008

N = size(CIJ,1);

% compare incoming connections only
Min = zeros(N,N);
for i=1:N-1
    for j=i+1:N
        c1 = CIJ(:,i);
        c2 = CIJ(:,j);
        use = ~(~c1&~c2);
        use(i) = 0;
        use(j) = 0;
        ncon = sum(c1(use))+sum(c2(use));
        if (ncon==0)
            Min(i,j) = 0;
        else
            Min(i,j) = 2*(sum(c1(use)&c2(use))/ncon);
        end;
    end;
end;

% compare outgoing connections only
Mout = zeros(N,N);
for i=1:N-1
    for j=i+1:N
        c1 = CIJ(i,:);
        c2 = CIJ(j,:);
        use = ~(~c1&~c2);
        use(i) = 0;
        use(j) = 0;
        ncon = sum(c1(use))+sum(c2(use));
        if (ncon==0)
            Mout(i,j) = 0;
        else
            Mout(i,j) = 2*(sum(c1(use)&c2(use))/ncon);
        end;
    end;
end;

% compare all (incoming+outgoing) connections
Mall = zeros(N,N);
for i=1:N-1
    for j=i+1:N
        c1 = [CIJ(:,i)' CIJ(i,:)];
        c2 = [CIJ(:,j)' CIJ(j,:)];
        use = ~(~c1&~c2);
        use(i) = 0;  use(i+N) = 0;
        use(j) = 0;  use(j+N) = 0;
        ncon = sum(c1(use))+sum(c2(use));
        if (ncon==0)
            Mall(i,j) = 0;
        else
            Mall(i,j) = 2*(sum(c1(use)&c2(use))/ncon);
        end;
    end;
end;

