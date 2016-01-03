function [MATreordered,MATindices,MATcost] = reorderMAT(MAT,H,cost);

% MAT       input matrix (can be anything...)
% H         number of reordering attempts
% cost      either 'line' or 'circ', for shape of lattice (linear or ring)

N = length(MAT);
diagMAT = diag(diag(MAT));
MAT = MAT-diagMAT;

% generate cost function
if (cost=='line');
    profil = fliplr(normpdf([1:N],0,N/2));
end;
if (cost=='circ');
    profil = fliplr(normpdf([1:N],N/2,N/4));
end;
COST = toeplitz(profil,profil);

% initialize lowCOST
lowMATcost = sum(sum(COST.*MAT));
%[x,y] = find(MAT);
%lowMATcost = sum(profil(abs(x-y)));

% keep track of starting configuration
MATstart = MAT;
starta = 1:N;
   
% reorder
for h=1:H
    a = 1:N;
    % choose two positions at random and flip them
    r = randperm(N);
    a(r(1)) = r(2);
    a(r(2)) = r(1);
    MATcostnew = sum(sum(MAT(a,a).*COST));
%    [x,y] = find(MAT(a,a));
%    MATcostnew = sum(profil(abs(x-y)));
    % if this reduced the overall cost
    if (MATcostnew < lowMATcost)
        MAT = MAT(a,a);
        r2 = starta(r(2));
        r1 = starta(r(1));
        starta(r(1)) = r2;
        starta(r(2)) = r1;
        lowMATcost = MATcostnew;
        am = a;
    end;
end;	% h

MATreordered = MATstart(starta,starta) + diagMAT(starta,starta);
MATindices = starta;
MATcost = lowMATcost;

