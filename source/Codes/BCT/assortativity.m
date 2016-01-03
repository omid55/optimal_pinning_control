function   r = assortativity(CIJ,flag)

% inputs:
%               CIJ       connection matrix
%               flag      1 = directed graph; 0 = non-directed graph
% outputs:
%               r         assortativity
%
% assortativity, computed after Newman (2002)
%
% Note: Weights are discarded, no edges on main diagonal
%
% Olaf Sporns, Indiana University, 2007/2008
% Vassilis Tsiaras, University of Crete, 2009
%
% =========================================================================

if (flag==0)
    [deg] = degrees_und(CIJ);
    [i,j] = find(triu(CIJ,1)>0);
    K = length(i);
    for k=1:K
        degi(k) = deg(i(k));
        degj(k) = deg(j(k));
    end;
end;
if (flag==1)
    [id,od,deg] = degrees_dir(CIJ);
    [i,j] = find(CIJ>0);
    K = length(i);
    for k=1:K
        degi(k) = deg(i(k));
        degj(k) = deg(j(k));
    end;
end;

% compute assortativity
r = (sum(degi.*degj)/K - (sum(0.5*(degi+degj))/K)^2)/(sum(0.5*(degi.^2+degj.^2))/K - (sum(0.5*(degi+degj))/K)^2);

