function E = eccentricity(G)
%E = eccentricity(G); eccentricity EBC, for a binary graph G
%Mahdi Jalili, Sharif, 2011.

SP = distance_bin(full(G)); SP(isfinite(SP) == 0) = 0;
E = max(SP.*(SP~=Inf),[],2);