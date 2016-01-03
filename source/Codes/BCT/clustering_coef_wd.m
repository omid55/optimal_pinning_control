function C=clustering_coef_wd(W)
%C=clustering_coef_wd(W); clustering coefficient C for weighted directed graph W.
%
%Reference: Fagiolo, 2007, Phys Rev E 76:026107
%(also see Onnela et al. 2005, Phys Rev E 71:065103);
%
%Mika Rubinov, UNSW, 2007 (last modified July 2008)

%See comments for clustering_coef_bd
%The weighted modification is as follows:
%- The numerator: adjacency matrix is replaced with weights matrix ^ 1/3
%- The denominator: no changes from the binary version
%
%The above reduces to symmetric and/or binary versions of the
%   clustering coefficient for respective graphs.

A=W~=0;                     %adjacency matrix
S=W.^(1/3)+(W.').^(1/3);	%symmetrized weights matrix ^1/3
K=sum(A+A.',2);            	%total degree (in + out)
cyc3=diag(S^3)/2;           %number of 3-cycles (ie. directed triangles)
K(cyc3==0)=inf;             %if no 3-cycles exist, make C=0 (via K=inf)
CYC3=K.*(K-1)-2*diag(A^2);	%number of all possible 3-cycles
C=cyc3./CYC3;               %clustering coefficient

