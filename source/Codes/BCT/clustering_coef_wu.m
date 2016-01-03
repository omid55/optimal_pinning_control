function C=clustering_coef_wu(W)
%C=clustering_coef_wu(W); clustering coefficient C for weighted undirected graph W.
%
%Reference: Onnela et al. 2005, Phys Rev E 71:065103
%
%Mika Rubinov, UNSW, 2007 (last modified July 2008)

K=sum(W~=0,2);            	
cyc3=diag((W.^(1/3))^3);           
K(cyc3==0)=inf;             %if no 3-cycles exist, make C=0 (via K=inf)
C=cyc3./(K.*(K-1));         %clustering coefficient
