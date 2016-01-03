function M=find_motif34(m,n)
%find_motif34; motif description function
%
%Returns all motif isomorphs for a given motif ID and class (3 or 4):
% MOTIF_MATRICES=find_motif34(MOTIF_ID,MOTIF_CLASS)
%
%Returns the motif id for a given motif matrix (e.g. [0 1 0; 0 0 1; 1 0 0])
% MOTIF_ID=find_motif34(MOTIF_MATRIX)
%
%Mika Rubinov, UNSW, 2007 (last modified July 2008)

persistent M3 ID3 M4 ID4

if isscalar(m)
    if n==3
        if isempty(ID3);
            load motif34lib M3 ID3;
        end
        ind=find(ID3==m).';
        M=zeros(3,3,length(ind));
        for i=1:length(ind)
            M(:,:,i)=reshape([0 M3(ind(i),1:3) 0 ...
                M3(ind(i),4:6) 0],3,3);
        end
    elseif n==4
        if isempty(ID4);
            load motif34lib M4 ID4;
        end
        ind=find(ID4==m).';
        M=zeros(4,4,length(ind));
        for i=1:length(ind)
            M(:,:,i)=reshape([0 M4(ind(i),1:4) 0 ...
                M4(ind(i),5:8) 0 M4(ind(i),9:12) 0],4,4);
        end
    end
else
    n=size(m,1);
    M=eval(['find(motif' int2str(n) 'count(m))']);
end