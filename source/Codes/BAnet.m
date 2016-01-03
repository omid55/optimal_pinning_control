function A=BAnet(N, m, B);

% input parameters: N = size of network
% m = degree of new vertices joining the network
% Starting with a fully connected graph with m+a nodes, the new edges are
% craeted with a propobability proportional to k+B, where k is the degree
% of the nodes.

% initialize the adjacency matrix

A=spalloc(N,N,2*2*N+50);

% create a seed network, here only 3 vertices all connected to each other

for i= 1:m
    for j = i+1:m+1
        A(i,j) = 1; A(j,i) = 1;
    end
end

distribution = ones(1, m+B);
for i = 2:m+1
    distribution = [distribution i*ones(1,m+B)];
end
        
                            % each element points to a vertex index, 
                            % # of elements = degree of vertex
                            % so, if an element is picked at random
                            % it points to a vertex of degree k 
                            % with the probability k/sum(k)

dlength=length(distribution);   % current length of the distribution vector = sum of vertex degrees
                                

for i=m+2:N,                % grow the network

    vertexlist=zeros(1,m);  % this vector will contain indices of vertices 
                            % chosen for attachment
    
    found=0;                % how many vertices we have found
    
    while (found<m),        % we want m different vertices
        
        candidate=distribution(ceil(rand(1,1)*dlength));  % pick a vertex with probab. proportional
                                                          % to its degree
                                                          
        if length(find(vertexlist==candidate))==0         % if not already on the list,
            
            found=found+1;                                % one more vertex chosen
            vertexlist(found)=candidate;                  % add it to the list
            
        end;
        
    end;
    
    % now we have m vertices
    % first, connect the new vertex (number i) to the m "parent" vertices
    
    for j=1:m,
        
        A(i,vertexlist(j))=1;
        A(vertexlist(j),i)=1;
        
    end;
    
    % next, add elements to the distribution vector
    % first, the chosen vertices: the degree of each
    % grows by one, so we add one
    
    for j=1:m,
        
        dlength=dlength+1;
        distribution(dlength)=vertexlist(j);
        
    end;
    
    % then add the newly born vertex i to the list
	% its degree is m, so we add m elements	

    for j=1:m+B,
        
        dlength=dlength+1;
        distribution(dlength)=i;
        
    end;
        
end;        % end of the network growth loop

% OK we have our network! Now just calculate the degree distribution.
% First, we'll find the degree of each vertex by summing over the 
% adjacency matrix (the number of entries in row i = # of i's neighbours)

k=full(sum(A,1)); % vector of N elements so that k(i) is the degree of vertex i
                  % note: full() makes the vector non-sparse
                  
kmax=max(k); % the maximum degree

Nk=zeros(1,kmax); % vector for degree counts, kmax elements

for i=1:N,        % go through all vertices
    
    Nk(k(i))=Nk(k(i))+1; % grow the number of vertices of degree k(i) by one
    
end;

Nk=Nk/sum(Nk); % normalize the distribution

% loglog(Nk,'ko'); % and finally plot it!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% next we will plot the theoretical distributions

for i=1:kmax,
    
    pk_asymptotic(i)=2*m*m/(i*i*i);
    pk_exact(i)=2*m*(m+1)/(i*(i+1)*(i+2));
    
end;

% hold on;
% 
% loglog(pk_asymptotic,'k-');
% loglog(pk_exact,'r-');
