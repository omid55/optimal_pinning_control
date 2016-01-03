function [t,X]=NetSimulator_Lorenz(s,r,b,c,OM,Tf,DT,cv,IC)

% [t,X]=NetSimulator_Lorenz(s,r,b,OM,Tf,DT,cv,IC)
%
% Simulate a Network of mutually coupled Lorenz oscillators. The oscillators
% are coupled on the state variable cv.
%
% INPUTS
%  s,r,b - standard Lorenz parameter vectors, {s,r,b}(i) is {s,r,b} of oscillator i
%      c - coupling strength (same for all oscillators)
%     OM - coupling matrix, OM(i,j) coupling from oscillator j to oscillator i
%     Tf - span time of simulation
%     DT - integration fixed step time
%     cv - coupling variable 1,2, or 3
%     IC - Initial condition of the oscillators (1x3*no)
%
% OUTPUTS
%    t - vector of the time
%    X - matrix of the state evolution -> t(i), X(i,:)
%
%
% EXAMPLES
% 
% no=5; 
% s=ones(no,1)*10;   s=s+s.*0.01.*(2*rand(size(s))-1);
% b=ones(no,1)*8/3;  b=b+b.*0.01.*(2*rand(size(b))-1);
% r=ones(no,1)*27;   r=r+r.*0.01.*(2*rand(size(r))-1);
% c=10;
% OM=15*diag(rand(1,no-1),-1); 
% Tf=100; DT=1e-2; cv=1; Xo=randn(1,3*no);
% [t,X]=NetSimulator_Lorenz(s,r,b,c,OM,Tf,DT,cv,Xo);
% figure(1); for n=1:3, subplot(3,1,n); plot(t,X(:,n:3:end)); end;
% figure(2); plot(X(:,1:3:end),X(:,3:3:end));
% 
%
% no=5; 
% s=ones(no,1)*10;   s=s+s.*0.01.*(2*rand(size(s))-1);
% b=ones(no,1)*8/3;  b=b+b.*0.01.*(2*rand(size(b))-1);
% r=ones(no,1)*27;   r=r+r.*0.01.*(2*rand(size(r))-1);
% c=10;
% OM=15*diag(rand(1,no-1),-1); 
% Tf=100; DT=1e-2; cv=2; Xo=randn(1,3*no);
% [t,X]=NetSimulator_Lorenz(s,r,b,c,OM,Tf,DT,cv,Xo);
% figure(1); for n=1:3, subplot(3,1,n); plot(t,X(:,n:3:end)); end;
% figure(2); plot(X(:,1:3:end),X(:,3:3:end));
% 
%
% no=5; 
% s=ones(no,1)*10;   s=s+s.*0.01.*(2*rand(size(s))-1);
% b=ones(no,1)*8/3;  b=b+b.*0.01.*(2*rand(size(b))-1);
% r=ones(no,1)*27;   r=r+r.*0.01.*(2*rand(size(r))-1);
% c=10;
% OM=15*diag(rand(1,no-1),-1); 
% Tf=100; DT=1e-2; cv=3; Xo=randn(1,3*no);
% [t,X]=NetSimulator_Lorenz(s,r,b,c,R,OM,Tf,DT,cv,Xo);
% figure(1); for n=1:3, subplot(3,1,n); plot(t,X(:,n:3:end)); end;
% figure(2); plot(X(:,1:3:end),X(:,3:3:end));
%
%

% Created by Mahdi Jalili based on a code by Oscar De Feo
% E-mail: mahdi.jalili@epfl.ch
% Date: July 2006

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Initialize variables and settings

  % Allocate results
  t=(0:ceil(Tf/DT))'*DT;
  np=size(t,1);
  no=size(OM,1);
  X=zeros(np,3*no);
  
  % Internal variables for simulation
  Xdot=zeros(2,3*no);
  Xnsc=zeros(2,3*no);

  % Lur'e rapresentation of system
  A=sparse(3*no,3*no);
  B=sparse(3*no,2*no);
  
  % Fill B
  stp=6*no+3;
  B(0*3*no+2:stp:end)=-1;
  B(1*3*no+3:stp:end)=+1;
  
  % Fill the diagonal blocks of A
  stp=9*no+3;
  A(0*3*no+1:stp:end)=-s;
  A(0*3*no+2:stp:end)=r;
  A(1*3*no+1:stp:end)=+s;
  A(1*3*no+2:stp:end)=-1;
  A(2*3*no+3:stp:end)=-b;
  
  % Correct the coupling
   OM=full(OM);
   for i=1:no
      OM(i,i)=-sum(OM(i,:));
   end
   OM(1:size(OM,1)+1:end)=0;
   OM(1:size(OM,1)+1:end)=-sum(OM,2);
   
  % Influencing of coupling strength
   OM=c*OM; 
  
  % Assign the coupling in A
  A(cv:3:end,cv:3:end)=A(cv:3:end,cv:3:end)+OM;
  
  % Adjust the orientation for row state
  A=A';
  B=B';
    
  % Initial state
  X(1,:)=IC(:)';
  
  % Initialize waitbar
  % h = waitbar(0,'Please wait...');
  
    
  
% Free run for np steps
for n=2:np,
        
    % Current point Step
    Xnsc(1,:)=X(n-1,:);
    
    % The two derivatives for the Heun method
    % F(T(n-1),X(n-1)) (previous step)
    % u=[Xnsc(1,1:3:end).*Xnsc(1,3:3:end); Xnsc(1,1:3:end).*Xnsc(1,2:3:end)];
    u=Xnsc(1,[1:3:end; 1:3:end]).*Xnsc(1,[3:3:end; 2:3:end]);
    Xdot(1,:)=Xnsc(1,:)*A+u*B;
    Xnsc(2,:)=Xnsc(1,:)+Xdot(1,:)*DT;
    % F(T(n),X(n-1)+DT*F(T(n-1),X(n-1))) (Half implicit step)
    u=Xnsc(2,[1:3:end; 1:3:end]).*Xnsc(2,[3:3:end; 2:3:end]);
    Xdot(2,:)=Xnsc(2,:)*A+u*B;
    % Simulation Step
    X(n,:)=Xnsc(1,:)+(Xdot(1,:)+Xdot(2,:))*DT/2;
    
    % Update waitbar
    % waitbar(n/np,h);
end;
    
% Close waitbar
% close(h);

% End
return;

