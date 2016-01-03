function [t,X]=NetSimulator_Rossler(a,b,c,w,OM,Tf,DT,cv,IC)

% [t,X]=NetSimulator_Rossler(a,b,c,R,OM,Tf,DT,cv,IC)
%
% Simulate a Network of mutually coupled Rossler oscillators. The oscillators
% are coupled on the state variable cv.
%
% INPUTS
%  a,b,c,w - standard Rossler parameter vectors, {a,b,c}(i) is {a,b,c} of oscillator i
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
% a=ones(no,1)*0.4;   a=a+a.*0.01.*(2*rand(size(a))-1);
% b=ones(no,1)*0.5;   b=b+b.*0.01.*(2*rand(size(b))-1);
% c=ones(no,1)*5.0;   c=c+c.*0.01.*(2*rand(size(c))-1);
% w=ones(no,1)*1.0;   w=w+w.*0.01.*(2*rand(size(w))-1);
% OM=0.5*diag(rand(1,no-1),-1); 
% Tf=500; DT=5e-2; cv=1; Xo=randn(1,3*no);
% [t,X]=NetSimulator_Rossler(a,b,c,w,OM,Tf,DT,cv,Xo);
% figure(1); for n=1:3, subplot(3,1,n); plot(t,X(:,n:3:end)); end;
% figure(2); plot(X(:,1:3:end),X(:,2:3:end));
% 
%
%
% Authors: O. De Feo
% E-mail: Oscar.DeFeo@epfl.ch
% Date: Feb 2004

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% NOTE: THE WAITBAR HAS BEEN COMMENTED OUT SINCE IT REPRESENTS THE 300% OF
% THE COMPUTING TIME
%
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
  B=sparse(3*no,no);
  BB=sparse(3*no,1);
  
  % Fill B
  stp=3*no+3;
  B(0*3*no+3:stp:end)=+1;

  % Fill BB
  stp=3;
  BB(0*3*no+3:stp:end)=b;
  
  % Fill the diagonal blocks of A
  stp=9*no+3;
  A(0*3*no+2:stp:end)=w;
  A(1*3*no+1:stp:end)=-w;
  A(1*3*no+2:stp:end)=a;
  A(2*3*no+1:stp:end)=-1;
  A(2*3*no+3:stp:end)=-c;
  
  % Correct the coupling
  OM(1:size(OM,1)+1:end)=0;
  OM(1:size(OM,1)+1:end)=-sum(OM,2);
  
  % Assign the coupling in A
  A(cv:3:end,cv:3:end)=A(cv:3:end,cv:3:end)+OM;
  
  % Adjust the orientation for row state
  A=A';
  B=B';
  BB=BB';
  
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
    u=Xnsc(1,[1:3:end]).*Xnsc(1,[3:3:end]);
    Xdot(1,:)=Xnsc(1,:)*A+u*B+BB;
    Xnsc(2,:)=Xnsc(1,:)+Xdot(1,:)*DT;
    % F(T(n),X(n-1)+DT*F(T(n-1),X(n-1))) (Half implicit step)
    u=Xnsc(2,[1:3:end]).*Xnsc(2,[3:3:end]);
    Xdot(2,:)=Xnsc(2,:)*A+u*B+BB;
    % Simulation Step
    X(n,:)=Xnsc(1,:)+(Xdot(1,:)+Xdot(2,:))*DT/2;
    
    % Update waitbar
    % waitbar(n/np,h);
end;
    
% Close waitbar
% close(h);

% End
return;

