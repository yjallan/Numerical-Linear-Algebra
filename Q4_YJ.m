clc
clear all
close all
T1=[];
N=[];

n=100;
while n<=1600
    disp(n);
    A=zeros(n); % create a matrix of zeros
    A(:,1)=randn(n,1); % random numbers for column 1
    A(1,:)=randn(1,n); % random numbers for row 1
    A(1:n+1:end) = diag(randn(n)); % random numbers for diagonal   
    
    % Givens Algorithm
    t=cputime;
    [Q,R]=Givens_YJ(A); % calling the function
    T1=[T1,cputime-t];   
       
    N=[N,n];
    n=n*2; % n=100, 200, 400, 800, 1600
end

hold on
f=figure('units','normalized','outerposition',[0 0 1 1]);
plot(N,T1,'-');
set(findall(gca, 'Type', 'Line'),'LineWidth',6);
set(findall(gcf,'-property','FontSize'),'FontSize',18)
xlabel('Dimension n-->');
ylabel('CPU Time-->');
title('Givens Rotation CPU Time vs Dimension');
saveas(f,'CPUTime_Givens.jpg');
hold off
close all


%% SCRATCH
clc
clear all
close all
n=400;

A=zeros(n); % create a matrix of zeros
A(:,1)=randn(n,1); % random numbers for column 1
A(1,:)=randn(1,n); % random numbers for row 1
A(1:n+1:end) = diag(randn(n)); % random numbers for diagonal   

n = size(A,1);
Q = eye(n);
R = A;

% QR FACTORIZATION USING GIVENS
for j = 1:n-1
    if j==1   % For 1st column, do rotations for for entire column
        range_i=n:-1:(j+1);
    else % For all other columns, just do rotation for just 1 element below diagonal
        range_i=j+1;        
    end
for i = range_i
  G = eye(n);  
  % GIVENS Rotation- Finding c and s
  a=R(i-1,j);
  b=R(i,j);
  if b == 0 % if alreaady zero, then just multiply by identity
    c = 1;
    s = 0;
  else
    if abs(b) > abs(a)
      r = a / b;
      s = 1 / sqrt(1 + r^2); % sin is positive
      c = s*r;
    else
      r = b / a;
      c = 1 / sqrt(1 + r^2); % cos is positive
      s = c*r;
    end
  end  
  G([i-1, i],[i-1, i]) = [c -s; s c];
  R = G'*R;
  Q = Q*G;
end
end

  
%% TEST OF ORTHOGONALITY

clc
clear all
close all

n=10;
A=zeros(n); % create a matrix of zeros
A(:,1)=randn(n,1); % random numbers for column 1
A(1,:)=randn(1,n); % random numbers for row 1
A(1:n+1:end) = diag(randn(n)); % random numbers for diagonal   

[Q,R]=Givens_YJ(A); % calling the function
    
T = Q; % assign whichever matrix you want to test.
T_inv = inv(T);
T_tran = T';

tolerance=10^(-9);

error=abs(T_inv-T_tran);
if sum(error(:))<tolerance
   disp('Unitary')
else
   disp('NOT Unitary')
end

%% check to see if the R matrix is indeed Upper triangular
test=R;
n=size(test,1);
tol=10^(-6);
er=0;
for j=1:n
    for i=j+1:n
        er=er+abs(test(i,j));        
    end
end
if er<tol
   disp('matrix is Upper triangular')
else
   disp('matrix is NOT Upper triangular')
end

%% Check to see if multiplying Q * R gives same A matrix
tol=10^(-6);
er=abs(A-Q*R);

if sum(er(:))<tol
   disp('Same')
else
   disp('Different')
end