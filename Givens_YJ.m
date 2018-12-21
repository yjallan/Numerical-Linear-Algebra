% This is special case of Givens for the type of matrix in Q4
function [Q,R] = Givens_YJ(A)
n = size(A,1);
Q = eye(n);
R = A;
% QR FACTORIZATION USING GIVENS
for j = 1:n-1
    if j==1   % For 1st column, do rotations for entire column
        range_i=n:-1:(j+1);
    else % For all other columns, just do rotation for just 1 element below diagonal
        range_i=j+1;        
    end
for i = range_i
  G = eye(n);  
  % GIVENS Rotation- Finding c and s
  a=R(i-1,j);
  b=R(i,j);
  if b == 0 % if already zero, then just multiply by identity
    c = 1;
    s = 0;
  else
    if abs(b) > abs(a)
      r = a / b;
      s = 1 / sqrt(1 + r^2); 
      c = s*r;
    else
      r = b / a;
      c = 1 / sqrt(1 + r^2); 
      s = c*r;
    end
  end  
  G([i-1, i],[i-1, i]) = [c -s; s c];
  R = G'*R;
  Q = Q*G;
end
end
end