%% Householder
function [final_Q,R] = Householder_YJ(A)
n=size(A,1);
W=zeros(n,n);
E=eye(n);
R=A;
for k = 1:n
    x=R(k:end,k);
    v=sign(x(1))*norm(x)*E(k:end,k)+x;
    v=v/norm(v);
    W(k:end,k)=v;
    R(k:end,k:end)=R(k:end,k:end)-2*v*(v'*R(k:end,k:end));    
end

final_Q=eye(n);
for k=1:n    
    v=W(k:end,k); % get the v-vector
    F=eye(size(v,1))-2*v*v'; % form the F matrix
    Q=eye(n); % initialize the Q matrix
    Q(k:end,k:end)=F; % update the lower F part in it
    final_Q=Q*final_Q;   %pre multiply Q_n*.....*Q_2*Q_1*A=R 
end
final_Q=inv(final_Q); % invert (Q_n*.....*Q_2*Q_1)
end