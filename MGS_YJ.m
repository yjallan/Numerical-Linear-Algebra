%% Modified gram scmidt
function [Q,R] = MGS_YJ(A)
n=size(A,1);
Q=zeros(n,n);
R=zeros(n,n);
V=A;
for i=1: n    
    R(i,i)=norm(V(:,i));
    Q(:,i)=V(:,i)/R(i,i);
    for j=i+1:n
        R(i,j)=Q(:,i)'*V(:,j);        
        V(:,j)=V(:,j)-R(i,j)*Q(:,i);
    end    
end
end