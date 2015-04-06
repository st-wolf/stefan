function [ rez ] = prog( M, V )

% INPUT:
% :M: matrix
% :V: right side values vector

% OUTPUT:
% :res: solution vector

A = diag(M, -1);
B = diag(M);
C = diag(M, 1);
D = V;
n = length(B);

Cz=zeros(1,n);
Dz=zeros(1,n);
Cz(1)=C(1)/B(1);
for j=1:n-2
    Cz(j+1)=C(j+1)/(B(j+1)-A(j)*Cz(j));
end
Dz(1)=D(1)/B(1);
for j=1:n-1
    Dz(j+1)=(D(j+1)-A(j)*Dz(j))/(B(j+1)-A(j)*Cz(j));
end

X(n)=Dz(n);
for j=n-1:-1:1
    X(j)=Dz(j)-Cz(j)*X(j+1);
end

rez = X';

end

