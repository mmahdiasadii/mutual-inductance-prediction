function MI = Mutual_Inductance (n1, n2, nt1, nt2, s1, s2, ri1, ri2,...
    phi, tetha, psi, l, w, h)
r1=zeros(nt1,1);
r2=zeros(nt2,1);
for i=0:nt1-1
    r1(i+1)=ri1+i*s1/cos(pi/n1);
end
for i=0:nt2-1
    r2(i+1)=ri2+i*s2/cos(pi/n2);
end
r1 = flipud(r1); r2 = flipud(r2);
A1=cell(1,nt1);
A2=cell(1,nt2);
for i=1:nt1
    for m=1:n1
        A1{1,i}(1,m)=r1(i)*cos(2*(m-1)*pi/n1);
        A1{1,i}(2,m)=r1(i)*sin(2*(m-1)*pi/n1);
        A1{1,i}(3,m)=0;
    end
end
for j=1:nt2
    for k=1:n2
        A2{1,j}(1,k)=r2(j)*cos(2*(k-1)*pi/n2);
        A2{1,j}(2,k)=r2(j)*sin(2*(k-1)*pi/n2);
        A2{1,j}(3,k)=0;
    end
end
%---------------------------
%Transport Vector
T=[l;w;h];

%Rotation Matrix
Rx=[1 0 0;0 cos(phi) -sin(phi);0 sin(phi) cos(phi)];
Ry=[cos(tetha) 0 sin(tetha);0 1 0;-sin(tetha) 0 cos(tetha)];
Rz=[cos(psi) -sin(psi) 0;sin(psi) cos(psi) 0;0 0 1];
Rt=Rz*Ry*Rx;
rA2=cell(1,nt2);

for j=1:nt2
    for k=1:n2
        rA2{1,j}=[rA2{1,j} Rt*A2{1,j}(:,k)+T];
    end
end
MI=0;
int=cell(i,j);
% syms t1 t2
for i=1:nt1
    for j=1:nt2
        for m=1:n1
            for k=1:n2
                    a1=A1{1,i}(1,m);
                    a2=A1{1,i}(2,m);
                    a3=A1{1,i}(3,m);
                if m==n1
                    b1=A1{1,i}(1,1);
                    b2=A1{1,i}(2,1);
                    b3=A1{1,i}(3,1);
                else
                    b1=A1{1,i}(1,m+1);
                    b2=A1{1,i}(2,m+1);
                    b3=A1{1,i}(3,m+1);
                end
                    c1=rA2{1,j}(1,k);
                    c2=rA2{1,j}(2,k);
                    c3=rA2{1,j}(3,k);
                if k==n2
                    d1=rA2{1,j}(1,1);
                    d2=rA2{1,j}(2,1);
                    d3=rA2{1,j}(3,1);
                else
                    d1=rA2{1,j}(1,k+1);
                    d2=rA2{1,j}(2,k+1);
                    d3=rA2{1,j}(3,k+1);
                end
                u = [b1-a1; b2-a2; b3-a3];
                v = [d1-c1; d2-c2; d3-c3];
%                 x1 = u(1)*t1+a1;
%                 y1 = u(2)*t1+a2;
%                 z1 = u(3)*t1+a3;
%                 x2 = v(1)*t2+c1;
%                 y2 = v(2)*t2+c2;
%                 z2 = v(3)*t2+c3;
                x1 = @(t1) u(1)*t1+a1;
                y1 = @(t1) u(2)*t1+a2;
                z1 = @(t1) u(3)*t1+a3;
                x2 = @(t2) v(1)*t2+c1;
                y2 = @(t2) v(2)*t2+c2;
                z2 = @(t2) v(3)*t2+c3;
                dl1dl2 = u(1).*v(1)+u(2).*v(2)+u(3).*v(3);
                if dl1dl2 == 0
                    PMI = 0;
                    int{i,j}(m,k)=PMI;
                    MI=MI+PMI;
                else
%                     R = sqrt((x2-x1).^2+(y2-y1).^2+(z2-z1).^2);
%                     fun = dl1dl2./R;
%                     f = matlabFunction(fun);
                R = @(t1, t2) sqrt((x2(t2)-x1(t1)).^2+...
                    (y2(t2)-y1(t1)).^2+...
                    (z2(t2)-z1(t1)).^2);
                f = @(t1, t2) dl1dl2./R(t1, t2);
                PMI=integral2(f, 0, 1, 0, 1);
                int{i,j}(m,k)=PMI;
                MI=MI+PMI;
                end
            end
        end
    end
end
MI=MI*1e+2;
end