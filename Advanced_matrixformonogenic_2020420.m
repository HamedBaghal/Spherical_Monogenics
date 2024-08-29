% date: 30-11-2020 Hamed Baghal Ghaffari
%I want to calculate the matrix G including the kernels of spherical
%harmonics
%please enter l as degree of homogeniouty 
%clear

close
clear all
clc


%you may change this k (non-negative integers)
% syms e
% syms f
% syms g
%k=2;
%G=zeros(2*k+1,2*k+1);
%G=zeros(2*k+1,2*k+1,'quaternion');
% a=zeros(2*k+1,3);
% while abs(det(G))==0
%     for i=1:2*k+1
%         theta=(pi)*rand(1);
%         phi=2*pi*rand(1);
%         x=sin(theta)*cos(phi);
%         y=sin(theta)*sin(phi); 
%         z=cos(theta);
%         a(i,1)=x;
%         a(i,2)=y;
%         a(i,3)=z;
% %       c(i)=[x;y;z];
%     end
%     for i=1:2*k+1
%         for j=1:2*k+1
%         G(i,j)=(k+1)*gegenbauerC(k,1/2,dot(a(i,:),a(j,:)));
%         end
%     end

k=2;
syms x
Gegenoneovertwo=gegenbauerC(k,1/2,x);
diffGegenoneovertwo(x)=diff(Gegenoneovertwo,x);


Gegenthreeovertwo=gegenbauerC(k-1,3/2,x);
diffGegenthreeovertwo(x)=diff(Gegenthreeovertwo,x);

%a=zeros(3,k+1);
%for i=1:k+1
%x=rand(1,3);
%a(:,i)=x'/norm(x);
%end

%a=[-.962110407959442, -.211417294179026, -.172180982375314;-.0917140620969262, .989002178900472, -.11603112059399;.410386993978472, -.725295785362625, -.552746360308321;.57065831715715, .00126270293388618, .821186635675029;.195841697823354, -.0962555385036475, -.97590004647273;.140281073335604,-.592624038592789,.793169571668127;-.840938560236222,.131538146854321,.524900041713751];
%a=[-.962110407959442, -.211417294179026, -.172180982375314;-.0917140620969262, .989002178900472, -.11603112059399;.410386993978472, -.725295785362625, -.552746360308321;.57065831715715, .00126270293388618, .821186635675029;.195841697823354, -.0962555385036475, -.97590004647273];

% G=zeros(k+1,k+1);
%a=zeros(k+1,3);

% for i=1:k+1
%     theta=(pi)*rand(1);
%     phi=2*pi*rand(1);
%     x=sin(theta)*cos(phi);
%     y=sin(theta)*sin(phi);
%     z=cos(theta);
%     a(i,1)=x;
%     a(i,2)=y;
%     a(i,3)=z;
% %       c(i)=[x;y;z];
% end

%a=[-.962110407959442, -.211417294179026, -.172180982375314;-.0917140620969262, .989002178900472, -.11603112059399;.410386993978472, -.725295785362625, -.552746360308321];

%a=[-0.078201764197916  -0.165470515326948   0.983109349275943; 0.293431722313200  -0.776791238412154   0.557219163585048; 0.572827045157400   0.338704317307587   0.746423848609785];
%
%a=[-0.078201764197916  -0.165470515326948   0.983109349275943;0.293431722313200  -0.776791238412154   0.557219163585048;0.572827045157400   0.338704317307587   0.746423848609785];

%a=[-0.078201764197916  -0.165470515326948   0.983109349275943;0.293431722313200  -0.776791238412154   0.557219163585048];

a=[-0.338697302058063   0.175744308003532  -0.924336559800027;  0.455569105113474   0.573207545549534  -0.681094633800023; 0.583522346782549  -0.250680170199832  -0.772438426720195];
% % for i=1:k+1
% %         for j=1:k+1
% %         G(i,j)=(2*k+1)*gegenbauerC(k,1/2,dot(a(i,:),a(j,:)))*;
% %         end
% % end


% objectiveFun=0;
% 
% for l=1:k+1
% for j=1:k+1
%     if (l<j) 
%     objectiveFun=objectiveFun+(k+1)^2*(gegenbauerC(k,1/2,dot(a(l,:),a(j,:))))^2+(1-(dot(a(l,:),a(j,:)))^2)*(gegenbauerC(k-1,3/2,dot(a(l,:),a(j,:))))^2;
%     end
% end
% end
% 
% objectiveFun

for p=1:5



for n=1:k+1



w=zeros(1,3);
for j=1:k+1
    if ~(j==n)  
    w=w-((k+1)^2*gegenbauerC(k,1/2,dot(a(n,:),a(j,:)))*diffGegenoneovertwo(dot(a(n,:),a(j,:))) ...
        -dot(a(n,:),a(j,:))*(gegenbauerC(k-1,3/2,dot(a(n,:),a(j,:))))^2 ...
        +(1-(dot(a(n,:),a(j,:)))^2)*gegenbauerC(k-1,3/2,dot(a(n,:),a(j,:)))*diffGegenthreeovertwo(dot(a(n,:),a(j,:)))).*(-a(j,:)+dot(a(n,:),a(j,:)).*a(n,:));
    end
end
%end


w=vpa(w/(norm(w)));

syms t
objectiveFun=0;
for j=1:k+1
    for l=1:k+1
        if ~(n==j) && ~(n==l) && (l<j)
            objectiveFun=objectiveFun+(k+1)^2*(gegenbauerC(k,1/2,dot(a(l,:),a(j,:))))^2+(1-(dot(a(l,:),a(j,:)))^2)*(gegenbauerC(k-1,3/2,dot(a(l,:),a(j,:))))^2;
        end
    end
end

for j=1:k+1
    if ~(n==j)
        objectiveFun=objectiveFun+(k+1)^2*(gegenbauerC(k,1/2,dot(cos(t)*a(n,:)+sin(t)*w,a(j,:))))^2+(1-(dot(cos(t)*a(n,:)+sin(t)*w,a(j,:)))^2)*(gegenbauerC(k-1,3/2,dot(cos(t)*a(n,:)+sin(t)*w,a(j,:))))^2;
    end
end

% 
% % Define the objective function as an anonymous function
f = @(t) double(subs(objectiveFun, t));
% 
% % Perform the optimization
[tmin, fval] = fminbnd(f, 0, 4*pi);

a(n,:)=cos(tmin)*a(n,:)+sin(tmin)*w;

objectiveFun=0;

for l=1:k+1
for j=1:k+1
    if (l<j) 
    objectiveFun=objectiveFun+(k+1)^2*(gegenbauerC(k,1/2,dot(a(l,:),a(j,:))))^2+(1-(dot(a(l,:),a(j,:)))^2)*(gegenbauerC(k-1,3/2,dot(a(l,:),a(j,:))))^2;
    end
end
end

objectiveFun;

end

n=randi([1,k+1]);
% 
w=zeros(1,3);
for j=1:k+1
    if ~(j==n)  
    w=w-((k+1)^2*gegenbauerC(k,1/2,dot(a(n,:),a(j,:)))*diffGegenoneovertwo(dot(a(n,:),a(j,:))) ...
        -dot(a(n,:),a(j,:))*(gegenbauerC(k-1,3/2,dot(a(n,:),a(j,:))))^2 ...
        +(1-(dot(a(n,:),a(j,:)))^2)*gegenbauerC(k-1,3/2,dot(a(n,:),a(j,:)))*diffGegenthreeovertwo(dot(a(n,:),a(j,:)))).*(-a(j,:)+dot(a(n,:),a(j,:)).*a(n,:));
    end
end
%end


w=vpa(w/(norm(w)));

syms t
objectiveFun=0;
for j=1:k+1
    for l=1:k+1
        if ~(n==j) && ~(n==l) && (l<j)
            objectiveFun=objectiveFun+(k+1)^2*(gegenbauerC(k,1/2,dot(a(l,:),a(j,:))))^2+(1-(dot(a(l,:),a(j,:)))^2)*(gegenbauerC(k-1,3/2,dot(a(l,:),a(j,:))))^2;
        end
    end
end

for j=1:k+1
    if ~(n==j)
        objectiveFun=objectiveFun+(k+1)^2*(gegenbauerC(k,1/2,dot(cos(t)*a(n,:)+sin(t)*w,a(j,:))))^2+(1-(dot(cos(t)*a(n,:)+sin(t)*w,a(j,:)))^2)*(gegenbauerC(k-1,3/2,dot(cos(t)*a(n,:)+sin(t)*w,a(j,:))))^2;
    end
end
% 
% % Define the objective function as an anonymous function
f = @(t) double(subs(objectiveFun, t));
% 
% % Perform the optimization
[tmin, fval] = fminbnd(f, 0, 4*pi);


a(n,:)=cos(tmin)*a(n,:)+sin(tmin)*w;

objectiveFun=0;

for l=1:k+1
for j=1:k+1
    if (l<j) 
    objectiveFun=objectiveFun+(k+1)^2*(gegenbauerC(k,1/2,dot(a(l,:),a(j,:))))^2+(1-(dot(a(l,:),a(j,:)))^2)*(gegenbauerC(k-1,3/2,dot(a(l,:),a(j,:))))^2;
    end
end
end

objectiveFun;

end
% 
% end
% 
% a
% 

% objectiveFun=0;
% 
% for l=1:k+1
% for j=1:k+1
%     if (l<j) 
%     objectiveFun=objectiveFun+(k+1)^2*(gegenbauerC(k,1/2,dot(a(l,:),a(j,:))))^2+(1-(dot(a(l,:),a(j,:)))^2)*(gegenbauerC(k-1,3/2,dot(a(l,:),a(j,:))))^2;
%     end
% end
% end
% 
% objectiveFun


%The best points 
a=[0.440661305911686  -0.115546017429835   0.890206004994525; -0.332170702395980  -0.752137836517007   0.569167198061586; 0.540708000691552  -0.251635411477915  -0.802692019194464];

objectiveFun=0;

for l=1:k+1
for j=1:k+1
    if (l<j) 
    objectiveFun=objectiveFun+(k+1)^2*(gegenbauerC(k,1/2,dot(a(l,:),a(j,:))))^2+(1-(dot(a(l,:),a(j,:)))^2)*(gegenbauerC(k-1,3/2,dot(a(l,:),a(j,:))))^2;
    end
end
end

objectiveFun;

syms e0 e1 e2 e12 e13 e23
G=sym(zeros(k+1));
for i=1:(k+1)
    G(i,i)=(3.0)*e0;
    for j=2:(k+1)
        if j>i
            G(i,j)=3*gegenbauerC(k,1/2,dot(a(i,:),a(j,:)))*e0+gegenbauerC(k-1,3/2,dot(a(i,:),a(j,:)))*[(a(i,1)*a(j,2)-a(i,2)*a(j,1))*e12 + (a(i,1)*a(j,3)-a(i,3)*a(j,1))*e13+ ...
                    +(a(i,2)*a(j,3)-a(i,3)*a(j,2))*e23];
            G(j,i)=3*gegenbauerC(k,1/2,dot(a(i,:),a(j,:)))*e0+gegenbauerC(k-1,3/2,dot(a(i,:),a(j,:)))*[-(a(i,1)*a(j,2)-a(i,2)*a(j,1))*e12 - (a(i,1)*a(j,3)-a(i,3)*a(j,1))*e13+ ...
                    -(a(i,2)*a(j,3)-a(i,3)*a(j,2))*e23];
        end
    end
end
G=vpa(G,4);

T_mat_of_G=vpa(subs(G, {e12 e13 e23}, {e1 e2 e12}),4);

T_mat_of_G_plus=vpa(subs(T_mat_of_G, {e1 e2}, {0 0}),4);

% 1st Method
T_mat_of_G_minus=vpa(subs(T_mat_of_G, {e0 e12}, {0 0}),4);
T_mat_of_G_tilde_plus=vpa(subs(T_mat_of_G_minus, {e1 e2}, {e0 e12}),4);

% 2nd Method
%T_mat_of_G_minus_times_e1=vpa(T_mat_of_G_minus*e1,4);
%T_mat_of_G_minus_times_e1=vpa(expand(T_mat_of_G_minus*e1),4)
%T_mat_of_G_tilde_plus=vpa(subs(T_mat_of_G_minus_times_e1, {e1*e1 e2*e1}, {1 e12}),4)

phi_T_mat_of_G_plus=vpa(subs(T_mat_of_G_plus, {'e12'}, {-e12}),4);

phi_T_mat_of_G_tilde_plus=vpa(subs(T_mat_of_G_tilde_plus, {'e12'}, {-e12}),4);

Chi_T_mat_of_G=[T_mat_of_G_plus,T_mat_of_G_tilde_plus;-phi_T_mat_of_G_tilde_plus,phi_T_mat_of_G_plus];

Chi_T_mat_of_G=vpa(Chi_T_mat_of_G,4);

T_mat_of_Chi_T_mat_of_G=vpa(subs(Chi_T_mat_of_G, {e0 e12}, {1 1i}),4);

% [V,D]=eig(BB)
% CC=V*D^(-1/2)*inv(V)

B=vpa((T_mat_of_Chi_T_mat_of_G)^(-1/2),4);

T_mat_inverse_of_B=vpa(real(B)*e0+subs(B-real(B), {'1i'}, {e12}),4);
 
[r, c] = size(T_mat_inverse_of_B);

r2 = r / 2;
c2 = c / 2;
  
T_mat_inverse_of_B_plus=vpa(T_mat_inverse_of_B(1 : r2 , 1 : c2),4);
T_mat_inverse_of_B_tilde_plus = T_mat_inverse_of_B(1 : r2 , c2 + 1 : c);
T_mat_inverse_of_B_tilde_plus=vpa(T_mat_inverse_of_B_tilde_plus,4);
T_mat_inverse_of_B_tilde_plus_times_e1=vpa(expand(T_mat_inverse_of_B_tilde_plus*e1),4);
T_mat_inverse_of_B_minus=vpa(subs(T_mat_inverse_of_B_tilde_plus_times_e1, {e0 e1*e12}, {1 e2}),4);


Chi_inverse_T_mat_inverse_of_B=vpa(T_mat_inverse_of_B_minus+T_mat_inverse_of_B_plus,4);

T_mat_inverse_Chi_inverse_T_mat_inverse_of_B=vpa(subs(Chi_inverse_T_mat_inverse_of_B, {'e1' 'e2' 'e12'},{e12 e13 e23}),4)

S_Without_Expansion=(T_mat_inverse_Chi_inverse_T_mat_inverse_of_B*G*T_mat_inverse_Chi_inverse_T_mat_inverse_of_B);

SS=0;
for k=1:3
        SS=SS+G(1,k)*T_mat_inverse_Chi_inverse_T_mat_inverse_of_B(k,1);
end

SS=vpa(expand(SS),4)
SS=vpa(subs(SS,{e0 e12^2 e13^2 e23^2},{1 -1 -1 -1}),4);
SS=vpa(subs(SS,{e12*e13 e12*e23 e13*e23},{e23 -e13 e12}),4);
SS=vpa(subs(SS,{e23^2},{-1}),4)

S=vpa(expand(T_mat_inverse_Chi_inverse_T_mat_inverse_of_B*G*T_mat_inverse_Chi_inverse_T_mat_inverse_of_B),4);

AA=[quaternion(3.0,0,0,0), quaternion(- 0.6,- 0.4962,0.7332,+ 0.8101),quaternion(- 0.6,0.06495,1.12, -0.425);...
quaternion(- 0.6,0.4962,- 0.7332,- 0.8101),quaternion(3.0,0,0,0),quaternion(- 0.6,- 0.6578,0.05517,- 1.002);...
quaternion(- 0.6,- 0.06495,- 1.12,0.425),quaternion(- 0.6,0.6578,- 0.05517,1.002),quaternion(3.0,0,0,0)];
BB=adjoint(AA);
CC=BB^(-1/2);
Final=unadjoint(CC);


Final(1,1);

Final(1,2);
Final(1,3);

Final(2,1);

Final(2,2);
Final(2,3);
Final(3,1);

Final(3,2);
Final(3,3);
% 
% AdointFinal=[quaternion(0.912,0,0,0),quaternion(0.1933,-0.1598,0.2362,0.2609),quaternion(0.1933,0.02091,0.3608,-0.1369); ...
%     quaternion(0.1933,0.1598,-0.2362,-0.2609),quaternion(0.9121,0,0,0),quaternion(0.1932,-0.2119,0.01776,-0.3228);...
%     quaternion(0.1933,-0.02091,-0.3608,0.1369),quaternion(0.1932,0.2119,-0.01776,0.3228),quaternion(0.912,0,0,0)];
% 
% 
% Adjoint_T_mat_inverse_Chi_inverse_T_mat_inverse_of_B=vpa(subs(T_mat_inverse_Chi_inverse_T_mat_inverse_of_B,{e12 e13 e23}, {-e12 -e13 -e23}),4);
% 
% S=0;
% for k=1:3
%     for l=1:3
%         S=S+Adjoint_T_mat_inverse_Chi_inverse_T_mat_inverse_of_B(1,k)*G(k,l)*T_mat_inverse_Chi_inverse_T_mat_inverse_of_B(l,2);
%     end
% end
S=vpa(subs(S,{e0 e12^2 e13^2 e23^2},{1 -1 -1 -1}),4);
S=vpa(subs(S,{e12^3 e13^3 e23^3},{-e12 -e13 -e23}),4);
S=vpa(subs(S,{e12*e13 e12*e23 e13*e23},{e23 -e13 e12}),4);
S=vpa(subs(S,{e23^2},{-1}),4);
% 
% S=vpa(subs(S,{e12^2 e13^2 e23^2},{-1 -1 -1}),4)
% S=vpa(subs(S,{e12^3 e13^3 e23^3},{-e12 -e13 -e23}),4)
% S=vpa(subs(S,{e12*e13 e12*e23 e13*e23},{e23 -e13 e12}),4)
% S=vpa(subs(S,{e23^2},{-1}),4)



a1=vpa(T_mat_inverse_Chi_inverse_T_mat_inverse_of_B(:,1),4)

vpa(G*a1,4)
vpa(expand(G*a1),4)
