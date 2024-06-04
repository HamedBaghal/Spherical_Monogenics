clear all
clc
close all
k=2;
a=[0.440661305911686  -0.115546017429835   0.890206004994525; -0.332170702395980  -0.752137836517007   0.569167198061586; 0.540708000691552  -0.251635411477915  -0.802692019194464];

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
G=vpa(G,4)

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

%B=vpa((T_mat_of_Chi_T_mat_of_G)^(-1/2),4);


T_mat_inverse_of_T_mat_of_Chi_T_mat_of_G=vpa(real(T_mat_of_Chi_T_mat_of_G)*e0+subs(T_mat_of_Chi_T_mat_of_G-real(T_mat_of_Chi_T_mat_of_G), {'1i'}, {e12}),4);
%T_mat_inverse_of_B=vpa(real(B)*e0+subs(B-real(B), {'1i'}, {e12}),4);

[r, c] = size(T_mat_inverse_of_T_mat_of_Chi_T_mat_of_G);

r2 = r / 2;
c2 = c / 2;
  
T_mat_inverse_of_T_mat_of_Chi_T_mat_of_G_plus=vpa(T_mat_inverse_of_T_mat_of_Chi_T_mat_of_G(1 : r2 , 1 : c2),4);
T_mat_inverse_of_T_mat_of_Chi_T_mat_of_G_tilde_plus = vpa(T_mat_inverse_of_T_mat_of_Chi_T_mat_of_G(1 : r2 , c2 + 1 : c),4);
% T_mat_inverse_of_B_tilde_plus=vpa(T_mat_inverse_of_B_tilde_plus,4);
T_mat_inverse_of_T_mat_of_Chi_T_mat_of_G_tilde_plus_times_e1=vpa(expand(T_mat_inverse_of_T_mat_of_Chi_T_mat_of_G_tilde_plus*e1),4);
T_mat_inverse_of_T_mat_of_Chi_T_mat_of_G_minus=vpa(subs(T_mat_inverse_of_T_mat_of_Chi_T_mat_of_G_tilde_plus_times_e1, {e0 e1*e12}, {1 e2}),4);
% 
% 
Chi_inverse_T_mat_inverse_of_T_mat_of_Chi_T_mat_of_G=vpa(T_mat_inverse_of_T_mat_of_Chi_T_mat_of_G_minus+T_mat_inverse_of_T_mat_of_Chi_T_mat_of_G_plus,4);
% 
T_mat_inverse_Chi=vpa(subs(Chi_inverse_T_mat_inverse_of_T_mat_of_Chi_T_mat_of_G, {'e1' 'e2' 'e12'},{e12 e13 e23}),4)
% 

