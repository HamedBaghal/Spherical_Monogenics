clear;
clc;
k = 2;

theta = linspace(0, 2*pi, 50);
phi = linspace(0, pi, 50);
[theta, phi] = meshgrid(theta, phi);

x = sin(phi) .* cos(theta);
y = sin(phi) .* sin(theta);
z = cos(phi);

A = cat(3, x, y, z);
c = [0.540708000691552  -0.251635411477915  -0.802692019194464];
Psi = quaternion(zeros(size(x)), zeros(size(x)), zeros(size(x)), zeros(size(x)));
for i = 1:size(A, 1)
    for j = 1:size(A, 2)
        dot_product = dot(squeeze(A(i, j, :)), c);
        Psi(i, j) = quaternion( ...
            3 * gegenbauerC(k, 1/2, dot_product), ...
            gegenbauerC(k-1, 3/2, dot_product) * (x(i, j)*c(2) - y(i, j)*c(1)), ...
            gegenbauerC(k-1, 3/2, dot_product) * (x(i, j)*c(3) - z(i, j)*c(1)), ...
            gegenbauerC(k-1, 3/2, dot_product) * (y(i, j)*c(3) - z(i, j)*c(2)) );
    end
end

RGB = cat(3, Psi.I, Psi.J, Psi.K);
RGB = rescale(RGB, 0, 1);

figure;
surf(x, y, z, RGB, 'FaceColor', 'interp', 'EdgeColor', 'interp');
axis square tight on;
set(gca, 'box', 'off');
xlabel('x'); ylabel('y'); zlabel('z');

title( ...
  ['Imaginary Part of Spherical Zonal Monogenics $k=2$, $K_{2}(x,\eta_{3})$, ' ...
   'the red point $\eta_{3}=(0.540708  -0.25163  -0.80269)$'], ...
  'Interpreter','latex', 'FontSize',20);

hold on
scatter3(c(1), c(2), c(3), 80, 'r', 'filled');

% ---- label eta_1 on the sphere ----
text(c(1), c(2), c(3), ' $\eta_{3}$', ...
     'Interpreter','latex', ...
     'FontSize',20, ...
     'VerticalAlignment','bottom');



AbsPsi = scalar(Psi);

figure;
surf(x, y, z, AbsPsi, 'FaceColor','interp', 'EdgeColor','none');
axis square tight on;
colormap(gray);
colorbar;

title( ...
  ['Scalar Part of Spherical Zonal Monogenics $k=2$, $K_{2}(x,\eta_{2})$, ' ...
   'the red point $\eta_{2}=(0.540708  -0.25163  -0.80269)$'], ...
  'Interpreter','latex', 'FontSize',20);
axis square tight on;
set(gca, 'box', 'off');
xlabel('x'); ylabel('y'); zlabel('z');

%%%%%%%%%%%
%%%%%%%%%%%
%%%%%%%%%%%
%%%%%%%%%%%

k = 2;

theta = linspace(0, 2*pi, 50);
phi = linspace(0, pi, 50);
[theta, phi] = meshgrid(theta, phi);

x = sin(phi) .* cos(theta);
y = sin(phi) .* sin(theta);
z = cos(phi);

A = cat(3, x, y, z);
a = [0.440661305911686 , -0.115546017429835 ,  0.890206004994525];

Psi1 = quaternion(zeros(size(x)), zeros(size(x)), zeros(size(x)), zeros(size(x)));

for i = 1:size(A, 1)
    for j = 1:size(A, 2)
        dot_product = dot(squeeze(A(i, j, :)), a);
        Psi1(i, j) = quaternion(0.1933,0.1599,-0.2363,-0.261) * ...
            quaternion( ...
            3 * gegenbauerC(k, 1/2, dot_product), ...
            gegenbauerC(k-1, 3/2, dot_product) * (x(i, j)*a(2) - y(i, j)*a(1)), ...
            gegenbauerC(k-1, 3/2, dot_product) * (x(i, j)*a(3) - z(i, j)*a(1)), ...
            gegenbauerC(k-1, 3/2, dot_product) * (y(i, j)*a(3) - z(i, j)*a(2)) );
    end
end

b = [-0.332170702395980  -0.752137836517007   0.569167198061586];
Psi2 = quaternion(zeros(size(x)), zeros(size(x)), zeros(size(x)), zeros(size(x)));

for i = 1:size(x,1)
    for j = 1:size(x,2)
        dot_product = dot([x(i,j),y(i,j),z(i,j)], b);
        Psi2(i,j) = quaternion(0.9123,0,0,0) * ...
            quaternion( ...
            3 * gegenbauerC(k, 1/2, dot_product), ...
            gegenbauerC(k-1, 3/2, dot_product) * (x(i,j)*b(2) - y(i,j)*b(1)), ...
            gegenbauerC(k-1, 3/2, dot_product) * (x(i,j)*b(3) - z(i,j)*b(1)), ...
            gegenbauerC(k-1, 3/2, dot_product) * (y(i,j)*b(3) - z(i,j)*b(2)) );
    end
end

c = [0.540708000691552  -0.251635411477915  -0.802692019194464];
Psi3 = quaternion(zeros(size(x)), zeros(size(x)), zeros(size(x)), zeros(size(x)));

for i = 1:size(x,1)
    for j = 1:size(x,2)
        dot_product = dot([x(i,j),y(i,j),z(i,j)], c);
        Psi3(i,j) = quaternion(0.1933,-0.0212,0.01778,-0.3229) * ...
            quaternion( ...
            3 * gegenbauerC(k, 1/2, dot_product), ...
            gegenbauerC(k-1, 3/2, dot_product) * (x(i,j)*c(2) - y(i,j)*c(1)), ...
            gegenbauerC(k-1, 3/2, dot_product) * (x(i,j)*c(3) - z(i,j)*c(1)), ...
            gegenbauerC(k-1, 3/2, dot_product) * (y(i,j)*c(3) - z(i,j)*c(2)) );
    end
end

Psii = Psi1 + Psi2 + Psi3;

RGB = cat(3, Psii.I, Psii.J, Psii.K);
RGB = rescale(RGB, 0, 1);

figure;
surf(x, y, z, RGB, 'FaceColor', 'interp', 'EdgeColor', 'interp');
axis square tight on;
set(gca, 'box', 'off');
xlabel('x'); ylabel('y'); zlabel('z');


title( ...
  ['Imaginary Part of Spherical Pseudo Zonal Monogenics $k=2$, ' ...
   '$Z_{3}(x)=\sum_{j=1}^{3} K_{2}(x,\eta_{j})\, a_{j2}$'], ...
  'Interpreter','latex', 'FontSize',20);


AbsPsi = scalar(Psi);

figure;
surf(x, y, z, AbsPsi, 'FaceColor','interp', 'EdgeColor','none');
axis square tight on;
colormap(gray);
colorbar;


title( ...
  ['Scalar Part of Spherical Pseudo Zonal Monogenics $k=2$, ' ...
   '$Z_{2}(x)=\sum_{j=1}^{3} K_{2}(x,\eta_{j})\, a_{j2}$'], ...
  'Interpreter','latex', 'FontSize',20);

axis square tight on;
set(gca, 'box', 'off');
xlabel('x'); ylabel('y'); zlabel('z');





