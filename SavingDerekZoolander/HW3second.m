clear all; close all; clc;

Img1 = imread('derek3', 'jpg');
Img2 = imread('derek4', 'jpg');

gaussmult = 0.0001;

Rash1 = Img1(135:165, 155:185, :);
Rash2 = Img2(135:165, 155:185);
Rash1 = double(Rash1);
Rash2 = double(Rash2);

[nx, ny, nz] = size(Rash1);
x = linspace(0,1,nx); y = linspace(0,1,ny);
dx = x(2)-x(1); dy = y(2)-y(1);

el = ones(nx,1);
Ax = spdiags([el -2*el el],[-1 0 1],nx,nx)/(dx^2); Ix = eye(nx);
Ay = spdiags([el -2*el el],[-1 0 1],ny,ny)/(dy^2); Iy = eye(ny);
L = kron(Iy, Ax) + kron(Ay, Ix);

tspan = [0 0.002 0.004 0.01]; D = 0.1;

Rash1r = reshape(Rash1, nx*ny, 3);

for ii = 1:3
    [t usol1] = ode113('RHSfunc', tspan, Rash1r(:, ii),[],L, D);
    u1(:,:,ii) = usol1;

end

for ii = 1:3
    
    U1(:,:,ii) = reshape(u1(4,:,ii), nx, ny);
end

Img1(140:160, 160:180, :) = U1(6:26, 6:26, :);

Rash2r = reshape(Rash2, nx*ny, 1);
[t usol2] = ode113('RHSfunc', tspan, Rash2r,[],L, D);
U2 = reshape(usol2(4,:),nx, ny);
Img2(140:160, 160:180, :) = U2(6:26, 6:26, :);

x = 1:253; y = 1:361; [Kx, Ky] = meshgrid(x, y);
gauss = exp(-gaussmult*(Kx-126.5).^2 - gaussmult*(Ky-180.5).^2);
Img2f = fftshift(fft2(Img2));
Img2f = Img2f.*gauss';
Img2clear = (abs((ifft2(Img2f))));




