clear all; close all; clc

Img1 = imread('derek1','jpg');
Img2 = imread('derek2', 'jpg');

Img1 = double(Img1);
Img2 = double(Img2);
[nx ny nz] = size(Img1);


gaussmult = 0.0001;
gaussmult2 = 0.0001;

x = 1:nx; y = 1:ny; [Kx, Ky] = meshgrid(x, y);

for ii = 1:3
    Img1f = fftshift(fft2(Img1(:,:,ii)));
    subplot(3,1,1), pcolor(log(abs(Img1f))), shading interp;
    colormap(hot)
    gauss = exp(-gaussmult*(Kx-nx/2).^2 - gaussmult*(Ky-ny/2).^2);
    Img1f = Img1f.*gauss';
    Imax = max(max(abs(Img1f)));
    Img1cleara(:,:,ii) = ifft2(Img1f);
    Img1f = fftshift(fft2(Img1(:,:,ii)));
    for jj = 1:nx
        for kk = 1:ny
            if abs(Img1f(jj, kk)) < (Imax/800)
                Img1f(jj, kk) = 0;
            end
        end
    end

    Img1clearb(:,:,ii) = ifft2(Img1f);
end

Img1cleara = uint8(abs(Img1cleara));
Img1clearb = uint8(abs(Img1clearb));

Img2f = fftshift(fft2(Img2));
gauss = exp(-gaussmult*(Kx-126.5).^2 - gaussmult*(Ky-180.5).^2);
Img2f = Img2f.*gauss';
Img2clear = (abs((ifft2(Img2f))));

