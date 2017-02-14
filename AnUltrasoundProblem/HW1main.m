clear all; close all; clc;
load Testdata

L=15; % spatial domain
n=64; % Fourier modes

% define grid vectors, create 3D cartesian grid
x2=linspace(-L,L,n+1); x=x2(1:n); y=x; z=x;
k=(2*pi/(2*L))*[0:(n/2-1) -n/2:-1]; ks=fftshift(k);

[X,Y,Z]=meshgrid(x,y,z);
[Kx,Ky,Kz]=meshgrid(ks,ks,ks);

Utave = zeros(n,n,n); % variable for averaged frequency 

% calculate average of frequencies over all 20 signals
for j=1:20
Un(:,:,:)=reshape(Undata(j,:),n,n,n);
Utave = Utave + fftn(Un);
end


Utave = abs(fftshift(Utave)); % shift and take abs value of 
                              % averaged frequency data

% find index of center frequency in Utave
maxi = 0;
for ii = 1:64
    for jj = 1:64
        for kk = 1:64
            if Utave(ii, jj, kk) > maxi
                maxi = Utave(ii, jj, kk);
                a = ii; b = jj, c = kk;
            end
        end
    end
end

% isosurface visualization of frequency signature
figure(1)

isosurface(Kx,Ky,Kz,abs(Utave)/max(Utave(:)),0.6, 'r');
set(gca, 'FontSize',18);
axis([ks(1) -ks(1) ks(1) -ks(1) ks(1) -ks(1)]), grid on;
xlabel('Kx'); ylabel('Ky'); zlabel('Kz');

% create Gauss filter
gaussfilter = exp(-((Kx-ks(b)).^2 + (Ky-ks(a)).^2 + (Kz-ks(c)).^2)/2);

% slice visualization of Gauss filter
slice(Kx,Ky, Kz, gaussfilter, 2,-1,0);
set(gca, 'FontSize', 20);
axis([ks(1) -ks(1) ks(1) -ks(1) ks(1) -ks(1)]);
xlabel('Kx'); ylabel('Ky'); zlabel('Kz');

% Apply Gauss filter to all 20 signals
maxi = 0; M = [];
for j = 1:20
    Un(:,:,:)=reshape(Undata(j,:),n,n,n);
    Utn = fftshift(fftn(Un));
    
    Newt = Utn.*gaussfilter;
    New= ifftn(Newt); 
    
    % find index of maximum element in each 3D matrix
    for ii = 1:64
        for jj = 1:64
            for kk = 1:64
                if abs(New(ii,jj,kk)) > maxi
                    maxi = abs(New(ii,jj,kk));
                    bb = X(1,ii,1); aa = Y(jj,1,1); cc = Z(1,1,kk);
                end
            end
        end
    end
    maxi = 0;
    M = [M ;[aa bb cc]]; % store indices of maximum elements 
                         % of all 20 time steps in M
    
    % isosurface visualization of trajectory of marble
    figure(2)
    isosurface(X, Y, Z, abs(New)/max(abs(New(:))), 0.7);
    axis([-L L -L L -L L]), grid on;
    set(gca, 'FontSize',20);
    xlabel('X'); ylabel('Y'); zlabel('Z');
    pause(0.1);
end


% visualization of trajectory of marble using plot3
figure(3)
plot3(M(:,1), M(:,2), M(:,3),'k--o', 'LineWidth', 5);
axis([-L L -L L -L L]), grid on;
set(gca, 'FontSize',20);
xlabel('X'); ylabel('Y'); zlabel('Z');

% isosurface visualization of position of marble in final time step
figure(4)
isosurface(X, Y, Z, abs(New)/max(abs(New(:))), 0.7);
set(gca, 'FontSize',20);
xlabel('X'); ylabel('Y'); zlabel('Z');
axis([-L L -L L -L L]), grid on;