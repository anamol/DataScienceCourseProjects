clear all; close all; clc;
load handel

a = 10;
ShaWid = 0.25;
b = 100;

v = y'/2;
L = length(v)/Fs;
k=(2*pi/(2*L))*[0:(length(v)-1)/2 -(length(v)-1)/2:-1]; ks=fftshift(k);

tfinal = length(v)/Fs;
t = (1:length(v))/Fs;

Sgt_spec = []; tslide = 0:0.1:tfinal;


figure(2)
for ii = 1:length(tslide)
    %g = exp(-a*(t-tslide(ii)).^10);
%     if tslide(ii) + ShaWid > 8.92
%         break;
%     end
%     g = 0 * [1:length(v)];
%     if tslide(ii) < 2 * ShaWid
%         tend = floor((tslide(ii) + ShaWid) * 8192); 
%         g(1:tend) = 1;
%    
%     elseif tfinal - tslide(ii) < 2 * ShaWid
%         tbeg = ceil((tslide(ii) * 8192)) + 1;
%         g(tbeg:end) = 1;
%     else
%         tbeg = floor((tslide(ii) - ShaWid) * 8192);
%         tend = floor((tslide(ii) + ShaWid) * 8192);
%         g(tbeg:tend) = 1;
%     end
    g = (1 - b * (t - tslide(ii)).^2) .* exp(-a*(t-tslide(ii)).^2);
    Sg = g.*v;
    Sgt = fft(Sg);
    Sgt_spec = [Sgt_spec; abs(fftshift(Sgt))];
    subplot(3,1,1), plot((1:length(v))/Fs,v, (1:length(v))/Fs, g, 'r');
    xlabel('Time (sec)'); ylabel('Amplitude');
    set(gca, 'Fontsize', 14);
    subplot(3,1,2), plot((1:length(v))/Fs, Sg);
    xlabel('Time (sec)'); ylabel('Amplitude');
    set(gca, 'Fontsize', 14);
    subplot(3,1,3), plot(ks/(2*pi), abs(fftshift(Sgt)));
    xlabel('Frequency (hz)'); ylabel('Amplitude');
    set(gca, 'Fontsize', 14);
    pause(0.3)
    
end

figure(5)
pcolor(tslide,ks/(2*pi),Sgt_spec.'), shading interp 
set(gca,'Ylim',[0 800],'Fontsize',[20])
xlabel('Time(sec)'); ylabel('Frequency(Hz)');

