clear all; close all; clc;

tr_piano=16; % record time in seconds
y=wavread('music1'); 
Fs=length(y)/tr_piano;

v = y'/2;
a = 500;

% L = length(v)/Fs;
% 
% k=(2*pi/(2*L))*[0:length(v)/2-1 -length(v)/2:-1]; ks=fftshift(k);
% 
% tfinal = length(v)/Fs;
% 
% t = (1:length(v))/Fs;
% 
% Sgt_spec = []; tslide = 0:0.05:tfinal;
% figure(2)
% for ii = 1:length(tslide)
%     g = exp(-a*(t-tslide(ii)).^2);
%     Sg = g.*v;
%     Sgt = fft(Sg);
%     Sgt_spec = [Sgt_spec; abs(fftshift(Sgt))];
%     subplot(3,1,1), plot((1:length(v))/Fs,v, (1:length(v))/Fs, g, 'r');
%     axis([0 length(v)/Fs -0.5 1]);
%     subplot(3,1,2), plot((1:length(v))/Fs, Sg);
%     axis([0 length(v)/Fs -0.5 1]);
%     subplot(3,1,3), plot(ks, abs(fftshift(Sgt)));
%     pause(0.5)
%  
%     
% end
% 
% figure(5)
% pcolor(tslide,ks/(2*pi),Sgt_spec.'), shading interp 
% set(gca,'Ylim', [100 600],'Fontsize',[14]) 
% 
