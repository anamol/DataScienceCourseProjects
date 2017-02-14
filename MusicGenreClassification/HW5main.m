clear all; close all; clc;
% x = waveread('Deodato')
y = wavread('Floyd');
z = wavread('Coltrane');
deodato_rec = 75;
floyd_rec = 75;
y = single(y); z = single(z); 

v = y'/2;
w = z'/2;

Fs = length(v)/floyd_rec;
Fs2 = length(w)/deodato_rec;
% plot((1:length(v))/Fs,v);
%p8 = audioplayer(v,Fs); playblocking(p8);
a = 100;

L = length(v)/Fs;
L2 = length(w)/Fs2;
k=(2*pi/(2*L))*[0:(length(v)-1)/2 -(length(v)-1)/2:-1]; ks=fftshift(k);

if rem(length(k), 2) > 0
    
    ks = ks(length(ks)/2:end);
else
    ks = ks(length(ks)/2-1:end);
end
k = single(k); ks = single(ks);
tfinal = length(v)/Fs;
t = single(1:length(v))/Fs;

Sgt_spec = []; tslide = 0:0.1:tfinal; tslide = single(tslide);
Sgt_spec2 = []; Spec1 = []; Spec2 = [];
% figure(2)
for ii = 1:length(tslide)
    g = exp(-a*(t-tslide(ii)).^2);

    %g = (1 - b * (t - tslide(ii)).^2) .* exp(-a*(t-tslide(ii)).^2);
    Sg = g.*v; Sg2 = g.*w;
    Sgt = fft(Sg); Sgt2 = fft(Sg2);
    Sgt = fftshift(Sgt); Sgt2 = fftshift(Sgt2);
    LSgt = length(Sgt);
    Sgt = Sgt(int64(LSgt/2):end); Sgt2 = Sgt2(int64(LSgt/2):end);
    Sgt = single(Sgt); Sgt2 = single(Sgt2);
    %Sgt_spec = [Sgt_spec; abs((Sgt))]; Sgt_spec2 = [Sgt_spec2; abs((Sgt2))];
    Sgt_spec = [Sgt_spec abs((Sgt))]; Sgt_spec2 = [Sgt_spec2 abs((Sgt2))];
    if rem(ii+50, 50) == 0
%         Sgt_spec = reshape(Sgt_spec, 1, length(Sgt)*50);
%         Sgt_spec2 = reshape(Sgt_spec2, 1, length(Sgt2)*50);
        Spec1 = [Spec1; Sgt_spec];
        Spec2 = [Spec2; Sgt_spec2];
        Sgt_spec = [];
        Sgt_spec2 = [];
    end
    
%     subplot(3,1,1), plot((1:length(v))/Fs,v, (1:length(v))/Fs, g, 'r');
%     xlabel('Time (sec)'); ylabel('Amplitude');
%     set(gca, 'Fontsize', 14);
%     subplot(3,1,2), plot((1:length(v))/Fs, Sg);
%     xlabel('Time (sec)'); ylabel('Amplitude');
%     set(gca, 'Fontsize', 14);
%     %subplot(3,1,3), plot(ks/(2*pi), abs((Sgt)));
%     subplot(3,1,3), plot( abs((Sgt)));
%     xlabel('Frequency (hz)'); ylabel('Amplitude');
%     set(gca, 'Fontsize', 14);
%     pause(0.3);
    
end

'done!'
%%

% figure(5)
% pcolor(tslide,ks/(2*pi),Sgt_spec2.'), shading interp 
% set(gca,'Ylim',[0 1200],'Fontsize',[20])
% xlabel('Time(sec)'); ylabel('Frequency(Hz)');
% 
% clear g; clear v; clear y; clear k;

% Sgt1 = reshape(Sgt_spec, 1, length(Sgt)*length(tslide));
% Sgt2 = reshape(Sgt_spec2, 1, length(Sgt2)*length(tslide));
% 
Sgtot = [Spec1; Spec2];
clear Sgt_spec; clear Sgt_spec2; clear Sgt1; clear Sgt2;
[u,s,v] = svd(Sgtot',0);
%%

rowbeg = 1;
rowend = 5;
trainp = 10;
testp = 15 - trainp;
allanswers = [ones(15,1); 2*ones(15,1)];
Efinal = 0;
for jj = 1:1000
    floydtrain = []; coltrain = [];
    ftest = []; coltest = [];
    fsampl = randsample(15,trainp); ft = setdiff(1:15, fsampl)';
    colsampl = randsample(15,trainp) + 15; colt = setdiff(16:30, colsampl)';
    for kk = 1:trainp
        ff = v(fsampl(kk), rowbeg:rowend);
        cc = v(colsampl(kk), rowbeg:rowend);
        floydtrain = [floydtrain ; ff];
        coltrain = [coltrain ; cc]; 
    end
    for kk = 1:testp
        ff = v(ft(kk), rowbeg:rowend);
        cc = v(colt(kk), rowbeg:rowend);
        ftest = [ftest ; ff];
        coltest = [coltest ; cc]; 
    end
    train = [floydtrain; coltrain];
    tests = [ftest; coltest];
%     floydtrain = v(1:7,rowbeg:rowend); deodatotrain = v(17:23,rowbeg:rowend);
%     train = [floydtrain; deodatotrain];
%     testF = v(11:15,rowbeg:rowend); testD = v(26:30,rowbeg:rowend);
%     tests = [testF; testD];
    correct = [allanswers(1:testp); allanswers(16:(16 + testp-1))];
    answer = [ones(trainp,1); 2*ones(trainp,1)];
%     
%     figure(1)
    
    ind = classify(tests, train, answer);
    Errlda = 0; Errnb = 0;
    for ii = 1:length(ind)
        if ind(ii) == correct(ii)
            Errlda = Errlda + 1;
        end
        
    end
    Errlda = Errlda/length(ind)*100;
    Efinal = Efinal + Errlda;
%     
%     figure(2)
%     nb = fitNaiveBayes(train, answer);
%     prednb = nb.predict(tests);
%     bar(prednb);
%     
%     for ii = 1:length(ind)
%         if ind(ii) == correct(ii)
%             Errnb = Errnb + 1;
%         end
%         
%     end
%     Errnb = Errnb/length(ind)*100
end
Efinal = Efinal/jj;
'After' 
jj 
'trials, the correct % is' 
Efinal