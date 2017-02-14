clear all; close all; clc;

y = wavread('Allthatjazz');
z = wavread('ClassicWrok');
x = wavread('Grunge');

deodato_rec = 75;
floyd_rec = 75;
y = single(y); z = single(z); x = single(x);

v = y'/2;
w = z'/2;
x = x'/2;
Fs = length(v)/floyd_rec;
Fs2 = length(w)/deodato_rec;
Fs3 = length(x)/deodato_rec;
a = 100;

L = length(v)/Fs;
L2 = length(w)/Fs2;
L3 = length(x)/Fs3;
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
Sgt_spec2 = []; Sgt_spec3 = []; Spec1 = []; Spec2 = [];
Spec3 = [];

for ii = 1:length(tslide)
    g = exp(-a*(t-tslide(ii)).^2);
    Sg = g.*v; Sg2 = g.*w; Sg3 = g.*x;
    Sgt = fft(Sg); Sgt2 = fft(Sg2); Sgt3 = fft(Sg3);
    Sgt = fftshift(Sgt); Sgt2 = fftshift(Sgt2); Sgt3 = fftshift(Sgt3);
    LSgt = length(Sgt);
    Sgt = Sgt(int64(LSgt/2):end); Sgt2 = Sgt2(int64(LSgt/2):end);
    Sgt3 = Sgt3(int64(LSgt/2):end);
    Sgt = single(Sgt); Sgt2 = single(Sgt2); Sgt3 = single(Sgt3);

    Sgt_spec = [Sgt_spec abs((Sgt))]; Sgt_spec2 = [Sgt_spec2 abs((Sgt2))];
    Sgt_spec3 = [Sgt_spec3 abs((Sgt3))];
    if rem(ii+50, 50) == 0
        Spec1 = [Spec1; Sgt_spec];
        Spec2 = [Spec2; Sgt_spec2];
        Spec3 = [Spec3; Sgt_spec3];
        Sgt_spec = [];
        Sgt_spec2 = [];
        Sgt_spec3 = [];
    end
end

%%


Sgtot = [Spec1; Spec2; Spec3];
clear Sgt_spec; clear Sgt_spec2; clear Sgt_spec3; 
clear Sgt1; clear Sgt2; clear Sgt3;
[u,s,v] = svd(Sgtot',0);
%%
nberr = []; ldaerr = [];
nbstd = []; ldastd = [];
for kkk = 1:1
    
    clc;
    rowbeg = 1;
    rowend = 4;
    trainp = 10;
    testp = 15 - trainp;
    allanswers = [ones(15,1); 2*ones(15,1); 3*ones(15,1)];
    Efinal = 0; Enbfinal = 0;
    Pf = 0; Col = 0; Deo = 0;
    Errlda = []; Errnb = [];
    for jj = 1:100
        floydtrain = []; coltrain = []; deotrain = [];
        ftest = []; coltest = []; deotest = [];
        fsampl = randsample(15,trainp); ft = setdiff(1:15, fsampl)';
        colsampl = randsample(15,trainp) + 15; colt = setdiff(16:30, colsampl)';
        deosampl = randsample(15,trainp) + 30; deot = setdiff(31:45, deosampl)';
        
        for kk = 1:trainp
            ff = v(fsampl(kk), rowbeg:rowend);
            cc = v(colsampl(kk), rowbeg:rowend);
            dd = v(deosampl(kk), rowbeg:rowend);
            floydtrain = [floydtrain ; ff];
            coltrain = [coltrain ; cc];
            deotrain = [deotrain ; dd];
        end
        for kk = 1:testp
            ff = v(ft(kk), rowbeg:rowend);
            cc = v(colt(kk), rowbeg:rowend);
            dd = v(deot(kk), rowbeg:rowend);
            ftest = [ftest ; ff];
            coltest = [coltest ; cc];
            deotest = [deotest ; dd];
        end
        
        
        train = [floydtrain; coltrain; deotrain];
        tests = [ftest; coltest; deotest];
        correct = [allanswers(11:15); allanswers(26:30); allanswers(41:45)];
        answer = [ones(trainp,1); 2*ones(trainp,1); 3*ones(trainp, 1)];
        
        
        [ind err] = classify(tests, train, answer);
        
        Err = 0;
        for ii = 1:length(ind)
            if ind(ii) == correct(ii)
                Err = Err + 1;
                if ii <=5
                    Pf = Pf +1;
                elseif ii <=10
                    Col = Col + 1;
                else
                    Deo = Deo + 1;
                end
            end
            
        end
        Err = Err/15*100;
        
        Errlda = [Errlda Err];
        
        
        figure(2)
        nb = fitNaiveBayes(train, answer);
        prednb = nb.predict(tests);
        
        figure(2)
        bar(ind);
        set(gca, 'FontSize',16);
        xlabel('Test Songs'); ylabel('Classification');
        
        Err = 0;
        for ii = 1:length(prednb)
            if prednb(ii) == correct(ii)
                Err = Err + 1;
            end
            
        end
        Err = Err/15*100;
        Errnb = [Errnb Err];
        
    end
    
    mean(Errnb)
    std(Errnb)
    
    mean(Errlda)
    std(Errlda)
    

    Pf = Pf/(testp*jj)*100;
    Col = Col/(testp*jj)*100;
    Deo = Deo/(testp*jj)*100;
    Banderr = [Pf Col Deo];
    Bands = ['Floyd', 'COl', 'ddd'];
    
    
    nberr = [nberr mean(Errnb)];
    nbstd = [nbstd std(Errnb)];
    ldaerr = [ldaerr mean(Errlda)];
    ldastd = [ldastd std(Errlda)];
    
    
end


