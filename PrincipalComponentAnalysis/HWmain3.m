close all; clear all; clc;

load('cam1_3.mat');
load('cam2_3.mat');
load('cam3_3.mat');
vidFrames1_1 = vidFrames1_3;
vidFrames2_1 = vidFrames2_3;
vidFrames3_1 = vidFrames3_3;

[m1,n1] = size(vidFrames1_1(:,:,1,1));
[m2,n2] = size(vidFrames2_1(:,:,1,1));
[m3,n3] = size(vidFrames3_1(:,:,1,1));

numFrames1=size(vidFrames1_1,4);
numFrames2=size(vidFrames2_1,4);
numFrames3=size(vidFrames3_1,4);
maxFrames = max([numFrames1 numFrames2 numFrames3]);

for k = 1:maxFrames
    if k <= numFrames1
        mov1(k).cdata = vidFrames1_1(:,:,:,k);
        mov1(k).colormap = [];
        
    end
    if k <= numFrames2
        mov2(k).cdata = vidFrames2_1(:,:,:,k);
        mov2(k).colormap = [];
    end
    if k <= numFrames3
         mov3(k).cdata = vidFrames3_1(:,:,:,k);
         mov3(k).colormap = [];
    end  
end

X1 = []; X2 = []; X3 = [];
Y1 = []; Y2 = []; Y3 = [];
for jj = 1:maxFrames
    
    if jj <= numFrames1
        X=frame2im(mov1(jj));
        a = rgb2gray(X);
        Vid1org(:,:,jj) = a;
        a(:,1:290) = 0; a(:,350:end) = 0;
        a(1:200,:) = 0;
        Vid1(:,:,jj) =  a;
        [Max Ind] = max(a(:));
        [x1 y1] = ind2sub(size(a), Ind);
        X1 = [X1 x1]; Y1 = [Y1 y1];
    end
    
    if jj <= numFrames2
        X=frame2im(mov2(jj));
        a = rgb2gray(X);
        Vid2org(:,:,jj) = a;
        a(:,1:230) = 0; a(:,350:end) = 0;
        Vid2(:,:,jj) =  a;
        [Max Ind] = max(a(:));
        [x2 y2] = ind2sub(size(a), Ind);
        X2 = [X2 x2]; Y2 = [Y2 y2];
    end 
 
    if jj <= numFrames3
        X=frame2im(mov3(jj));
        a = rgb2gray(X);
        Vid3org(:,:,jj) = a;
        a(1:200,:) = 0; a(300:end,:) = 0;
        a(:, 1:260) = 0;
        Vid3(:,:,jj) =  a; 
        [Max Ind] = max(a(:));
        [x3 y3] = ind2sub(size(a), Ind);
        X3 = [X3 x3]; Y3 = [Y3 y3];
    end

end

%%
[Min1 I1]=min(X1(1:50)); X1=X1(I1:I1+200); Y1=Y1(I1:I1+200);
[Min2 I2]=min(X2(1:50)); X2=X2(I2:I2+200); Y2=Y2(I2:I2+200);
[Min3 I3]=min(Y3(1:50)); Y3=Y3(I3:I3+200); X3=X3(I3:I3+200);

X1f = fft(X1); X2f = fft(X2); X3f = fft(X3);
Y1f = fft(Y1); Y2f = fft(Y2); Y3f = fft(Y3);
X1f(10:end-10) = 0; X2f(10:end-10) = 0; X3f(10:end-10) = 0;
Y1f(10:end-10) = 0; Y2f(10:end-10) = 0; Y3f(10:end-10) = 0;
Xo1 = abs(ifft(X1f)); Yo1 = abs(ifft(Y1f));
Xo2 = abs(ifft(X2f)); Yo2 = abs(ifft(Y2f));
Xo3 = abs(ifft(X3f)); Yo3 = abs(ifft(Y3f));
%% 
figure (1)
subplot(3,2,1), plot(X1, 'k', 'LineWidth', 1.5); hold on;
plot(Xo1, 'b--', 'LineWidth', 0.25)
axis([0 200 0 480]); ylabel('Position in Y');
subplot(3,2,2), plot(Y1, 'k', 'LineWidth', 1.5); hold on;
plot(Yo1, 'b--', 'LineWidth', 0.25)
axis([0 200 0 640]); ylabel('Position in X');
subplot(3,2,3), plot(X2, 'k', 'LineWidth', 1.5); hold on;
plot(Xo2, 'b--', 'LineWidth', 0.25)
axis([0 200 0 480]); ylabel('Position in Y');
subplot(3,2,4), plot(Y2, 'k', 'LineWidth', 1.5); hold on;
plot(Yo2, 'b--', 'LineWidth', 0.25)
axis([0 200 0 640]); ylabel('Position in X');
subplot(3,2,5), plot(X3, 'k', 'LineWidth', 1.5); hold on;
plot(Xo3, 'b--', 'LineWidth', 0.25)
axis([0 200 0 480]); ylabel('Position in Y'); xlabel('Frame');
subplot(3,2,6), plot(Y3, 'k', 'LineWidth', 1.5); hold on;
plot(Yo3, 'b--', 'LineWidth', 0.25)
axis([0 200 0 640]); ylabel('Position in X'); xlabel('Frame');

Amat = [X1; Y1; X2; Y2; X3; Y3];
Aomat = [Xo1; Yo1; Xo2; Yo2; Xo3; Yo3];

[u s v] = svd(Amat);
[uo so vo] = svd(Aomat);
%%
% figure(2)
% plot(diag(s)*100/sum(abs(diag(s))),'ro');
% axis([1 6 0 100]);

Anew = u(:,1:2)*s(1:2,1:2)*v(:,1:2)';
Aonew = uo(:,1:2)*so(1:2,1:2)*vo(:,1:2)';
%Anew = u*s*v';


Xnew1 = Anew(1,:); Ynew1 = Anew(2,:);
Xnew2 = Anew(3,:); Ynew2 = Anew(4,:);
Xnew3 = Anew(5,:); Ynew3 = Anew(6,:);
Xonew1 = Aonew(1,:); Yonew1 = Aonew(2,:);
Xonew2 = Aonew(3,:); Yonew2 = Aonew(4,:);
Xonew3 = Aonew(5,:); Yonew3 = Aonew(6,:);

figure (2)
subplot(3,2,1), plot(Xnew1, 'k', 'LineWidth', 1.5); hold on;
plot(Xonew1, 'b--', 'LineWidth', 0.25)
axis([0 200 0 480]); ylabel('Position in Y');
subplot(3,2,2), plot(Ynew1, 'k', 'LineWidth', 1.5); hold on;
plot(Yonew1, 'b--', 'LineWidth', 0.25)
axis([0 200 0 640]); ylabel('Position in X');
subplot(3,2,3), plot(Xnew2, 'k', 'LineWidth', 1.5); hold on;
plot(Xonew2, 'b--', 'LineWidth', 0.25)
axis([0 200 0 480]); ylabel('Position in Y');
subplot(3,2,4), plot(Ynew2, 'k', 'LineWidth', 1.5); hold on;
plot(Yonew2, 'b--', 'LineWidth', 0.25)
axis([0 200 0 640]); ylabel('Position in X');
subplot(3,2,5), plot(Xnew3, 'k', 'LineWidth', 1.5); hold on;
plot(Xonew3, 'b--', 'LineWidth', 0.25)
axis([0 200 0 480]); xlabel('Frame'); ylabel('Position in Y');
subplot(3,2,6), plot(Ynew3, 'k', 'LineWidth', 1.5); hold on;
plot(Yonew3, 'b--', 'LineWidth', 0.25)
axis([0 200 0 640]); xlabel('Frame'); ylabel('Position in X');

figure(3)
plot(diag(s)*100/sum(diag(s)), 'ko--', 'MarkerSize', 8'); hold on;
axis([1 6 0 100]); xlabel('Modes'); ylabel('% of Energy');




% figure(4)
% for jj = 1:200
%     A = Vid1org(:, :, jj+I1);
%     A(floor(X1(jj))-5:floor(X1(jj))+5, floor(Y1(jj))-5:floor(Y1(jj))+5) = 0;
%     imshow(A); drawnow;
% 
% end



