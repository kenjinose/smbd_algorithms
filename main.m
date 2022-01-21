%%% 

clc;
clear;
close all;

load synthetic2

%catches the default stream used by the random generator
defaultStream = RandStream.getGlobalStream;

%establishes the random seed
stream_exp = RandStream('mt19937ar','seed',1); defaultStream.State = stream_exp.State;

r   = Data(1:301,1:20);
dt  = 2e-3;
t   = 0:dt:(length(r)-1)*dt;
cdp = 1:size(r,2);

% Ricker wavelet
fc = 40;
tw = -0.03:dt:0.03;
h  = (1-2*pi^2*fc^2*tw.^2).*exp(-pi^2*fc^2*tw.^2);
hh = hilbert(h);
th = -50*pi/180;
h  = cos(th)*real(hh)+sin(th)*imag(hh);
h  = h';

% Numerical mixed phase wavelet
% h  = [0.36257276552;0.26004577809;-0.29809191996;-0.57265987622;0.06845741412;0.90000000000;0.09833009553;-0.30119915255;-0.10461525146;0.16908296900;0.08821194769;-0.08750787106;-0.06779517513;0.02762973173;0.03568423704;-0.00398180992;-0.01294844851;0.00040330489;0.00499510078;0.00074562673;-0.00148674943;-0.00066787272;0.00013898780;0.00018550224;0.00004195727;-0.00002480371;-0.00002472178;0.00000062429;0.00001311489;0.00000540320];

h  = h/norm(h,'inf');
x  = filter(h,1,r);

%% Aditive White Gaussian Noise
n              = randn(size(r));
x              = x./std(x(:));
n              = n./std(n(:));
SNR_dB         = 1000
SNR            = 10^(SNR_dB/10)
alpha          = 1./sqrt(SNR);
n              = alpha*n;
x_awgn         = x + n;
[x_awgn,~]     = delag1(r,x_awgn,length(h));

[Nsamples,Ntraces] = size(x_awgn);

%% C-PEF

Niter   = 1000;
mu      = 0.1;
d       = 1;
K       = 25;
Ff      = [1; zeros(K,1);];
Fb      = [zeros(K,1); 1];

tic
[ref_pef,Jcost_pef,~,~,~] = c_pef(x_awgn,mu,d,K,Niter,Ff,Fb);
toc
[ref_pef,~] = delag1(r,ref_pef,length(h));

sfont = 14;

figure(1),
plot(Jcost_pef), grid on
title('Cost Function x Iteration - C-PEF','FontSize',sfont,'Fontname','Arial')
xlabel('Iteration','FontSize',sfont,'Fontname','Arial'), ylabel('Cost Function','FontSize', sfont,'Fontname','Arial')

% Calculate wavelet from inverse filter

%% AM-SMBD

Niter   = 200;   % Number of iterations
mu      = 0.05;  % Step size
lambda  = 0.05;  % Regularization parameter (reflectivity)
sigma   = 0.10;  % Regularization parameter (wavelet)
L       = length(h); % Length of the wavelet

tic
[ref_am,wavelet_am,wavelet0_am,JcostMSE_am,JcostL1_am] =  am_smbd(x_awgn,mu,lambda,sigma,Niter,L);
toc
[ref_am,~] = delag1(r,ref_am,L);

figure(2),
title('Cost Function x Iteration - AM-SMBD','FontSize',sfont,'Fontname','Arial')
subplot(131), plot(JcostMSE_am), grid on
xlabel('Iteration','FontSize',sfont,'Fontname','Arial'), ylabel('J_{MSE}','FontSize', sfont,'Fontname','Arial')
subplot(132), plot(JcostL1_am), grid on
xlabel('Iteration','FontSize',sfont,'Fontname','Arial'), ylabel('J_{l1}','FontSize', sfont,'Fontname','Arial')
subplot(133), plot(JcostMSE_am+lambda*JcostL1_am), grid on
xlabel('Iteration','FontSize',sfont,'Fontname','Arial'), ylabel('J_{MSE}+\lambdaJ_{l1}','FontSize', sfont,'Fontname','Arial')

%% FISTA-SMBD

Niter  = 1000;      % Number of iterations
mu     = 0.05;      % Step size
lambda = 0.001;     % Regularization parameter
L      = length(h); % Length of the wavelet

tic
[ref_smbd,JcostMSE_smbd,JcostL1_smbd] = smbd_fista(x_awgn,mu,lambda,Niter,L);
toc
[ref_smbd,~] = delag1(r,ref_smbd,L);

figure(3),
title('Cost Function x Iteration - FISTA-SMBD','FontSize',sfont,'Fontname','Arial')
subplot(131), plot(JcostMSE_smbd), grid on
xlabel('Iteration','FontSize',sfont,'Fontname','Arial'), ylabel('J_{MSE}','FontSize', sfont,'Fontname','Arial')
subplot(132), plot(JcostL1_smbd), grid on
xlabel('Iteration','FontSize',sfont,'Fontname','Arial'), ylabel('J_{l1}','FontSize', sfont,'Fontname','Arial')
subplot(133), plot(JcostMSE_smbd+lambda*JcostL1_smbd), grid on
xlabel('Iteration','FontSize',sfont,'Fontname','Arial'), ylabel('J_{MSE}+\lambdaJ_{l1}','FontSize', sfont,'Fontname','Arial')

wavelet0_am = delag1(h,wavelet0_am,L);
wavelet_am  = delag1(h,wavelet_am,L);

figure(4),
subplot(256),
colormap(gray),
imagesc(cdp,t,r),
title('True Reflectivity','FontSize',sfont,'Fontname','Arial')
xlabel('CMP','FontSize',sfont,'Fontname','Arial'), ylabel('Time(s)','FontSize', sfont,'Fontname','Arial')
axis(set(gca,'FontSize',sfont,'Fontname','Arial'));

subplot(651),
plot((0:length(h)-1)*dt , h./max(abs(h)),'k','LineWidth',2), grid on
title('True wavelet','FontSize',sfont,'Fontname','Arial')
xlabel('Time (s)','FontSize',sfont,'Fontname','Arial'), ylabel('Amplitude','FontSize', sfont,'Fontname','Arial')
axis(set(gca,'FontSize',sfont,'Fontname','Arial'));

subplot(257),
colormap(gray),
imagesc(cdp,t,x_awgn),
title('Synthetic Trace','FontSize',sfont,'Fontname','Arial')
xlabel('CMP','FontSize',sfont,'Fontname','Arial'), ylabel('Time(s)','FontSize', sfont,'Fontname','Arial')
axis(set(gca,'FontSize',sfont,'Fontname','Arial'));

subplot(652),
plot((0:length(h)-1)*dt , wavelet0_am./max(abs(wavelet0_am)),'k','LineWidth',2), grid on
title('Zero phase wavelet','FontSize',sfont,'Fontname','Arial')
xlabel('Time (s)','FontSize',sfont,'Fontname','Arial'), ylabel('Amplitude','FontSize', sfont,'Fontname','Arial')
axis(set(gca,'FontSize',sfont,'Fontname','Arial'));


subplot(258),
colormap(gray),
imagesc(cdp,t,ref_pef),
title('C-PEF','FontSize',sfont,'Fontname','Arial')
xlabel('CMP','FontSize',sfont,'Fontname','Arial'), ylabel('Time(s)','FontSize', sfont,'Fontname','Arial')
axis(set(gca,'FontSize',sfont,'Fontname','Arial'));

subplot(259),
colormap(gray),
imagesc(cdp,t,ref_am),
title('AM-SMBD','FontSize',sfont,'Fontname','Arial')
xlabel('CMP','FontSize',sfont,'Fontname','Arial'), ylabel('Time(s)','FontSize', sfont,'Fontname','Arial')
axis(set(gca,'FontSize',sfont,'Fontname','Arial'));

subplot(654),
plot((0:length(h)-1)*dt , wavelet_am./max(abs(wavelet_am)),'k','LineWidth',2), grid on
title('Estimated wavelet','FontSize',sfont,'Fontname','Arial')
xlabel('Time (s)','FontSize',sfont,'Fontname','Arial'), ylabel('Amplitude','FontSize', sfont,'Fontname','Arial')
axis(set(gca,'FontSize',sfont,'Fontname','Arial'));

subplot(2,5,10)
colormap(gray)
imagesc(cdp,t,ref_smbd),
title('FISTA-SMBD','FontSize',sfont,'Fontname','Arial')
xlabel('CMP','FontSize',sfont,'Fontname','Arial'), ylabel('Time(s)','FontSize', sfont,'Fontname','Arial')
axis(set(gca,'FontSize',sfont,'Fontname','Arial'));

position1 = [10 51 800 400];
position2 = [10 51 1600 800];

set(1,'Position',position1)
set(2,'Position',position1)
set(3,'Position',position1)
set(4,'Position',position2)