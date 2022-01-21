function [ref,wavelet,wavelet0,JcostMSE,JcostL1] = am_smbd(data,mu,lambda,sigma,Niter,L)

%%% AM-SMBD algorithm
%%% Inputs
%%% data        - matrix containing data N x J (N-Number of time samples, J-Number of traces)
%%% mu          - step size
%%% lambda      - regularization parameter (l1-norm, reflectivity)
%%% sigma       - regularization parameter (l2-norm, wavelet)
%%% Niter       - number of iterations
%%% L           - length of the wavelet
%%% author: Kenji Nose Filho
%%% email:  kenji.nose@ufabc.edu.br
%%% last modification: 19/01/2022

[Nsamples,Ntraces] = size(data);

data = [data; zeros(L,Ntraces)];

[Nsamples,Ntraces] = size(data);

% Initialize reflectivity series
ref = data;

% Temporal window
wind            = ones(Nsamples,1);
Lini            = floor(L/2);
Lfim            = Nsamples-floor(L/2);
wind(Lini:Lfim) = 0;

% Initialize wavelet
ampX      = abs(fft(data));
ampX      = mean(ampX,2);

wavelet   = real(ifft(ampX)).*wind;
wavelet0  = fftshift(wavelet);

W  = fft(wavelet);
R  = fft(ref);
X  = fft(data);
A0 = abs(X);

for kk = 1:Ntraces
    JcostMSE(1,kk) = (1/Nsamples)*((X(:,kk)-R(:,kk).*W)'*(X(:,kk)-R(:,kk).*W));
    JcostL1(1,kk)  = sum(abs(data(:,kk)));
end

Ry = R;
R0 = R;
t  = ones(Ntraces,1);
t0 = ones(Ntraces,1);

for ii = 1:Niter
    
    for kk = 1:Ntraces
        
        t0(kk)  = t(kk);
        gradRy  = (1/Nsamples)*(2*(Ry(:,kk).*W.*conj(W)-X(:,kk).*conj(W)));
        r       = wthresh(real(ifft(Ry(:,kk) - 2*mu*gradRy)),'s',lambda*mu);
        R(:,kk) = fft(r);
        t(kk)   = (1 + sqrt(1+4*t0(kk)^2))/2;
        alpha   = (t0(kk)-1)/t(kk);
        Ry(:,kk) = R(:,kk)+alpha*(R(:,kk)-R0(:,kk));
        R0(:,kk) = R(:,kk);
        
        % Compute gradient of W
        
        vecW(:,kk) = conj(R(:,kk)).*X(:,kk)./(conj(R(:,kk)).*R(:,kk)+sigma);
        
        JcostMSE(ii+1,kk) = (1/Nsamples)*((X(:,kk)-R(:,kk).*W)'*(X(:,kk)-R(:,kk).*W));
        rT(:,kk)          = real(ifft(R(:,kk)));
        JcostL1(ii+1,kk)  = sum(abs(rT(:,kk)));
        
    end
    
    % Update W
    
    W         = mean(vecW,2);
    wavelet   = real(ifft(W)).*wind;
    W         = fft(wavelet);
    
end

ref     = real(ifft(R));
wavelet = fftshift(real(ifft(W)));

wavelet0 = wavelet0(round(length(wavelet0)/2+1)-Lini:round(length(wavelet0)/2)-Lini+L);
wavelet  = wavelet(round(length(wavelet)/2+1)-Lini:round(length(wavelet)/2)-Lini+L);

ref = ref(1:end-L,:);

end

