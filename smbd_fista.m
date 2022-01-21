function [ref,JcostMSE,JcostL1] = smbd_fista(data,mu,lambda,Niter,L)

%%% FISTA-SMBD algorithm
%%% Inputs
%%% data        - matrix containing data N x J (N-Number of time samples, J-Number of traces)
%%% mu          - step size
%%% lambda      - regularization parameter (l1-norm, reflectivity)
%%% Niter       - number of iterations
%%% L           - length of the wavelet
%%% author: Kenji Nose Filho
%%% email:  kenji.nose@ufabc.edu.br
%%% last modification: 19/01/2022

gain = norm(data,'fro');
data = data/gain;

[Ns, Nx] = size(data);
N        = Ns-L+1;
ref      = data(1:N,:);
Df       = zeros(2*Ns+1,Nx);
for ii=1:Nx
    Df(:,ii)=fft(data(:,ii),2*Ns+1);
end

[idx3] = select_traces(ref);

refy = ref;
ref0 = ref;
t    = 1;
ll   = 1;
kk   = 1;
JcostMSE = zeros(Niter,1);
JcostL1  = zeros(Niter,1);

for kk = 1:Niter
    t0    = t;
    grad  = multiplyR(refy,Df,idx3,Nx,N,Ns); %%% A'Ax 
    ref   = wthresh(refy-2*mu*grad,'s',lambda*mu);
    ref   = ref/(norm(ref(:))+eps);
    t     = (1 + sqrt(1+4*t0^2))/2;
    alpha = (t0-1)/t;
    refy  = ref + alpha*(ref - ref0);    
    ref0  = ref;
    mu    = 0.99*mu;

    JcostMSE(kk) = sum(sum(ref.*multiplyR(ref,Df,idx3,Nx,N,Ns)));
    JcostL1(kk) = sum(abs(ref(:)));
  
end

ref = gain*[ref; zeros(L-1,Nx)];

end

function [idx3] = select_traces(X)

nt = size(X,2);

ll = 1;
for ii = 1:nt-2
    idx2(ll,1)   = ii;
    idx2(ll,2)   = -ii-1;
    idx2(ll+1,1) = ii;
    idx2(ll+1,2) = -ii-2;
    ll = ll+2;
end
ii = nt-1;
idx2(ll,1) = ii;
idx2(ll,2) = -ii-1;

ll = 1;
for kk = 1:size(idx2,1)
    idx3(ll,:)   = [kk  idx2(kk,1)  idx2(kk,1) abs(idx2(kk,2)) abs(idx2(kk,2))];
    idx3(ll+1,:) = [kk  idx2(kk,2)  idx2(kk,2) abs(idx2(kk,1)) abs(idx2(kk,1))];
    idx3(ll+2,:) = [kk  idx2(kk,1)  idx2(kk,2) abs(idx2(kk,2)) abs(idx2(kk,1))];
    idx3(ll+3,:) = [kk  idx2(kk,2)  idx2(kk,1) abs(idx2(kk,1)) abs(idx2(kk,2))];
    ll           = ll+4;
end
end

function [Xout] = multiplyR(Xk,Df,idx3,Nx,L,Ns)
Xf   = zeros(2*Ns+1,Nx);
Xout = zeros(L,Nx);
for jj=1:Nx
    Xf(:,jj) = fft(Xk(:,jj),2*Ns+1);
end

for ll=1:size(idx3,1)
    jj         = idx3(ll,2);
    ii         = idx3(ll,3);
    mm         = idx3(ll,4);
    nn         = idx3(ll,5);
    aux        = sign(ii)*sign(jj)*ifft(conj(Df(:,abs(ii))).*Df(:,abs(jj)).*Xf(:,mm));
    Xout(:,nn) = Xout(:,nn)+aux(1:L,1);
end

end
