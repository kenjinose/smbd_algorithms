function [ref,Jcost,Ff,Fb,Fc] = c_pef(data,mu,d,K,Niter,Ff,Fb)

%%% C-PEF algorithm
%%% Inputs
%%% data        - matrix containing data N x J (N-Number of time samples, J-Number of traces)
%%% mu          - step size
%%% d           - delay of the predictors
%%% K           - number of coefficients for each predictor
%%% Niter       - number of iterations
%%% Ff, Fb      - vectors contating the initial values for the coefficents
%%% of the predictors
%%% author: Kenji Nose Filho
%%% email:  kenji.nose@ufabc.edu.br
%%% last modification: 19/01/2022

zTaps   = d-1;
nzTaps  = K - zTaps;
[~,Ne]  = size(data);

I      = eye(K+1);
Wp     = -Ff(2:end);
Wr     = -Fb(1:end-1);
Fc     = conv(Ff,Fb);

for ii=1:Niter
    for kk=1:Ne
        x         = [zeros(2*K+1,1); data(:,kk); zeros(2*K+1,1)];
        [Nsamples,~] = size(x);
        y1        = filter(Ff,1,x);
        y2        = filter(Fb,1,x);
        e         = filter(Fc,1,x);
        Jcost(ii) = mean(abs(e));
        for jj=1:K+1
            Y1(:,jj) = filter(I(:,jj),1,y1);
            Y2(:,jj) = filter(I(:,jj),1,y2);
        end
        gradp(:,kk)  = -(1/Nsamples)*Y2(:,2:end)'*sign(e);
        gradr(:,kk)  = -(1/Nsamples)*Y1(:,1:end-1)'*sign(e);
    end
    gradpp = mean(gradp,2);
    gradrr = mean(gradr,2);
    Wp     = Wp - mu*gradpp;
    Wr     = Wr - mu*gradrr;
    Ff     = [1; zeros(zTaps,1); -Wp(zTaps+1:end)];
    Fb     = [-Wr(1:nzTaps); zeros(zTaps,1); 1];
    Fc     = conv(Ff,Fb);
end

for kk =1:Ne
    ref(:,kk)     = conv2(data(:,kk),Fc,'same');
end

end
