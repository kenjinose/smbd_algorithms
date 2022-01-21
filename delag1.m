function [z,lag] = delag1(ref,a,nlag)

%%% Remove o atrasao e corrige a polaridade inserida pelo filtro de
%%% desconvolução através da correlação com o traço ou a refletividade.

%%% Entradas:
%%% ref:  Entrada de referencia (na x ne)
%%% a:    Entrada a ser corrigida (na x ne)
%%% nlag: 2*nlag+1 corresponde ao numero de lags a serem testados

%%% Saída:
%%% z:    Saída corrigida (Na,Ne)

%%% Autor:              Kenji Nose Filho
%%% Email:              kenjinose@yahoo.com.br
%%% Ultima modificação: 27/05/2013

[na,ne] = size(a);

ganho   = ones(1,ne);
A       = zeros(na,ne,2*nlag+1);

for jj=1:ne
    if norm(a(:,jj))>0
        ganho(1,jj) = norm(a(:,jj),2);
        a(:,jj)     = a(:,jj)./ganho(1,jj);
    end
    if norm(ref(:,jj))>0
        ref(:,jj) = ref(:,jj)./norm(ref(:,jj),2);
    end
end

kk=1;
for ii = nlag:-1:1
    A(1:end-ii,:,kk) = a(ii+1:end,:);
    kk               = kk+1;
end

A(:,:,kk) = a;
kk        = kk+1;

for ii=1:1:nlag
    A(ii+1:end,:,kk) = a(1:end-ii,:);
    kk               = kk+1;
end

for ii = 1:1:2*nlag+1;
    C = A(:,:,ii)'*ref;
    for jj=1:ne
        D(jj) = C(jj,jj);
    end
    E(ii) = mean(D);
end

lag   = find(abs(E)==max(abs(E)));
if numel(lag)==0
    sinal = 0;
    z     = a;
else
    sinal = sign(E(lag(1)));
    z     = sinal*A(:,:,lag(1));
end

for jj=1:ne
    if norm(z(:,jj))>0
        z(:,jj)     = z(:,jj).*ganho(1,jj);
    end
end

end

