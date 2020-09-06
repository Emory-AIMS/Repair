function [ H, HforX, UL, U, GG, itr ] = fastADMM( Y, H, HforX, U, UL, d, GG, Constraints,PARFOR_FLAG)
% ADMM iterates to solve one mode of tensor Y. 
%based on formula 8 in the paper.
%% tune
betaL = [1e-3 1e-4 1e-4];

K=size(Y,1);
itr=0;
[ ~, k ] = size(H{d});
G = ones(k,k); prod = [ 1:d-1, d+1:length(GG) ];
for dd = prod
    G = G .* GG{dd}; 
end


rho =min( 1e-3,trace(G)/k);%% inverse of step size for internal fast-admm

%rho = rho * 10;
rhoL = rho; %%Tune it
L = chol( G + (rho)*eye(k), 'lower' );

%F = mttkrp( Y, H, d );
F=mttkrp_for_parafac2(Y, K, H, d, PARFOR_FLAG);
tol = 1e-2;
Hd = H{d}; Ud = U{d}; Hl = HforX{d}; Ul=UL{d};
for itr = 1:1 %you can change the number of inner iterations
    H0 = Hd;

    Ht   = L'\ ( L\ ( F + rho*(Hd+Ud) +rhoL*(Hl+Ul))' );
    Hd = proxr( Ht'-Ud, Constraints, d, rho);
    Hl = proxnn(Ht'-Ul, betaL, d, rhoL);
    Ud = Ud + Hd - Ht';
    Ul = Ul + Hl - Ht';
end

U{d} = Ud;
UL{d} = Ul;

    
H{d} = Hd; %% original
GG{d} = Hd'*Hd; %% original
% H{d} = Ht';
HforX{d} = Hl;  
% GG{d} = Ht*Ht';
end


function H = proxr( Hb, Constraints, d, rho )
    switch Constraints{d}
        case 'nonnegative'
            H = max( 0, Hb );
        case 'l1'
            l1_regul=0.0000085;
            H = sign( Hb ) .* max( 0, abs(Hb) - (l1_regul/rho) ); 
        case 'l0'
            l0_regul=8;
            Hb(Hb <= l0_regul) = 0; 
            H=Hb;
    end
end

function H = proxnn(Hb, betaL,d,rhoL)

    threshS = betaL(d) / rhoL;
    [U,S,V] = svd(Hb,'econ');
    sv = diag(S);
    nnzsv = length(find(sv>threshS));
    if  nnzsv>= 1
        sv = sv(1:nnzsv) - threshS;
        H = U(:,1:nnzsv)*diag(sv)*V(:,1:nnzsv)';
    else
        H = zeros(size(Hb));
    end

end
