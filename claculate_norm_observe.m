function [ normX ] = claculate_norm_observe(X,K,PARFOR_FLAG,Omega)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    normX=0;
    if (PARFOR_FLAG)
        parfor k = 1:K
         Xk = X{k};
         indexK = Omega{k};
         normX = normX + sum(sum( (Xk(indexK)).^2));
        end
    else
        for k = 1:K
         Xk = X{k};
         indexK = Omega{k};
         normX = normX + sum(sum( (Xk(indexK)).^2));
        end
    end


end

