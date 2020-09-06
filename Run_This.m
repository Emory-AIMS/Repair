clc;
clear all;
 
%Import neccessary packages
addpath(genpath('./N-Way'));
addpath(genpath('./tensor_toolbox'));
addpath(genpath('./mmttkrp_parafac2'));
 
load('mimicdata.mat'); % load the data set
%load('CMS_500K.mat','X');

 
X = mimicdata;
K = max(size(X)); %number of subjects (K)
R = 4; %number of factors or components
Basis = 7;
PARFOR_FLAG=1;
%normX=claculate_norm(X,K,PARFOR_FLAG,Omega); %Calculate the norm of the input X
conv_tol=1e-3; %converegance tolerance
Smoothness=1; %if smoothness is 1, the U_k will be smooth otherwise not.
 
Constraints={'nonnegative', 'l1','nonnegative'};
seed=2;
 
%% index of observed entries
Omega = cell(K,1);
X_obs = cell(K,1);
 
 
%%add error and missing entries
missing_p = 0.3;
error_p = 0.3;
for k = 1:K
    [rr,ll] = size(X{k});
    X_obs{k} = zeros(rr,ll);
    X_tmp = full(X{k});
    % [~,error_pos] = datasample(1:rr*ll,floor(nnz(X{k})*error_p));
    error_pos = randi([1 rr*ll],1,floor(nnz(X{k})*error_p));
    X_obs{k}(error_pos) = 4;
%     error_val = randi(4);
%     if error_val==X_obs{k}(error_pos) 
%         switch error_val
%             case 1
%                 X_obs{k}(error_pos) = 2;
%             otherwise
%                 X_obs{k}(error_pos) = 1;
%         end 
%     else
%         X_obs{k}(error_pos) = error_val;
%     end
    true_pos = setdiff(1:rr*ll,error_pos);
    % observed_pos = datasample(true_pos,floor(length(true_pos)*(1-missing_p)));
    observed_pos_ind = randi([1 numel(true_pos)],1,floor(length(true_pos)*(1-missing_p)));
    observed_pos = true_pos(observed_pos_ind);
    X_obs{k}(observed_pos) = X{k}(observed_pos);
    Omega{k} = union(observed_pos,error_pos);
    X_obs{k} = sparse(X_obs{k});
% 
%     [rr,ll] = size(X{k});
%   X_obs{k} = zeros(rr,ll);
%   X_tmp = full(X{k});
%   [~,error_pos] = datasample(1:rr*ll,floor(nnz(X{k})*error_p));
%   X_obs{k}(error_pos) = randi([1,3],1,1);
%   true_pos = find(X_tmp);
%     Omega_gt{k} = find(X_tmp);
%     observed_pos = datasample(true_pos,floor(length(true_pos)*(1-missing_p)));
%   X_obs{k}(observed_pos) = X{k}(observed_pos);
%   Omega{k} = find(X_obs{k});
%   X_obs{k} = sparse(X_obs{k});
end
 
 
 
normXgt =   claculate_norm(X,K,PARFOR_FLAG,Omega); %Calculate the norm of the input Xgt
normXobs = claculate_norm_observe(X_obs,K,PARFOR_FLAG,Omega);
 
if Smoothness==1
     
    GAP=0;
      %[fit1,U,V,W]=Smooth_COPA(X,X_obs,Omega,R,conv_tol,seed,PARFOR_FLAG,normXgt,normXobs,Constraints,GAP,Basis);
 

        [fit2,U,V,W]=Robust_Smooth_COPA(X,X_obs,Omega,R,conv_tol,seed,PARFOR_FLAG,normXgt,normXobs,Constraints,GAP,Basis);
 
else
    tStart=tic;
    [fit,FIT_TIME]=COPA(X_obs,R,conv_tol,seed,PARFOR_FLAG,normX,Constraints);
    tEnd = toc(tStart);
end
 
 
 




