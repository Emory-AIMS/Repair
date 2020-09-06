
clc;
clear all;
 
%Import neccessary packages
addpath(genpath('./N-Way'));
addpath(genpath('./tensor_toolbox'));
addpath(genpath('./mmttkrp_parafac2'));
 
load('mimicdata.mat'); % load the data set
%load('CMS_500K.mat','X');
% load('parafac2_problem.mat','X');
 
X = mimicdata;
%X = X(1:5000,:);
K = max(size(X)); %number of subjects (K)
R = 4; %number of factors or components
Basis = 7;
PARFOR_FLAG=1;
%normX=claculate_norm(X,K,PARFOR_FLAG,Omega); %Calculate the norm of the input X
conv_tol=1e-3; %converegance tolerance
Smoothness=1; %if smoothness is 1, the U_k will be smooth otherwise not.
 
%Defining the constraints here:
%you can use the folowing constraints in COPA:
%1-non-negativity on H, S_k, and V
%2-Smoothness on U_k
%3-sparsity constraint (l1 and l0) on H, S_k, V
 
%the first element impose a constraint on H
%the second element impose a constraint on V
%the third element impose a constraint on W (S_k)
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
    %X_obs{k} = X_obs{k};
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
 
% save('X_obs_50_0202.mat','X_obs');
% save('Omega_50_0202.mat','Omega');


 
normXgt =   claculate_norm(X,K,PARFOR_FLAG,Omega); %Calculate the norm of the input Xgt
normXobs = claculate_norm_observe(X_obs,K,PARFOR_FLAG,Omega);
X_obs_even = cell(K,1);
Omega_even = cell(k,1);
for h = 1:K
    Omega_full = zeros(size(X_obs{h}));
    Omega_full(Omega{h,1}) = 1;
    Xbasis = X_obs{h};
    prompt = ones(41-min(size(Xbasis)),198);
    Xafter = [Xbasis;prompt];
    %Omegabasis = Omega{h};
    Omegaafter = [Omega_full;prompt];
    X_obs_even{h} = Xafter;
    Omega_even{h} = Omegaafter;    
end

Xobs3D = permute(cat(3,X_obs_even{:}),[3 1 2]);
Omega3D = permute(cat(3,Omega_even{:}),[3 1 2]);

Xobstensor = tensor(Xobs3D);
Omegatensor = tensor(Omega3D);

M_init = create_guess('Data', Xobstensor, 'Num_Factors', R, ...
    'Factor_Generator', 'nvecs');
[P, P0, output] = cp_wopt(Xobstensor,Omegatensor,R,'init', M_init)

result = full(P);
    
fit = 0;
for k = 1:K
    resultslice = result(1,:,:);
    count = min(size(X{k}));
    M = resultslice(1:count,:);
    fit = fit + sum(sum( (full(X{k}) - double(M)).^2));
end
fit=1-(fit/normXgt);
        



