function [robust_fit_gt,U,V,W] = Robust_Smooth_COPA(X_gt, X, Omega, R,conv_tol,seed,PARFOR_FLAG,normXgt,normXobs,Constraints,GAP,Basis )
%Implementation of smooth parafac 2 where smoothness impose on mode U_k
% If GAP=1 then the smoothness considers the gap between two time stamps

%% tune 
rhoO =3; %% inverse of step size for recovry XClean
alphaE = 1e-3;
rhoOMax = 50;


%tStart=tic;
%FIT_TIME=[];

J=size(X{1}, 2); %  number of features
K = max(size(X));% number of subjects

E = cell(K,1);
UX = cell(K,1);
Xclean = cell(K,1);
OmegaBot = cell(K,1);

if(PARFOR_FLAG)
    parfor k=1:K
        E{k} = sparse(zeros(size(X{k})));
        Xclean{k} = sparse(X{k});
        [rr,cc] = size(X{k});
        % Xclean{k} = zeros(size(X{k}));
        UX{k} = sparse(zeros(size(X{k})));
        OmegaBot{k} = setdiff(1:(rr*cc),Omega{k});
    end
else
    for k=1:K
        E{k} = sparse(zeros(size(X{k})));
        Xclean{k} = sparse(X{k});
        [rr,cc] = size(X{k});
        % Xclean{k} = zeros(size(X{k}));
        UX{k} = sparse(zeros(size(X{k})));
        OmegaBot{k} = setdiff(1:(rr*cc),Omega{k});
    end
end

if(GAP==1)
    %fid = fopen('gaps_per_subject.csv','rt');
    % the file format should be like this:
    %each line is related to a subject. (4,7,12,19,39,45,......)
    %each number is a time stamp.
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F={};%containts all basis functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ra=cell(K,1);
num_of_basis_fun=Basis;
if(PARFOR_FLAG)
    parfor k=1:K
        %create B spline for each subject
        if(GAP==1) % Incorporating the observational gap for creating the smooth function.
            knots = str2num(fgetl(fid)) %read the next line of the days of a patient
            knots=knots-(knots(1)-1);
            patient_dist=knots/knots(end);
            knots=knots*size(X{k},1);
            
        else
            knots=1:size(X{k},1);
        end
        
        F{k}=MSplineBasis([knots], num_of_basis_fun,3, [knots(1) knots(end)] ); 
        [u,~,~]=svd(F{k},'econ');
        Ra{k}=u*u';

    end
else
    for k=1:K
        if(GAP==1)
            knots = str2num(fgetl(fid)); %read the next line of the days of a patient
            knots=knots-(knots(1)-1);
            knots=knots/knots(end);
            knots=knots*size(X{k},1);
            
        else
            knots=1:size(X{k},1);
        end
        F{k}=MSplineBasis([knots], num_of_basis_fun,3, [knots(1) knots(end)] ); 
        [u,~,~]=svd(F{k},'econ');
        Ra{k}=u*u';

    end
end

rng(seed)
H=rand(R,R);
V=rand(J,R);
W=rand(K,R);

HforX=H;
VforX=V;
WforX=W;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
itr=0;
Rknew=cell(K,1);
Xtilde=cell(K,1);
prev_fit=0; robust_fit_gt=1;

%myfilepath = './RobustCopa.txt';
%fid = fopen( sprintf( 'file%i.txt',j) );
%fid = fopen(myfilepath, 'wt');
%fid = fopen( ['file' num2str(j) '.txt'] );
%fprintf(fid, '%f\n',alphaE); 
while(abs(robust_fit_gt-prev_fit)>conv_tol*prev_fit)
%while(1)
    prev_fit=robust_fit_gt;
    if(PARFOR_FLAG)
        parfor k=1:K

            [u,~,v]=svd(Ra{k}*Xclean{k}*V*(diag(W(k,:)))*H','econ');
            Rknew{k}=u*v';
            Xtilde{k} = sparse(Rknew{k}'*Xclean{k});

        end
    else
        for k=1:K

            [u,~,v]=svd(Ra{k}*Xclean{k}*V*(diag(W(k,:)))*H','econ');
            Rknew{k}=u*v';
            Xtilde{k} = sparse(Rknew{k}'*Xclean{k});

        end
    end


        
    [ Tensor,TensorX]  = Robust_COPA_optimizer( Xtilde, R,  'maxiters', 1 ,'init', {H, V, W}, 'initX', {HforX, VforX, WforX}, 'Constraints',Constraints,'PARFOR_FLAG',PARFOR_FLAG);
    H=Tensor{1};
    V=Tensor{2};
    W=Tensor{3};

    HforX=TensorX{1};
    VforX=TensorX{2};
    WforX=TensorX{3};

    if(PARFOR_FLAG)
        parfor k=1:K
            % Tmp = (Rknew{k}*HforX)*diag(WforX(k,:))*(VforX)';
            Xclean{k} = ((Rknew{k}*HforX)*diag(WforX(k,:))*(VforX)'+rhoO*(UX{k}-E{k}+X{k})) / (1.0+rhoO);
            Xclean{k} = min(max(0,round(Xclean{k})),4); %% optional
            % Ektmp = E{k};
            E{k} = max(X{k}+UX{k}-Xclean{k}-alphaE/rhoO,0);
            % DiffClean = X{k} - Xclean{k};
            E{k}(OmegaBot{k}) = X{k}(OmegaBot{k}) - Xclean{k}(OmegaBot{k});
            % E{k} = Ektmp;
            UX{k} = UX{k} + (Xclean{k}+E{k}-X{k});
        end
        
        
    else
        for k=1:K
        % Tmp = (Rknew{k}*HforX)*diag(WforX(k,:))*(VforX)';
            Xclean{k} = ((Rknew{k}*HforX)*diag(WforX(k,:))*(VforX)'+rhoO*(UX{k}-E{k}+X{k})) / (1.0+rhoO);
            Xclean{k} = max(0,round(Xclean{k})); %% optional
            % Ektmp = E{k};
            E{k} = max(X{k}+UX{k}-Xclean{k}-alphaE/rhoO,0);
            % DiffClean = X{k} - Xclean{k};
            E{k}(OmegaBot{k}) = X{k}(OmegaBot{k}) - Xclean{k}(OmegaBot{k});
            % E{k} = Ektmp;
            UX{k} = UX{k} + (Xclean{k}+E{k}-X{k});
%         for k=1:K
%             Tmp = (Rknew{k}*HforX)*diag(WforX(k,:))*(VforX)';
%             Xclean{k} = (Tmp+rhoO*(UX{k}-E{k}+X{k})) / (1.0+rhoO);
%             %Xclean{k} = max(0,round(Xclean{k})); %% optional
%             Ektmp = E{k};
%             Ektmp = max(X{k}+UX{k}-Xclean{k}-alphaE/rhoO,0);
%             DiffClean = X{k} - Xclean{k};
%             Ektmp(OmegaBot{k}) = DiffClean(OmegaBot{k});
%             E{k} = Ektmp;
%             UX{k} = UX{k} + (Xclean{k}+E{k}-X{k});
        end
    end

     
    itr=itr+1;
    
    if(itr > 10)
       rhoO = 40;
    end
    
    robust_fit_gt=calculate_fit(X_gt,Rknew,H,W,V,normXgt,K,PARFOR_FLAG)
    %fprintf(fid, '%d\n',itr); 
    %fprintf(fid, 'RobustCopa = %f\n ',robust_fit_gt); 
    %robust_fit_ob=calculate_fit_obs(X,Rknew,H,W,V,normXobs,K,PARFOR_FLAG,Omega);
    
    %diffnorm = calculate_diffnorm(X_gt,Xclean,Omega,PARFOR_FLAG,normXgt,K)
    %tEnd = toc(tStart);
    %FIT_TIME(itr,1)=tEnd;
    %FIT_TIME(itr,2)=robust_fit_gt;
    %if(itr==50)
    %    break;
    %end
end
%fprintf(fid, 'RobustCopa = %f\n ',robust_fit_gt); 
%fclose(fid);


U=cell(K,1);
if(PARFOR_FLAG)
    parfor k=1:K
        U{k}=Rknew{k}*H;
    end
else
    for k=1:K
        U{k}=Rknew{k}*H;
    end
end

%plot U_k for a random subject w/o gap.
% subject_number=1;
% if(GAP==1)
%     M = csvread('gaps_per_subject.csv');
%     temp=M(subject_number,:);
%     temp(temp==0) = [];
%     for r=1:R
%       plot(temp,U{subject_number}(:,r))
%       hold on;
%     end

% else
%     for r=1:R
%       plot([1:size(U{subject_number},1)],U{subject_number}(:,r))
%       hold on;
%     end
% end
% ylabel("Value")
% xlabel("Number of observations")
% title("plot U{k} where k is 1")



end
