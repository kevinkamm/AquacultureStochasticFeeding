function [params,negLogLikelihood]=EMschwartz(yData,r,dt,panel,varEpsMode)
    % fminsearch does not work
    % fmincon changes with different bounds
    % fminunc with exp(params)->positive and rho=tanh(params)\in [0,1]
    %   similar to fmincon with reasonable bounds. During optimization many
    %   warnings with rcond
    % lsqnonlin does not make sense here
    % params=[mu,sigma11,sigma12,kappa1,alpha1,lambda1,rho]

    [nObservations,nContracts] = size(yData);
    
    switch varEpsMode
        case 'single'
            n=1;
        case 'diag'
            n=nContracts;
        case 'tridiag'
            n=nContracts*(nContracts+1)/2;
        otherwise
            error('Unrecognized mode for variance of epsilon')
    end

    % lb=zeros(1,7+n);
    lb=.01.*ones(1,7+n);lb(8:end)=0.001;
    lb(5)=-5;%alpha
    lb(7)=-1;
    ub=5.*ones(1,7+n);
    ub(7:end)=1;%ub(5)=1;
    x0=.1.*ones(1,7+n);
    x0(7)=.5;

    optionsF = optimoptions('fmincon','Display','off','UseParallel',false,'StepTolerance',1e-12,'OptimalityTolerance',1e-9,'MaxFunctionEvaluations',1e6);
    optionsU = optimoptions('fminunc','Display','off','UseParallel',false,'StepTolerance',1e-12,'OptimalityTolerance',1e-9,'MaxFunctionEvaluations',1e6);
    % optionsGA = optimoptions('ga','Display','iter','UseParallel',true);
    
    % x0 = ga(@(params) objective(params,yData,r,dt,panel),length(lb),[],[],[],[],lb,ub./10,@(params)mycon(params,yData,r,dt,panel),[],optionsGA);
    % x0 = ga(@(params) objective(params,yData,r,dt,panel),length(lb),[],[],[],[],lb,ub,[],[],optionsGA);
    [params,negLogLikelihood]= fmincon(...
        @(params) objective(params,yData,r,dt,panel,varEpsMode),...
        x0,[],[],[],[],lb,ub,[],optionsF);
    % [params,negLogLikelihood]= fminunc(...
    %     @(params) objective(params,yData,r,dt,panel,varEpsMode),...
    %     x0,optionsU);

end
function logL=objective(params,yData,r,dt,panel,varEpsMode)
    [P0,a0,c,T,varEta,dtau,Ztau,varEps,yData2]=schwartzSSM(params,yData,r,dt,panel,varEpsMode);
    logL=kalmanFilter(P0,a0,c,T,varEta,dtau,Ztau,varEps,yData2);
end
function [logL,ceq] = mycon(params,yData,r,dt,panel,varEpsMode)
    logL=-objective(params,yData,r,dt,panel,varEpsMode);
    ceq=[];
end