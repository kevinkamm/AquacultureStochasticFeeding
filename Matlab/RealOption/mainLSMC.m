clear all; close all; fclose('all'); seed=0;rng(seed);
pool=gcp('nocreate');
if isempty(pool)
%     pool=parpool('local'); % multiprocessing
    pool=parpool('threads'); % multithreading
end
availableGPUs = gpuDeviceCount('available');
if availableGPUs > 0
    gpuDevice([]); % clears GPU
    gpuDevice(1); % selects first GPU, change for multiple with spmd
end

%% Figure options
backgroundColor='w';
textColor='k';

%% Parameters
% Monte Carlo parameters
simFactor=10;
N=3*24*simFactor; % points in time grid
M=100000; % number of simulations
antithetic=false;
useSVM=false;

% Parameters for fish farming, see Table 4 Ewald2017
T=3; % time horizon
m=.1; % mortality rate
cr=1.1; % conversion rate
n0=10000; % number of recruits
hc=3; % variable harvesting cost per kg (NOK/kg)
fc=10; % variable feeding cost per kg per year (NOK/kg)
wInf=6; % asymptotic weight (kg)

% Bertalanffyo's growth function, see Footnote 20 Ewald2017
a=1.113;
b=1.097;
c=1.43;

r=0.0303;

gamma=0.0;

% Salmon
% mu, sigma1, sigma2, kappa, alpha, lambda, rho, delta0, P0
salmonParam=[0.428409, 0.998737, 1.22314, 0.862705, 0.0100041, 0.322524, 0.998149, 0.0303, 45.2295];
fc=salmonParam(end).*.5.*.25;
salmonParam(end)=salmonParam(end).*.5+fc+hc;

% Soy
% mu, sigma1, sigma2, kappa, alpha, lambda, rho, delta0, P0
soyParam=[1,1.1636    0.22166    0.27045    2.2616    0.78656    0.53143, r, 1];

for stochFeeding=[false, true]
    rng(seed);
    %% Simulate processes
    %% Salmon
    % Brownian motions, directly under Q
    ticBM=tic;
    [W1,W2,t]=brownianMotions(T,N,M,salmonParam(7),'antithetic',antithetic);
    ctimeBM=toc(ticBM);
    fprintf('Elapsed time BM %g s.\n',ctimeBM)
    
    % Schwartz 2-factor model under Q
    ticSchwartz=tic;
    [salmonP,salmonDelta]=schwartzTwoFactor(salmonParam(8),salmonParam(9),r,salmonParam(2),salmonParam(3),salmonParam(4),salmonParam(5),salmonParam(6),W1,W2,t);
    ctimeSchwartz=toc(ticSchwartz);
    fprintf('Elapsed time Schwartz %g s.\n',ctimeSchwartz)
    
    
    %% Soy
    % Brownian motions, directly under Q
    ticBM=tic;
    [W1,W2,t]=brownianMotions(T,N,M,soyParam(7),'antithetic',antithetic);
    ctimeBM=toc(ticBM);
    fprintf('Elapsed time BM %g s.\n',ctimeBM)
    
    % Schwartz 2-factor model under Q
    ticSchwartz=tic;
    [soyP,soyDelta]=schwartzTwoFactor(soyParam(8),soyParam(9),r,soyParam(2),soyParam(3),soyParam(4),soyParam(5),soyParam(6),W1,W2,t);
    soyP=soyP./soyParam(9);
    ctimeSchwartz=toc(ticSchwartz);
    fprintf('Elapsed time Schwartz %g s.\n',ctimeSchwartz)


    %% Anticipative Stopping
    dt=T./(N-1);
    % Bertalanffy’s growth function
    wt=wInf.*(a-b.*exp(-c.*t)).^3;
    dwt=diff(wt,1,1)./dt;
    % number of fish
    nt=n0.*exp(-m.*t);
    
    % total biomass (kg)
    Xt=nt.*wt;
    
    % harvesting cost
    CH=Xt.*hc;
    
    % feeding cost
    CF=zeros(size(Xt));
    CF(2:end)=dwt.*nt(2:end).*cr.*fc;
    if stochFeeding
        CF=CF.*soyP;
    else
        CF=CF.*mean(soyP,2);
        % CF=CF;
        % CF=CF.*soyP;
    end
    CF=CF.^(1-gamma)./(1-gamma);
    
    ICFdt=zeros(size(CF));
    ICFdt(2:end,:)=exp(-r.*t(2:end)).*CF(2:end,:);
    ICFdt=cumsum(ICFdt,1).*dt;
    
    
    %% sim points -> eval points
    indCoarse = unique([1:simFactor:N,N]);
    tCoarse = t(indCoarse);
    salmonP=salmonP(indCoarse,:);
    salmonDelta=salmonDelta(indCoarse,:);
    soyP=soyP(indCoarse,:);
    soyDelta=soyDelta(indCoarse,:);
    Xt=Xt(indCoarse,:);
    CH=CH(indCoarse,:);
    ICFdt=ICFdt(indCoarse,:);
    CF=CF(indCoarse,:);
    
    %% LSMC
    VH=(salmonP.*Xt-CH).^(1-gamma)./(1-gamma);
    VC=zeros(size(VH));
    exercise=zeros(size(VH));
    V=zeros(size(VH));
    
    n=size(VH,1);
    dtn=T/(n-1);

    if stochFeeding
        exerciseRegion=cell(n,4);
        exerciseRegionOpt=cell(n,4);
        contRegion=cell(n,4);
        contRegionOpt=cell(n,4);
    else
        exerciseRegion=cell(n,2);
        exerciseRegionOpt=cell(n,2);
        contRegion=cell(n,2);
        contRegionOpt=cell(n,2);
    end
    
    V(end,:)=VH(end,:);
    exercise(end,:)=1;
    
    ticLSMC=tic;
    for ti=n-1:-1:2
        if stochFeeding
            VC(ti,:)=-CF(ti,:).*dtn + exp(-r.*dtn).*basis([salmonP(ti,:);salmonDelta(ti,:);soyP(ti,:);soyDelta(ti,:)],V(ti+1,:)')';
        else
            VC(ti,:)=-CF(ti,:).*dtn + exp(-r.*dtn).*basis([salmonP(ti,:);salmonDelta(ti,:)],V(ti+1,:)')';
        end

        ind=VC(ti,:)<=VH(ti,:);
        exercise(ti,:)=ind;
        exerciseRegion{ti,1}=salmonP(ti,ind);
        exerciseRegion{ti,2}=salmonDelta(ti,ind);
        contRegion{ti,1}=salmonP(ti,~ind);
        contRegion{ti,2}=salmonDelta(ti,~ind);
        if stochFeeding
            exerciseRegion{ti,3}=soyP(ti,ind);
            exerciseRegion{ti,4}=soyDelta(ti,ind);
            contRegion{ti,3}=soyP(ti,~ind);
            contRegion{ti,4}=soyDelta(ti,~ind);
        end
        % Glassermann p. 461
        V(ti,:)=exp(-r.*dtn).*V(ti+1,:)-CF(ti,:).*dtn; %Longstaff-Schwartz
        % V(ti,:)=VC(ti,:); %Tsitsiklis and Van Roy
        V(ti,ind)=VH(ti,ind);
    end
    ctimeLSMC=toc(ticLSMC);
    fprintf('Elapsed time LSMC in case stoch=%d: %g s\n',stochFeeding,ctimeLSMC)
    value=exp(-r.*dt).*mean(V(2,:),2);
    
    tau=zeros(M,1);
    for wi=1:M
        tau(wi)=find(exercise(:,wi),1,'first');
    end
    if stochFeeding
        tauStoch=tau;
    else
        tauDeterm=tau;
    end
    tau=tCoarse(tau);
    

    decisionSVMdeterm=cell(n,1);
    decisionSVMstoch=cell(n,1);
    remainingPaths=cell(n,1);
    paths=1:M;
    totalPaths=1:M;
    ticTau=tic;
    for ti=1:1:n
        ind=tau==tCoarse(ti);
        paths=setdiff(paths,totalPaths(ind));
        remainingPaths{ti}=paths;
        exerciseRegionOpt{ti,1}=salmonP(ti,ind);
        exerciseRegionOpt{ti,2}=salmonDelta(ti,ind);
        contRegionOpt{ti,1}=salmonP(ti,paths);
        contRegionOpt{ti,2}=salmonDelta(ti,paths);
        if stochFeeding
            exerciseRegionOpt{ti,3}=soyP(ti,ind);
            exerciseRegionOpt{ti,4}=soyDelta(ti,ind);
            contRegionOpt{ti,3}=soyP(ti,paths);
            contRegionOpt{ti,4}=soyDelta(ti,paths);
           
        end
    end

    ctimeTau=toc(ticTau);
    fprintf('Elapsed time for tau %g \n',ctimeTau);
    
    % %%
    objective=exp(-r.*tCoarse).*(salmonP.*Xt-CH).^(1-gamma)./(1-gamma)-ICFdt;
    
    [anticipativeMax,I]=max(objective,[],1);
    stoppingTime=tCoarse(I);
    
    if stochFeeding
        % title('Objective stoch feeding')
        fprintf('Mean anticipative max (stochastic feeding):\n\t %g at mean time %g\n',mean(anticipativeMax,2),mean(stoppingTime,1));
        fprintf('Mean LSMC max (stochastic feeding):\n\t %g at mean time %g\n',value,mean(tau));
    else
        % title('Objective determ feeding')
        fprintf('Mean anticipative max (determ feeding):\n\t %g at mean time %g\n',mean(anticipativeMax,2),mean(stoppingTime,1));
        fprintf('Mean LSMC max (determ feeding):\n\t %g at mean time %g\n',value,mean(tau));
    end

    %% Exercise Region
    % figure();tiledlayout(4,4);
    % tn=find(tCoarse>=mean(tau),1,'first');
    % for j=-5:1:5
    %     if tn+j<=n
    %         ax=nexttile;hold on;
    %         scatter(exerciseRegion{tn+j,1},exerciseRegion{tn+j,2},'blue','.');hold on;
    %         scatter(contRegion{tn+j,1},contRegion{tn+j,2},'cyan','.');hold on;
    %         scatter(exerciseRegionOpt{tn+j,1},exerciseRegionOpt{tn+j,2},'red','.');
    %         title(ax,sprintf('at time t=%g, i=%d',tCoarse(tn+j),tn+j))
    %     end
    % end
    figure();tiledlayout(4,4);
    tn=find(tCoarse>=mean(tau),1,'first');
    for j=-5:1:5
        if tn+j<=n
            ax=nexttile;hold on;
            scatter(contRegionOpt{tn+j,1},contRegionOpt{tn+j,2},'yellow','.');hold on;
            scatter(exerciseRegionOpt{tn+j,1},exerciseRegionOpt{tn+j,2},'red','.');
            title(ax,sprintf('at time t=%g, i=%d',tCoarse(tn+j),tn+j))
        end
    end
end

vDeterm=zeros(M,1);
vStoch=zeros(M,1);
for wi=1:M
    vDeterm(wi)=objective(tauDeterm(wi),wi);
    vStoch(wi)=objective(tauStoch(wi),wi);
end
fprintf('Determ stopping: Value %g at mean %g\n',mean(vDeterm),mean(tCoarse(tauDeterm)));
fprintf('Stoch stopping: Value %g at mean %g\n',mean(vStoch),mean(tCoarse(tauStoch)));
fprintf('Stoch revenue = %g * Determ Revenue\n',mean(vStoch)/mean(vDeterm));