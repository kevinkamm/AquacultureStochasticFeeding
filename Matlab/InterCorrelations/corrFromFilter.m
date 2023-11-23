function [Rho,xRho,Lag]=corrFromFilter(ss_att,mode)
    if nargin<2
        mode='none';
    end
    switch mode
        %de-trending
        case 'differences'
            ss_att=diff(ss_att,1,1); 
        case 'returns'
            ss_att=(ss_att(2:end,:)-ss_att(1:end-1,:))./ss_att(1:end-1,:);
    end
    
    % Cross Correlation
    [rho12,lag12]=xcorr(ss_att(:,1),ss_att(:,2),'normalized');
    [rho13,lag13]=xcorr(ss_att(:,1),ss_att(:,3),'normalized');
    [rho14,lag14]=xcorr(ss_att(:,1),ss_att(:,4),'normalized');
    [rho23,lag23]=xcorr(ss_att(:,2),ss_att(:,3),'normalized');
    [rho24,lag24]=xcorr(ss_att(:,2),ss_att(:,4),'normalized');
    [rho34,lag34]=xcorr(ss_att(:,3),ss_att(:,4),'normalized');
    
    xRho=repmat(eye(4,4),[1,1,length(rho12)]);
    Lag=zeros(4,4,length(rho12));
    xRho(1,2,:)=rho12;xRho(2,1,:)=rho12;
    xRho(1,3,:)=rho13;xRho(3,1,:)=rho13;
    xRho(1,4,:)=rho14;xRho(4,1,:)=rho14;
    xRho(2,3,:)=rho23;xRho(3,2,:)=rho23;
    xRho(2,4,:)=rho24;xRho(4,2,:)=rho24;
    xRho(3,4,:)=rho34;xRho(4,3,:)=rho34;
    
    Lag(1,2,:)=lag12;Lag(2,1,:)=lag12;
    Lag(1,3,:)=lag13;Lag(3,1,:)=lag13;
    Lag(1,4,:)=lag14;Lag(4,1,:)=lag14;
    Lag(2,3,:)=lag23;Lag(3,2,:)=lag23;
    Lag(2,4,:)=lag24;Lag(4,2,:)=lag24;
    Lag(3,4,:)=lag34;Lag(4,3,:)=lag34;
    
    % Correlation
    rho12=corr(ss_att(:,1),ss_att(:,2));
    rho13=corr(ss_att(:,1),ss_att(:,3));
    rho14=corr(ss_att(:,1),ss_att(:,4));
    rho23=corr(ss_att(:,2),ss_att(:,3));
    rho24=corr(ss_att(:,2),ss_att(:,4));
    rho34=corr(ss_att(:,3),ss_att(:,4));
    
    Rho=eye(4,4);
    Rho(1,2)=rho12;Rho(2,1)=rho12;
    Rho(1,3)=rho13;Rho(3,1)=rho13;
    Rho(1,4)=rho14;Rho(4,1)=rho14;
    Rho(2,3)=rho23;Rho(3,2)=rho23;
    Rho(2,4)=rho24;Rho(4,2)=rho24;
    Rho(3,4)=rho34;Rho(4,3)=rho34;
end
