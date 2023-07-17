function [logL,varargout]=kalmanFilter(P0,a0,c,T,varEta,d,Z,varEps,y)
%%
% THE TRANSITION EQUATION 
% a(t)=c+T*a(t-1)+eta(t)   eta~N(0,Q), Q=varEta
%
% THE MEASUREMENT EQUATION 
% y(t,tau)=d(tau)+Z(tau)*a(t)+epsilon(t)  epsilon~N(0,H), H=varEps


Ptt = P0;
att = a0; 
[nobs,N]=size(y);
save_att=zeros(size(y,1),length(a0));
save_dFtt_1=zeros(size(y,1),1);
save_vFv=zeros(size(y,1),1);
for t = 1:size(y,1)
    Ptt_1   = T*Ptt*T'+varEta;
    Ftt_1   = Z*Ptt_1*Z'+varEps;
    dFtt_1  = det(Ftt_1);
    invFtt_1 = inv(Ftt_1);

    att_1   = T*att + c;
    yt      = y(t,:)';
    ytt_1   = Z*att_1+d;
    vt      = yt-ytt_1;

    att = att_1 + Ptt_1*Z'*invFtt_1*(vt);
    Ptt = Ptt_1 - Ptt_1*Z'*invFtt_1*Z*Ptt_1;
    
    % ytt = Z*att+d;
    % vtt  = yt-ytt;

    % save_vtt(t,:) = vtt';
    % save_vt(t,:)    = (vt)';
    save_att(t,:)   = att';
    % save_Ptt_1(t,:) = [Ptt_1(1,1), Ptt_1(1,2), Ptt_1(2,1), Ptt_1(2,2)]; 
    % save_Ptt(t,:)   = [Ptt(1,1), Ptt(1,2), Ptt(2,1), Ptt(2,2)];

    save_dFtt_1(t,:)= dFtt_1;
    save_vFv(t,:)   = vt'*invFtt_1*vt;
    
end

save_dFtt_1(save_dFtt_1<0)=0;%prevent complex number for unreasonable parameters
logL = -(-(N*nobs/2)*log(2*pi)-0.5*sum(log(save_dFtt_1))-0.5*sum(save_vFv));

if nargout>1
    varargout{1}=save_att;
end
end