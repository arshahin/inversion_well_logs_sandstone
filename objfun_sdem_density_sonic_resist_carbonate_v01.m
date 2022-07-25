function [objfun,logr,logk,logmu]=objfun_sdem_density_sonic_resist_carbonate_v01...
         (resist,dtp,dts,den,phit,phic,phii,phiv,swi,swv,lamdai,lamdav,tmpc,salppk,wetcase,...
          Lik,Limu,Lvk,Lvmu,...
          swt,sot,sgt,satype,kw,ko,kg,rhow,rhoo,rhog,ks,mus,rhos)

%%% Cx is the salinity in Kppm

%%%%%%%%%%% compute resistivity

[resist_sim,ffac_sim,rw_sim,logr]=sdem_resist_carbonate(phit,phii,phiv,swi,swv,lamdai,lamdav,tmpc,salppk,wetcase);

%%%%%%%%%%% compute velocities and density


[dtp_sim,dts_sim,vp_sim,vs_sim,ai_sim,si_sim,vr_sim,pr_sim,ksat_sim,musat_sim,den_sim,logk,logmu]=...
         sdem_sonic_carbonate_v01...
         (phit,phii,phiv,phic,Lik,Limu,Lvk,Lvmu,...
         swt,sot,sgt,satype,kw,ko,kg,rhow,rhoo,rhog,ks,mus,rhos);

   
%%%%%%%%%%% compute normalized error function 
sum0=sum((resist_sim-resist).^2);  sum00=sum((resist_sim+resist).^2);
sum1=sum((dtp_sim-dtp).^2);        sum11=sum((dtp_sim+dtp).^2);
sum2=sum((dts_sim-dts).^2);        sum22=sum((dts_sim+dts).^2);
sum3=sum((den_sim-den).^2);        sum33=sum((den_sim+den).^2);
objfun=(2*sum0)/(sum00+eps)+(2*sum1)/(sum11+eps)+(2*sum2)/(sum22+eps)+(2*sum3)/(sum33+eps);
objfun=objfun/4; %%% normalize error function


%%%%%%%%%%%% use density porosity assumiung water-filled and known matrix
%%%%%%%%%%%% as a proir data 
% % % % rho_ma=2.71; %%%% come from sdem_sonic_carbonate.m
% % % % rhow=1.05;
% % % % phitw_den=(den-rho_ma)./(rhow-rho_ma);
% % % % sum0=sum((resist_sim-resist).^2);  sum00=sum((resist_sim+resist).^2);
% % % % sum1=sum((dtp_sim-dtp).^2);        sum11=sum((dtp_sim+dtp).^2);
% % % % sum2=sum((dts_sim-dts).^2);        sum22=sum((dts_sim+dts).^2);
% % % % sum3=sum((den_sim-den).^2);        sum33=sum((den_sim+den).^2);
% % % % sum4=sum((phitw_den-phit).^2);        sum44=sum((phitw_den+phit).^2);
% % % % w0=0.15;
% % % % w1=0.15;
% % % % w2=0.15;
% % % % w3=0.15;
% % % % w4=1-(w0+w1+w2+w3);
% % % % objfun=w0.*(2*sum0)/(sum00+eps)+w1.*(2*sum1)/(sum11+eps)+w2.*(2*sum2)/(sum22+eps)+w3.*(2*sum3)/(sum33+eps)+w4.*(2*sum4)/(sum44+eps);
% % % % objfun=objfun/5; %%% normalize error function


end













