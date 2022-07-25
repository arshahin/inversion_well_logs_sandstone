function [objfun,logr,logk,logmu]=objfun_gammaray_sdem_density_sonic_resist_sandstone_v01...
         (gr,resist,dtp,dts,den,phi_sh,phit,phic,phii,phim,swi,swm,lamdai,lamdam,tmpc,salppk,wetcase,...
          Lik,Limu,Lmk,Lmmu,...
          GR_clay,csh,GRma,...
          swt,sot,sgt,satype,kw,ko,kg,rhow,rhoo,rhog,rho_clay,k_clay,mu_clay,rho_ma,k_ma,mu_ma)

%%% Cx is the salinity in Kppm

%%%%%%%%%%% compute resistivity

[resist_sim,ffac_sim,rw_sim,logr]=sdem_resist_sandstone(phit,phii,phim,swi,swm,lamdai,lamdam,tmpc,salppk,wetcase);

%%%%%%%%%%% compute velocities and density


[dtp_sim,dts_sim,vp_sim,vs_sim,ai_sim,si_sim,vr_sim,pr_sim,ksat_sim,musat_sim,den_sim,logk,logmu]=...
         sdem_sonic_sandstone_v01...
         (phi_sh,phit,phii,phim,phic,Lik,Limu,Lmk,Lmmu,...
         swi,swt,sot,sgt,satype,kw,ko,kg,rhow,rhoo,rhog,rho_clay,k_clay,mu_clay,rho_ma,k_ma,mu_ma);
         
     
[gr_sim,logr1]=GR_Sand_v01(phi_sh,rhow,den,GR_clay,csh,rho_clay,rho_ma,GRma,phii);     
     
     
%%%%%%%%%%% compute normalized error function 
sum0=sum((resist_sim-resist).^2);  sum00=sum((resist_sim+resist).^2);
sum1=sum((dtp_sim-dtp).^2);        sum11=sum((dtp_sim+dtp).^2);
sum2=sum((dts_sim-dts).^2);        sum22=sum((dts_sim+dts).^2);
sum3=sum((den_sim-den).^2);        sum33=sum((den_sim+den).^2);
sum4=sum((gr_sim-gr).^2);          sum44=sum((gr_sim+gr).^2);

objfun=(2*sum0)/(sum00+eps)+(2*sum1)/(sum11+eps)+(2*sum2)/(sum22+eps)+(2*sum3)/(sum33+eps)+(2*sum4)/(sum44+eps);
objfun=objfun/5; %%% normalize error function

end













