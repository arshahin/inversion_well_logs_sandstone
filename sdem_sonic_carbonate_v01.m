function [dtp,dts,vp,vs,ai,si,vr,pr,ksat,musat,den,logk,logmu]=...
         sdem_sonic_carbonate_v01...
         (phit,phii,phiv,phic,Lik,Limu,Lvk,Lvmu,...
         swt,sot,sgt,satype,kw,ko,kg,rhow,rhoo,rhog,ks,mus,rhos)


% % nc=size(phit,1);%%%%% number of data points (either core plugs or well data)


phit=phit+eps;


%%%%%%%%% properties of matrix(grain/inclusions)or soild rhos ks and mus are given  
% % % % fcalc=1.0;
% % % % fdol=0.0;
% % % % fsydr=0.0;
% % % % 
% % % % comp=fdol+fcalc+fsydr;
% % % % if comp~=1
% % % %     sprintf('The fraction of matrix components is greater than 1')
% % % %     logk=-1;
% % % %     logmu=-1;
% % % %     return;
% % % % end
% % % %     
% % % % rho_calc=2.71; k_calc=77;  mu_calc=32;
% % % % rho_dol=2.87;  k_dol=95;   mu_dol=45; 
% % % % rho_sydr=3.96; k_sydr=124;  mu_sydr=51;
% % % % 
% % % % rhos=fcalc.*rho_calc+fdol.*rho_dol+fsydr.*rho_sydr;  %%%%% matrix or soild density
% % % % 
% % % % k_ruess=sum(fcalc./k_calc+fdol./k_dol+fsydr./k_sydr).^-1;
% % % % mu_ruess=sum(fcalc./mu_calc+fdol./mu_dol+fsydr./mu_sydr).^-1;
% % % % 
% % % % k_voight=sum(fcalc.*k_calc+fdol.*k_dol+fsydr.*k_sydr);
% % % % mu_voight=sum(fcalc.*mu_calc+fdol.*mu_dol+fsydr.*mu_sydr);
% % % % 
% % % % %%hill average 
% % % % ks=0.5.*(k_ruess+k_voight); %%%% matrix bulk modulus 
% % % % mus=0.5.*(mu_ruess+mu_voight); %%%% matrix shear modulus 



%%%%%% bulk density
den=(1-phit).*rhos+phit.*(swt.*rhow+sot.*rhoo+sgt.*rhog);


%%%%%%%%% fluid properties 
if strcmp(satype,'uniform')==1; kf=(swt./kw+sot./ko+sgt./kg).^(-1);%Wood's mixing formula----uniform(homogenous)or fine scale saturation(Reuss lower bound)                               
elseif strcmp(satype,'patchy')==1; kf=(kw.*swt+ko.*sot+kg.*sgt);%Patchy saturation(Voight upper bound)hetrogenous
else
    sprintf('The saturation type is not in library')
    logk=-1;
    logmu=-1;
    return;
end

%%%%%%% determine mc (moduli at critical properties phic)

kc=((1-phic)./ks+phic./kw).^-1; %%%%% bulk modulus brine saturated sample at critical porosity
muc=10^-20; %%%%% shear modulus brine saturated sample at critical porosity


%%%%%%%%%% determine moduli for intergranular brine saturated data (mi)
a=(phic-phii)./phic;
bk=(ks-kc)./kc; 
bmu=(mus-muc)./muc;
ck=(ks-kc)./ks;
cmu=(mus-muc)./mus;

ki=kc.*(1+a.*bk.*Lik)./(1-a.*ck.*(1-Lik));  

mui=muc.*(1+a.*bmu.*Limu)./(1-a.*cmu.*(1-Limu)); 


%%%%%%%%%% determine moduli for brine saturated data (mbr)
%%%%%%%%%% Here Vugs come to play 
%%%%%%%%%% Reducing moduli using a linear function of porosity
kbr=ki-Lvk.*(phiv);
mubr=mui-Lvmu.*(phiv);

%%%%%% Gasmman from brine to partial saturations (determine msat)
gass=kbr./(ks-kbr)+kf./(phit.*(ks-kf))-kw./(phit.*(ks-kw));
ksat=(ks.*gass)./(1+gass);
musat=mubr;


%%%%%%%%%%%%% unit conversion
den=den.*1000;%unit conversion from gr/cc to Kg/m3
ksat=ksat*10.^9;%unit conversion from GPa to Pa for velocity calculation
musat=musat.*10.^9;%unit conversion from GPa to Pa for velocity calculation
vp=((ksat+4.*musat./3)./den).^0.5.*0.001;%km/s
vs=(musat./den).^0.5.*0.001;%km/s
ksat=ksat*10.^-9;%unit conversion from Pa to GPa for velocity calculation
musat=musat.*10.^-9;%unit conversion from Pa to GPa for velocity calculation
den=den./1000;%unit conversion from Kg/m3 to gr/cc
ai=vp.*den;
si=vs.*den;
vr=vp./vs;
pr=(0.5.*vr.^2-1)./(vr.^2-1);
dtp=10.^3./vp;
dts=10.^3./vs;

logk=1;
logmu=1;
end


