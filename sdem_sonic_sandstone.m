function [dtp,dts,vp,vs,ai,si,vr,pr,ksat,musat,den,rho_sh,logk,logmu]=...
         sdem_sonic_sandstone...
         (phi_sh,phit,phii,phim,phic,Lik,Limu,Lmk,Lmmu,...
         swi,swt,sot,sgt,satype,kw,ko,kg,rhow,rhoo,rhog,csh)


% % nc=size(phit,1);%%%%% number of data points (either core plugs or well data)


phit=phit+eps;

%%%%%%%%% properties of matrix(grain/inclusions)or soild rhos ks and mus 
fclay=csh;
fq=1-phii-fclay;  %%%% quartz volume concentration 
   
rho_q=2.65; k_q=37;  mu_q=44;
rho_clay=2.41;  k_clay=21.0;   mu_clay=7.0; 

rhos=fq.*rho_q+fclay.*rho_clay;  %%%%% matrix or soild density



cdc=1-phi_sh; %%%% volume concentration of dry clay in shale
rhow_sh=rhow;   %%%% salinity 120kppm
rho_sh=rhow_sh.*phi_sh+cdc.*rho_clay;


k_ruess=sum(fq./k_q+fclay./k_clay).^-1;
mu_ruess=sum(fq./mu_q+fclay./mu_clay).^-1;

k_voight=sum(fq.*k_q+fclay.*k_clay);
mu_voight=sum(fq.*mu_q+fclay.*mu_clay);

%%hill average 
ks=0.5.*(k_ruess+k_voight); %%%% matrix bulk modulus 
mus=0.5.*(mu_ruess+mu_voight); %%%% matrix shear modulus 



%%%%%% bulk density
den=(1-phii-csh).*rhos+phii.*(swi.*rhow+sot.*rhoo+sgt.*rhog)+rho_sh.*csh;

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
kbr=ki-Lmk.*(phim);
mubr=mui-Lmmu.*(phim);

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


