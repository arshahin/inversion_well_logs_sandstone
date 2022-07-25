function [dtp,dts,vp,vs,ai,si,vr,pr,ksat,musat,den,logk,logmu]=...
         sdem_sonic_sandstone_v01...
         (phi_sh,phit,phii,phim,phic,Lik,Limu,Lmk,Lmmu,...
         swi,swt,sot,sgt,satype,kw,k_oil,k_gas,rhow,rho_oil,rho_gas,rho_clay,k_clay,mu_clay,rho_ma,k_ma,mu_ma)

csh=phim./phi_sh;     
cma=1-phii-csh;
%%%%%% bulk density
rho_sh=phi_sh.*rhow+(1-phi_sh).*rho_clay; %%% shale density
den=cma.*rho_ma+phii.*(swi.*rhow+sot.*rho_oil+sgt.*rho_gas)+rho_sh.*csh;

%%%%%%%%%%% bulk moduli (integrate matrix and clay using Hill average 
phi0=eps;%%%% indication of zeros porosity at matrix point 
csum=csh+cma+phi0;
%%%%%% upper bound 

k_voit=(csh.*k_clay+cma.*k_ma)./csum;
mu_voit=(csh.*mu_clay+cma.*mu_ma)./csum;
%%%%%% lower bound 
k_ruess=((csh./k_clay+cma./k_ma)./csum).^-1;
mu_ruess=((csh./mu_clay+cma./mu_ma)./csum).^-1;
%%%%%% Hill average bound 
khill_ma=0.5.*(k_ruess+k_voit); %%%% effective bulk modulus of clay and non-clay minerals  
muhill_ma=0.5.*(mu_ruess+mu_voit); %%%% effective shear modulus of clay and non-clay minerals 

%%%%%%%%% fluid properties 
if strcmp(satype,'uniform')==1; kf=(swt./kw+sot./k_oil+sgt./k_gas).^(-1);%Wood's mixing formula----uniform(homogenous)or fine scale saturation(Reuss lower bound)                               
elseif strcmp(satype,'patchy')==1; kf=(kw.*swt+k_oil.*sot+k_gas.*sgt);%Patchy saturation(Voight upper bound)hetrogenous
else
    sprintf('The saturation type is not in library')
    logk=-1;
    logmu=-1;
    return;
end

%%%%%%% determine mc (moduli at critical properties phic)

kc=((1-phic)./khill_ma+phic./kw).^-1; %%%%% bulk modulus brine saturated sample at critical porosity
muc=eps; %%%%% shear modulus brine saturated sample at critical porosity


%%%%%%%%%% determine moduli for intergranular brine saturated data (mi)
a=(phic-phii)./phic;
bk=(khill_ma-kc)./kc; 
bmu=(muhill_ma-muc)./muc;
ck=(khill_ma-kc)./khill_ma;
cmu=(muhill_ma-muc)./muhill_ma;

ki=kc.*(1+a.*bk.*Lik)./(1-a.*ck.*(1-Lik));  

mui=muc.*(1+a.*bmu.*Limu)./(1-a.*cmu.*(1-Limu)); 


%%%%%%%%%% determine moduli for brine saturated data (mbr)
%%%%%%%%%% Here micro-porosity come to play 
%%%%%%%%%% Reducing moduli using a linear function of porosity
kbr=ki-Lmk.*(phim);
mubr=mui-Lmmu.*(phim);

%%%%%% Gasmman from brine to partial saturations (determine msat)
gass=kbr./(khill_ma-kbr)+kf./(phit.*(khill_ma-kf))-kw./(phit.*(khill_ma-kw));
ksat=(khill_ma.*gass)./(1+gass);
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


