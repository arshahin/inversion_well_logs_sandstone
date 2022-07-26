%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% version 1
%%%%%%%% sandstone model for log and seismic inversion
%%%%%%%% sandstone containing caly and quartz with two pore system
%%%%%% 1- microporosity(isolated) porosity associated with Micro : phim
%%%%%%% 2- intergranular porosity (matrix porosity) phii
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Cleaning 

close all
clear
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% VFSA varying parameters
% % % % pwn_data=0; %% white noise normal percentage (added to data)

reproduce=1; % 0 for random, 1 for reproducable results

if reproduce==1
    rand('seed',31415927)
    randn('seed',3111113)
end

nruns=1; %%%%% Number of runs(cases) 
niter=2000; % number of iteration per run
error_tresh=0.0001; %%%%error threshould to stop VFSA
%%%% Parameter definition for basic VFSA
nmov=3;

%%%%%%%%%%%%%%%%% # of iteration after which local errors go into effect
%%%%%%%%%%%%%%%%%%%iter_loc=1 mean from iter=1 local update go into effect
iter_loc=round(niter/20);
freq_loc=5;
%%%% weighting for local and global errors
wloc=0.75; wglob=1-wloc;
reannealing=2;  %%%% zero for no reannealing, one for Ingber sensitivity method and two for Armando method
iter_reann=iter_loc+50; %%% works only if reannealing~=0 (iter# after which local temp update goes into effect iter_reann=1 mean from iter=1 local temp go into effect )
freq_reann=50; %%% works only if reannealing~=0 (after tests 50 and 100 is optimum for model matching)
gama=1.20;
temp0=1.5;%%%common initial temp for global and local
decay=0.999;%%%common decay for global and local temp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOADING EDITTED LAS FILE
filename_slb01='slb_log01_edt.txt'; filename_slb02='slb_log02_edt.txt'; filename_slb03='slb_log03_edt.txt'; filename_slb04='slb_log04_edt.txt';
log_slb01=load(filename_slb01);     log_slb02=load(filename_slb02);   log_slb03=load(filename_slb03);   log_slb04=load(filename_slb04);
%%%%%%%%%%%%%%%%%%%%% ASSIGNING VARIBALE NAMEES, SIMPLY COPY AND PASTE FROM LAS FILE HEADER
NDPHI=log_slb01(:,9)./100;  
NPHI=log_slb01(:,10)./100;  
DPHI=log_slb01(:,11)./100;  
m=log_slb03(:,3);  
phit=NDPHI+0.07;
nz1=length(phit);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Creating Synthetic Logs 
swirr=0.10; %%% irreducable water saturation assiciated with CBW and Capilary bound water 
sores=0.03; %%% residual oil saturation

phi_sh=0.20; %%%%%% shale (Inllite with density 2.52 GR=275, Vp=4.32 and Vs=2.54)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% random generation of shale volume
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% concentration (Csh)
csh=zeros(size(phit));
for jz=1:nz1
    csh(jz)=monte_carlo_sample_matlab(0,0.5);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% smoothing csh
phim=csh.*phi_sh; %%%% isolated Micro-porosity (micro-porosity)

phii=phit-phim; %%%%%% intergranular porosity 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% making swi and swm
swi=zeros(nz1,1); %%swm=zeros(np,1);
%%%%%%%%%%%%%%%%%%%%%%%%% zone A oil zone 
for jz=1:74
    swi(jz)=swirr;
end

%%%%%%%%%%%%%%%%%%%%%%%%% zone B (transition zone) 
for jz=75:95
    slop_swi=(1-sores-swirr)/(95-75);
    intersect_swi=swirr-slop_swi*75;
    swi(jz)=slop_swi*jz+intersect_swi;
end

%%%%%%%%%%%%%%%%%%%%%%%%% zone C (brine)   

for jz=96:nz1
    swi(jz)=1-sores;
end

%%%%%%%%%%%%%%%%%%%%%% adding shale later up and down 
nsh=30;
ksh=25;
sma=eps;

csh=[ones(nsh,1);csh;ones(nsh,1)];
phii=[zeros(nsh,1)+sma;phii;zeros(nsh,1)+sma];
swi=[zeros(nsh,1)+sma;swi;zeros(nsh,1)+sma];

nz=length(phii);
% % dz=0.5; %%%% feet
% % dz=dz*0.3048; %% meter 
dz=1; 
zlas=(nz-1)*dz;
DEPTH=0:dz:zlas;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%linear interpolation at edges 

%%%%%%%%% csh
slop_csh1=(csh(nsh+1)-csh(ksh))/(nsh+1-ksh);
intersect_csh1=csh(ksh)-slop_csh1*ksh;
for jz=ksh:nsh    
    csh(jz)=jz*slop_csh1+intersect_csh1;   
end

slop_csh2=(csh(nz-nsh-1)-csh(nz-ksh))/(nz-nsh-1-nz+ksh);
intersect_csh2=csh(nz-ksh)-slop_csh2*(nz-ksh);

for jz=(nz-nsh):(nz-ksh)    
    csh(jz)=jz*slop_csh2+intersect_csh2;   
end
%%%%%%%%%%phii
slop_phii1=(phii(nsh+1)-phii(ksh))/(nsh+1-ksh);
intersect_phii1=phii(ksh)-slop_phii1*ksh;
for jz=ksh:nsh    
    phii(jz)=jz*slop_phii1+intersect_phii1;   
end

slop_phii2=(phii(nz-nsh-1)-phii(nz-ksh))/(nz-nsh-1-nz+ksh);
intersect_phii2=phii(nz-ksh)-slop_phii2*(nz-ksh);

for jz=(nz-nsh):(nz-ksh)    
    phii(jz)=jz*slop_phii2+intersect_phii2;   
end
%%%%%%%%% swi
slop_swi1=(swi(nsh+1)-swi(ksh))/(nsh+1-ksh);
intersect_swi1=swi(ksh)-slop_swi1*ksh;
for jz=ksh:nsh    
    swi(jz)=jz*slop_swi1+intersect_swi1;   
end

slop_swi2=(swi(nz-nsh-1)-swi(nz-ksh))/(nz-nsh-1-nz+ksh);
intersect_swi2=swi(nz-ksh)-slop_swi2*(nz-ksh);

for jz=(nz-nsh):(nz-ksh)    
    swi(jz)=jz*slop_swi2+intersect_swi2;   
end
%%%%%%%%%%%%%%%%%%%%%%%%


csh_smth1=smth3_handy(csh);
csh_smth2=smth3_handy(csh_smth1);
csh_smth3=smth3_handy(csh_smth2);
csh_smth4=smth3_handy(csh_smth3);
csh=csh_smth4';

phii_smth1=smth3_handy(phii);
phii_smth2=smth3_handy(phii_smth1);
phii_smth3=smth3_handy(phii_smth2);
phii_smth4=smth3_handy(phii_smth3);
phii_smth5=smth3_handy(phii_smth4);
phii=phii_smth5';

swi_smth1=smth3_handy(swi);
swi_smth2=smth3_handy(swi_smth1);
swi_smth3=smth3_handy(swi_smth2);
swi_smth4=smth3_handy(swi_smth3);
swi=swi_smth4';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% adding random noise 
pwn_add=0.01;                     % Standard Deviation Of Random Noise Vector
noisy=pwn_add*randn(nz,1);         % Create Normally-Distributed Random Noise Vector With Standard Deviation = �sd�
csh_noisy=csh+noisy;      % Add Noise To �signal
swi_noisy=swi+noisy;
phii_noisy=phii+noisy;
% % % % 
% % % % plot(csh_noisy,'LineWidth',2)
% % % % hold on
% % % % plot(csh,'r','LineWidth',2)
% % % % 
csh=csh_noisy;
phii=phii_noisy;
swi=swi_noisy;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% QC bounds 
for jz=1:nz
    if csh(jz)>1.0
        csh(jz)=1.0;
    end
    if swi(jz)<0.0
        swi(jz)=0.001;
    end
    if swi(jz)>1.0
        swi(jz)=1.0;
    end
    if phii(jz)<0.0
        phii(jz)=0.001;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phim=csh.*phi_sh; %%%% isolated Micro-porosity (micro-porosity)

phit=phii+phim; %%%%%% total porosity 

swm=(csh.*phi_sh)./phit; %%%%%%%%%% bound water saturation 

phiiw=swi.*phii;
swt=phiiw./phit+swm;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% QC bounds 
for jz=1:nz
    if phim(jz)>phi_sh
        phim(jz)=phi_sh;
    end
    if phim(jz)<0.0
        phim(jz)=0.001;
    end
    if phit(jz)<0.0
        phit(jz)=0.001;
    end
        
    if swm(jz)<0.0
        swm(jz)=0.001;
    end
    if swm(jz)>1.0
        swm(jz)=1.0;
    end
    if swt(jz)<0.0
        swt(jz)=0.001;
    end
    if swt(jz)>1.0
        swt(jz)=1.0;
    end
    if phiiw(jz)<0.0
        phiiw(jz)=0.001;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Resistivity Modeling 
lamdai=1.75;
lamdam=0.75;
tmpc=60;%%%%%% tempreture is known 
Cx=120;
wetcase=1;
[resist,ffac,rw,logr]=sdem_resist_sandstone(phit,phii,phim,swi,swm,lamdai,lamdam,tmpc,Cx,wetcase);

%%%%%%% manipulation for shale 
pwn_add_sh=0.01;
[resist,loge]=edge_log_intepol_noise(resist,1.0,nz,nsh,ksh,pwn_add_sh);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Velocity& Density Modeling 
%%%%%%SDEM model parameters
phic=0.45;
Lik=0.20;
Limu=0.20;
Lmk=1.0;
Lmmu=1.0;
%%%%%%%%%%% fluids 
sot=1-swt;
sgt=0;
satype='uniform';
kw=2.5; %% Gpa
k_oil=0.75; %% Gpa
k_gas=0.1; %% Gpa
rhow=1.05;
rho_oil=0.75;
rho_gas=0.04;%%gas

%%%%%%%%% properties of matrix(grain/inclusions)or soild rho_ma k_ma and mu_ma 
rho_q=2.65;      k_q=37;           mu_q=44;
rho_fp=2.52;     k_fp=29;          mu_fp=15;    %%%%%%% Orthoclase (Alkali Feldespar) with den=2.52 DTCO=69; and GR=220)

rho_clay=2.52;   k_clay=25.3519;   mu_clay=16.2580; %%%% (Inllite with density 2.52 GR=275, Vp=4.32 and Vs=2.54 and phi_sh=0.20 set at the biggening of the code 



cma=1-phii-csh;  %%%% matrix volume concentration

for jz=1:nz
    if cma(jz)<0.0
        cma(jz)=0.0001;        
    end
end

fcomp=0.75; %%%% relative percentage of mineral constituent 

cqz=fcomp.*cma; %%% percentage of quartz in matrix (scalar value)
cfp=(1-fcomp).*cma;%%%% percentage of feldespar in matrix (scalar value)

rho_ma=(cqz.*rho_q+cfp.*rho_fp)./cma;  %%%%% matrix or soild density(non-clay minerals)
rho_ma=mean(rho_ma(nsh:nz-nsh));

% % % % % cdc=1-phi_sh; %%%% volume concentration of dry clay in shale
% % % % % rhow_sh=rhow;   %%%% assume density of bound water is the same as non-shale water
% % % % % rho_sh=rhow_sh*phi_sh+cdc*rho_clay;

% % % %%%%% moduli upper bound for non-clay matrix 
k_voit_ma=(cqz.*k_q+cfp.*k_fp)./cma;
mu_voit_ma=(cqz.*mu_q+cfp.*mu_fp)./cma;

k_voit_ma=mean(k_voit_ma(nsh:nz-nsh));
mu_voit_ma=mean(mu_voit_ma(nsh:nz-nsh));
% % % %%%%% moduli lower bound for non-clay matrix 
k_ruess_ma=((cqz./k_q+cfp./k_fp)./cma).^-1;
mu_ruess_ma=((cqz./mu_q+cfp./mu_fp)./cma).^-1;

k_ruess_ma=mean(k_ruess_ma(nsh:nz-nsh));
mu_ruess_ma=mean(mu_ruess_ma(nsh:nz-nsh));

%%%%% moduli hill average for non-clay matrix
k_ma=0.5.*(k_ruess_ma+k_voit_ma); %%%% matrix bulk modulus 
mu_ma=0.5.*(mu_ruess_ma+mu_voit_ma); %%%% matrix shear modulus 

%%%%%%%%%%%%%%%%%matrix 
[dtco,dtsm,vp,vs,ai,si,vr,pr,ksat,musat,den,logk,logmu]=...
         sdem_sonic_sandstone_v01...
         (phi_sh,phit,phii,phim,phic,Lik,Limu,Lmk,Lmmu,...
         swi,swt,sot,sgt,satype,kw,k_oil,k_gas,rhow,rho_oil,rho_gas,csh,rho_clay,k_clay,mu_clay,rho_ma,k_ma,mu_ma); 
     
     
%%%%%%% manipulation for shale
%%%compute shale dtco and dtsm using illite properties (Vp=4.32 & Vs=2.54)
vp_sh=4.32;
vs_sh=2.54;
dtco_sh=1000/vp_sh; %%%msec/m
dtsm_sh=1000/vs_sh;
ai_sh=rho_clay.*vp_sh;
si_sh=rho_clay.*vs_sh;
vr_sh=vp_sh/vs_sh;
r=vr_sh;
pr_sh=(0.5*r^2-1)/(r^2-1);

[den,loge]=edge_log_intepol_noise(den,rho_clay,nz,nsh,ksh,pwn_add_sh);
[dtco,loge]=edge_log_intepol_noise(dtco,dtco_sh,nz,nsh,ksh,pwn_add_sh);
[dtsm,loge]=edge_log_intepol_noise(dtsm,dtsm_sh,nz,nsh,ksh,pwn_add_sh);
[vp,loge]=edge_log_intepol_noise(vp,vp_sh,nz,nsh,ksh,pwn_add_sh);
[vs,loge]=edge_log_intepol_noise(vs,vs_sh,nz,nsh,ksh,pwn_add_sh);
[ai,loge]=edge_log_intepol_noise(ai,ai_sh,nz,nsh,ksh,pwn_add_sh);
[si,loge]=edge_log_intepol_noise(si,si_sh,nz,nsh,ksh,pwn_add_sh);
[vr,loge]=edge_log_intepol_noise(vr,vr_sh,nz,nsh,ksh,pwn_add_sh);
[pr,loge]=edge_log_intepol_noise(pr,pr_sh,nz,nsh,ksh,pwn_add_sh);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GR modeling 
%%%% assume Illite with density 2.52 GR=275, Vp=4.32 and Vs=2.54 and phi_sh=0.20 set at the biggening of the code 
GR_clay=275;
%%%%%%%%%%%%%%% non-clay matrix 
GR_q=0.0;%%%% GR of Quartz
GR_fp=220; %%%%% Orthoclase (Alkali Feldespar) with den=2.52 DTCO=69; and GR=220)

GRma=(cqz.*rho_q.*GR_q+cfp.*rho_fp.*GR_fp)./(rho_ma.*cma);

[GR,logr]=GR_Sand_v01(phi_sh,rhow,den,GR_clay,csh,rho_clay,rho_ma,GRma,phii);
[GR,loge]=edge_log_intepol_noise(GR,GR_clay,nz,nsh,ksh,pwn_add_sh);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PLOTTING SETUP
depth_lim_all=[0 zlas];
fst=16; %%% title font size
fsx=16; %%% xlabel font size
fsy=16; %%% ylabel font size
fset=16; %%% set font size
lw=2;
fs=12;
phi_xlim=[-0.05 max(phit)+0.02];
sw_xlim=[-0.1 1.1];
resist_xlim=[0.1 100];
den_xlim=[min(den)-0.02 max(den)+0.02];
dtp_xlim=[min(dtco)-10 max(dtco)+10];
dts_xlim=[min(dtsm)-10 max(dtsm)+10];
GR_xlim=[min(GR)-20 max(GR)+20];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% final plotting of model 
figure(3)
subplot(1,4,1)
plot(phiiw,DEPTH,'b','LineWidth',lw)
hold on
plot(phii,DEPTH,'r','LineWidth',lw)
legend('Water-filled','Total')
xlabel('(v/v)','FontSize',fsx)
ylabel('Relative Depth(m)','FontSize',fsy)
title('Intg. porosity','FontSize',fst)
xlim(phi_xlim);
ylim(depth_lim_all);
grid MINOR
grid on
set(gca,'YDir','reverse','FontSize',fset)       
subplot(1,4,2)
plot(phim,DEPTH,'r','LineWidth',lw)
xlabel('(v/v)','FontSize',fsx)
title('Microporosity','FontSize',fst)
xlim(phi_xlim);
ylim(depth_lim_all);
grid MINOR
grid on
set(gca,'YDir','reverse','FontSize',fset)       
subplot(1,4,3)
plot(phiiw,DEPTH,'b',phit,DEPTH,'r','LineWidth',lw)
hold on
xlabel('(v/v)','FontSize',fsx)
title('Total porosity','FontSize',fst)
legend('Water-filled Intg.','Total')
xlim(phi_xlim);
ylim(depth_lim_all);
grid MINOR
grid on
set(gca,'YDir','reverse','FontSize',fset)       
subplot(1,4,4)
plot(swt,DEPTH,'k','LineWidth',lw)
hold on
plot(swi,DEPTH,'r','LineWidth',lw)
hold on
plot(swm,DEPTH,'b','LineWidth',lw)
legend('Total','Intg.','Bound')
xlabel('(v/v)','FontSize',fsx)
title('Water saturation','FontSize',fst)
xlim(sw_xlim);
ylim(depth_lim_all);
grid MINOR
grid on
set(gca,'YDir','reverse','FontSize',fset)
%%%%%%%%%%%%%%%%%%%%%%%QC 
nhist=40;
figure(4)
subplot(2,2,1)
hist(phii,nhist)
ylabel('Frequency','FontSize',fst)
title('Intergranular porosity(v/v)','FontSize',fst)
grid MINOR
grid on
set(gca,'FontSize',fset)
subplot(2,2,2)
hist(phim,nhist)
ylabel('Frequency','FontSize',fst)
title('Micro-porosity(v/v)','FontSize',fst)
grid MINOR
grid on
set(gca,'FontSize',fset)
subplot(2,2,3)
hist(swi,nhist)
ylabel('Frequency','FontSize',fst)
title('Intergranular water saturation(v/v)','FontSize',fst)
grid MINOR
grid on
set(gca,'FontSize',fset)
subplot(2,2,4)
hist(swm,nhist)
ylabel('Frequency','FontSize',fst)
title('Bound water saturation(v/v)','FontSize',fst)
grid MINOR
grid on
set(gca,'FontSize',fset)
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Data plotting 
figure(7)
subplot(1,8,1)
plot(phiiw,DEPTH,'b','LineWidth',lw)
hold on
plot(phii,DEPTH,'r','LineWidth',lw)
legend('Water-filled','Total')
xlabel('(v/v)','FontSize',fsx)
ylabel('Relative Depth(m)','FontSize',fsy)
title('Intg. porosity','FontSize',fst)
xlim(phi_xlim);
ylim(depth_lim_all);
grid MINOR
grid on
set(gca,'YDir','reverse','FontSize',fset)       
subplot(1,8,2)
plot(phim,DEPTH,'r','LineWidth',lw)
xlabel('(v/v)','FontSize',fsx)
title('Microporosity','FontSize',fst)
xlim(phi_xlim);
ylim(depth_lim_all);
grid MINOR
grid on
set(gca,'YDir','reverse','FontSize',fset)       
subplot(1,8,3)
plot(phiiw,DEPTH,'b',phit,DEPTH,'r','LineWidth',lw)
hold on
xlabel('(v/v)','FontSize',fsx)
title('Porosity','FontSize',fst)
legend('Water-filled Intg.','Total')
xlim(phi_xlim);
ylim(depth_lim_all);
grid MINOR
grid on
set(gca,'YDir','reverse','FontSize',fset)  
subplot(1,8,4)
plot(GR,DEPTH,'b','LineWidth',lw)
xlabel('gAPI','FontSize',fsx)
title('GR','FontSize',fst)
xlim(GR_xlim);
ylim(depth_lim_all);
grid MINOR
grid on
set(gca,'YDir','reverse','FontSize',fset)
subplot(1,8,5)
semilogx(resist,DEPTH,'b','LineWidth',lw)
xlabel('ohm-m','FontSize',fsx)
title('Resistivity','FontSize',fst)
xlim(resist_xlim);
ylim(depth_lim_all);
grid MINOR
grid on
set(gca,'YDir','reverse','FontSize',fset)       
hold on
subplot(1,8,6)
plot(den,DEPTH,'r','LineWidth',lw)
xlabel('gr/cm3','FontSize',fsx)
title('Density','FontSize',fst)
xlim(den_xlim);
ylim(depth_lim_all);
grid MINOR
grid on
set(gca,'YDir','reverse','FontSize',fset)      
subplot(1,8,7)
plot(dtco,DEPTH,'k','LineWidth',lw)
xlabel('microsec/m','FontSize',fsx)
title('DTCO','FontSize',fst)
xlim(dtp_xlim);
ylim(depth_lim_all);
grid MINOR
grid on
set(gca,'YDir','reverse','FontSize',fset)       
hold on
subplot(1,8,8)
plot(dtsm,DEPTH,'m','LineWidth',lw)
xlabel('microsec/m','FontSize',fsx)
title('DTSM','FontSize',fst)
xlim(dts_xlim);
ylim(depth_lim_all);
grid MINOR
grid on
set(gca,'YDir','reverse','FontSize',fset)  

%%%%%%%%%%%%%%%%%%%%%%% QC SDEM 
kks=1;
kkf=nz; %%% oil zone kk=1 to kk=72 and brine 100 to 156=np

mksz=5; %%%%% MarkerSize in x-plot

figure(8)
subplot(2,3,1)
plot(phit(kks:kkf),dtco(kks:kkf),'bo','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',mksz)
hold on
plot(phit(kks:kkf),dtsm(kks:kkf),'ro','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',mksz)
xlabel('Total Porosity (v/v)','FontSize',fs)
ylabel('Sonic Travel Time (microsec/m)','FontSize',fs)
legend('DTCO','DTSM')
xlim(phi_xlim);
grid on
set(gca,'FontSize',fset)      
subplot(2,3,2)
plot(phii(kks:kkf),dtco(kks:kkf),'bo','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',mksz)
hold on
plot(phii(kks:kkf),dtsm(kks:kkf),'ro','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',mksz)
xlabel('Intergranular Porosity (v/v)','FontSize',fs)
ylabel('Sonic Travel Time (microsec/m)','FontSize',fs)
legend('DTCO','DTSM')
xlim(phi_xlim);
grid on
set(gca,'FontSize',fset)      
subplot(2,3,3)
plot(phim(kks:kkf),dtco(kks:kkf),'bo','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',mksz)
hold on
plot(phim(kks:kkf),dtsm(kks:kkf),'ro','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',mksz)
xlabel('Micro-porosity (v/v)','FontSize',fs)
ylabel('Sonic Travel Time (microsec/m)','FontSize',fs)
legend('DTCO','DTSM')
xlim(phi_xlim);
grid on
set(gca,'FontSize',fset)      
subplot(2,3,6)
plot(ai(kks:kkf),vr(kks:kkf),'bo','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',mksz)
xlabel('Acoustic Impedance[Vp*Den](Km/sec*gr/cm3)','FontSize',fs)
ylabel('Velocity Ratio (Vp/Vs)','FontSize',fs)
grid on
set(gca,'FontSize',fset)      
subplot(2,3,4)
plot(phit(kks:kkf),vp(kks:kkf),'bo','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',mksz)
hold on
plot(phit,vs,'ro','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',mksz)
xlabel('Total Porosity (v/v)','FontSize',fs)
ylabel('Velocity (Km/sec)','FontSize',fs)
legend('P-wave','S-wave')
xlim(phi_xlim);
grid on
set(gca,'FontSize',fset) 
subplot(2,3,5)
plot(phii(kks:kkf),vp(kks:kkf),'bo','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',mksz)
hold on
plot(phii,vs,'ro','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',mksz)
xlabel('Intergranular Porosity (v/v)','FontSize',fs)
ylabel('Velocity (Km/sec)','FontSize',fs)
legend('P-wave','S-wave')
xlim(phi_xlim);
grid on
set(gca,'FontSize',fset)
subplot(2,3,6)
plot(phim(kks:kkf),vp(kks:kkf),'bo','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',mksz)
hold on
plot(phim(kks:kkf),vs(kks:kkf),'ro','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',mksz)
xlabel('Micro-porosity (v/v)','FontSize',fs)
ylabel('Velocity (Km/sec)','FontSize',fs)
legend('P-wave','S-wave')
xlim(phi_xlim);
grid on
set(gca,'FontSize',fset)

%%%%%%%%%%%%%%%%%%%%%%% QC RPT 

mksz=70; %%%%% MarkerSize in x-plot
figure(9)
subplot(2,2,1)
scatter(ai(kks:kkf),vr(kks:kkf),mksz,phit(kks:kkf),'filled')
xlabel('Acoustic Impedance[Vp*Den](Km/sec*gr/cm3)','FontSize',fs)
ylabel('Velocity Ratio [Vp/Vs]','FontSize',fs)
grid on
cmin=0.0; cmax=max(phit(kks:kkf));
c = colorbar; c.Label.String = 'Total Porosity (v/v)';
colormap(jet); caxis([cmin cmax])
set(gca,'FontSize',fset)
subplot(2,2,2)
scatter(ai(kks:kkf),vr(kks:kkf),mksz,phii(kks:kkf),'filled')
xlabel('Acoustic Impedance[Vp*Den](Km/sec*gr/cm3)','FontSize',fs)
ylabel('Velocity Ratio [Vp/Vs]','FontSize',fs)
grid on
cmin=0.0; cmax=max(phii(kks:kkf));
c = colorbar; c.Label.String = 'Intergranular Porosity (v/v)';
colormap(jet); caxis([cmin cmax])
set(gca,'FontSize',fset)
subplot(2,2,3)
scatter(ai(kks:kkf),vr(kks:kkf),mksz,phim(kks:kkf),'filled')
xlabel('Acoustic Impedance[Vp*Den](Km/sec*gr/cm3)','FontSize',fs)
ylabel('Velocity Ratio [Vp/Vs]','FontSize',fs)
grid on
cmin=0.0; cmax=max(phim(kks:kkf));
c = colorbar; c.Label.String = 'Micro-porosity (v/v)';
colormap(jet); caxis([cmin cmax])
set(gca,'FontSize',fset)
subplot(2,2,4)
scatter(ai(kks:kkf),vr(kks:kkf),mksz,swt(kks:kkf),'filled')
xlabel('Acoustic Impedance[Vp*Den](Km/sec*gr/cm3)','FontSize',fs)
ylabel('Velocity Ratio [Vp/Vs]','FontSize',fs)
grid on
cmin=0.0; cmax=max(swt(kks:kkf));
c = colorbar; c.Label.String = 'Water Saturation (v/v)';
colormap(jet); caxis([cmin cmax])
set(gca,'FontSize',fset)


mksz=70; %%%%% MarkerSize in x-plot
figure(10)
subplot(2,2,1)
scatter(dtco(kks:kkf),log10(resist(kks:kkf)),mksz,swt(kks:kkf),'filled')
xlabel('DTCO(microsec/m)','FontSize',fs)
ylabel('Resistivity [log10](ohm-m)','FontSize',fs)
grid on
cmin=0.0; cmax=max(swt(kks:kkf));
c = colorbar; c.Label.String = 'Water Saturation (v/v)';
colormap(jet); caxis([cmin cmax])
set(gca,'FontSize',fset)
subplot(2,2,2)
scatter(dtco(kks:kkf),log10(resist(kks:kkf)),mksz,phim(kks:kkf),'filled')
xlabel('DTCO(microsec/m)','FontSize',fs)
ylabel('Resistivity [log10](ohm-m)','FontSize',fs)
grid on
cmin=0.0; cmax=max(phim(kks:kkf));
c = colorbar; c.Label.String = 'Micro-porosity (v/v)';
colormap(jet); caxis([cmin cmax])
set(gca,'FontSize',fset)
subplot(2,2,3)
scatter(dtco(kks:kkf),log10(resist(kks:kkf)),mksz,phii(kks:kkf),'filled')
xlabel('DTCO(microsec/m)','FontSize',fs)
ylabel('Resistivity [log10](ohm-m)','FontSize',fs)
grid on
cmin=0.0; cmax=max(phii(kks:kkf));
c = colorbar; c.Label.String = 'Intergranular Porosity (v/v)';
colormap(jet); caxis([cmin cmax])
set(gca,'FontSize',fset)
subplot(2,2,4)
scatter(dtsm(kks:kkf),log10(resist(kks:kkf)),mksz,phim(kks:kkf),'filled')
xlabel('DTSM(microsec/m)','FontSize',fs)
ylabel('Resistivity [log10](ohm-m)','FontSize',fs)
grid on
cmin=0.0; cmax=max(phim(kks:kkf));
c = colorbar; c.Label.String = 'Micro-porosity (v/v)';
colormap(jet); caxis([cmin cmax])
set(gca,'FontSize',fset)

%%%%%%%%%%%%%%%%%%%%%%% QC SDEM 
mksz=70; %%%%% MarkerSize in x-plot
figure(11)
kks=1;
kkf=nz; %%% oil zone kk=1 to kk=72 and brine 100 to 156=np
subplot(2,3,1)
scatter(ai(kks:kkf),vr(kks:kkf),mksz,swt(kks:kkf),'filled')
xlabel('Acoustic Impedance(Km/sec*gr/cm3)','FontSize',fs)
ylabel('Velocity Ratio [Vp/Vs]','FontSize',fs)
grid on
cmin=0.0; cmax=max(swt(kks:kkf));
c = colorbar; c.Label.String = 'Water Saturation (v/v)';
colormap(jet); caxis([cmin cmax])
set(gca,'FontSize',fset)
subplot(2,3,2)
scatter(ai(kks:kkf),vr(kks:kkf),mksz,phii(kks:kkf),'filled')
xlabel('Acoustic Impedance(Km/sec*gr/cm3)','FontSize',fs)
% % ylabel('Velocity Ratio [Vp/Vs]','FontSize',fs)
grid on
cmin=0.0; cmax=max(phii(kks:kkf));
c = colorbar; c.Label.String = 'Intergranular Porosity (v/v)';
colormap(jet); caxis([cmin cmax])
set(gca,'FontSize',fset)
subplot(2,3,4)
scatter(dtco(kks:kkf),log10(resist(kks:kkf)),mksz,swt(kks:kkf),'filled')
xlabel('DTCO(microsec/m)','FontSize',fs)
ylabel('Resistivity [log10](ohm-m)','FontSize',fs)
grid on
cmin=0.0; cmax=max(swt(kks:kkf));
c = colorbar; c.Label.String = 'Water Saturation (v/v)';
colormap(jet); caxis([cmin cmax])
set(gca,'FontSize',fset)
subplot(2,3,5)
scatter(dtco(kks:kkf),log10(resist(kks:kkf)),mksz,phii(kks:kkf),'filled')
xlabel('DTCO(microsec/m)','FontSize',fs)
% % ylabel('Resistivity [log10](ohm-m)','FontSize',fs)
grid on
cmin=0.0; cmax=max(phii(kks:kkf));
c = colorbar; c.Label.String = 'Intergranular Porosity (v/v)';
colormap(jet); caxis([cmin cmax])
set(gca,'FontSize',fset)
subplot(2,3,3)
scatter(ai(kks:kkf),vr(kks:kkf),mksz,GR(kks:kkf),'filled')
xlabel('Acoustic Impedance(Km/sec*gr/cm3)','FontSize',fs)
% % ylabel('Velocity Ratio [Vp/Vs]','FontSize',fs)
grid on
cmin=0.0; cmax=max(GR(kks:kkf));
c = colorbar; c.Label.String = 'GR (gAPI)';
colormap(jet); caxis([cmin cmax])
set(gca,'FontSize',fset)
subplot(2,3,6)
scatter(dtco(kks:kkf),log10(resist(kks:kkf)),mksz,GR(kks:kkf),'filled')
xlabel('DTCO(microsec/m)','FontSize',fs)
% % % ylabel('Resistivity [log10](ohm-m)','FontSize',fs)
grid on
cmin=0.0; cmax=max(GR(kks:kkf));
c = colorbar; c.Label.String = 'GR (gAPI)';
colormap(jet); caxis([cmin cmax])
set(gca,'FontSize',fset)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mksz=70; %%%%% MarkerSize in x-plot
kks=1;
kkf=nz; %%% oil zone kk=1 to kk=72 and brine 100 to 156=npfigure(111)
figure(111)
subplot(2,2,1)
scatter(log10(resist(kks:kkf)),vr(kks:kkf),mksz,GR(kks:kkf),'filled')
xlabel('Resistivity [log10](ohm-m)','FontSize',fs)
ylabel('Velocity Ratio [Vp/Vs]','FontSize',fs)
grid on
cmin=0.0; cmax=max(GR(kks:kkf));
c = colorbar; c.Label.String = 'GR (gAPI)';
colormap(jet); caxis([cmin cmax])
set(gca,'FontSize',fset)
subplot(2,2,2)
scatter(log10(resist(kks:kkf)),vr(kks:kkf),mksz,swt(kks:kkf),'filled')
xlabel('Resistivity [log10](ohm-m)','FontSize',fs)
ylabel('Velocity Ratio [Vp/Vs]','FontSize',fs)
grid on
cmin=0.0; cmax=max(swt(kks:kkf));
c = colorbar; c.Label.String = 'Total water saturation (v/v)';
colormap(jet); caxis([cmin cmax])
set(gca,'FontSize',fset)
subplot(2,2,3)
scatter(log10(resist(kks:kkf)),vr(kks:kkf),mksz,phii(kks:kkf),'filled')
xlabel('Resistivity [log10](ohm-m)','FontSize',fs)
ylabel('Velocity Ratio [Vp/Vs]','FontSize',fs)
grid on
cmin=0.0; cmax=max(phii(kks:kkf));
c = colorbar; c.Label.String = 'Intergranular porosity (v/v)';
colormap(jet); caxis([cmin cmax])
set(gca,'FontSize',fset)
subplot(2,2,4)
scatter(log10(resist(kks:kkf)),vr(kks:kkf),mksz,phim(kks:kkf),'filled')
xlabel('Resistivity [log10](ohm-m)','FontSize',fs)
ylabel('Velocity Ratio [Vp/Vs]','FontSize',fs)
grid on
cmin=0.0; cmax=max(phim(kks:kkf));
c = colorbar; c.Label.String = 'Micro-porosity (v/v)';
colormap(jet); caxis([cmin cmax])
set(gca,'FontSize',fset)
%%%%%%%%%%%%%%%%%%%%%%%QC 
nhist=40;
figure(12)
subplot(2,2,1)
hist(resist,nhist)
ylabel('Frequency','FontSize',fst)
title('Resistivity (ohm-m)','FontSize',fst)
grid MINOR
grid on
set(gca,'FontSize',fset)
subplot(2,2,2)
hist(den,nhist)
ylabel('Frequency','FontSize',fst)
title('Density (gr/cm3)','FontSize',fst)
grid MINOR
grid on
set(gca,'FontSize',fset)
subplot(2,2,3)
hist(dtco,nhist)
ylabel('Frequency','FontSize',fst)
title('DTCO(microsec/m)','FontSize',fst)
grid MINOR
grid on
set(gca,'FontSize',fset)
subplot(2,2,4)
hist(dtsm,nhist)
ylabel('Frequency','FontSize',fst)
title('DTSM (microsec/m)','FontSize',fst)
grid MINOR
grid on
set(gca,'FontSize',fset)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Model Parametrization
%%%%%%%% Final constant & depth-varying model parameters to create PPD 
%%%%%%%% after each VFSA run, the best model will be collected
%%%%%%%%% and will be stored in column jrun=1 to nruns
% % np=length(phit); %%%%% number of depth intervals 
phii_opt=zeros(nruns,nz); 
phit_opt=zeros(nruns,nz); 
swt_opt=zeros(nruns,nz);
lamdai_opt=zeros(1,nruns); %%%% matrix lithology exponent 
lamdam_opt=zeros(1,nruns); %%%% Micro  lithology exponent
Lik_opt=zeros(1,nruns); %%%% L parameter for intergranular porosity associated with bulk modulus
Limu_opt=zeros(1,nruns); %%%% L parameter for intergranular porosity associated with shear modulus
Cx_opt=zeros(1,nruns); %%%% salinity of flushed zone
phic_opt=zeros(1,nruns); %%%% critical porosity 
ks_opt=zeros(1,nruns); %%%% matrix bulk modulus 
mus_opt=zeros(1,nruns); %%%% matrix shear modulus  
rhos_opt=zeros(1,nruns); %%%% matrix density
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% VFSA secondary parameters
%%%% range of model parameters for optimization

lamdai_max=1.65;     lamdai_min=1.55;      dlamdai=0.005;       nlamdai=round(abs(lamdai_max-lamdai_min)./dlamdai);
Cx_max=140.0;       Cx_min=100.0;          dCx=1;              nCx=round(abs(Cx_max-Cx_min)./dCx);
Lik_max=0.20;        Lik_min=0.05;         dLik=0.01;          nLik=round(abs(Lik_max-Lik_min)./dLik);
Limu_max=0.20;       Limu_min=0.05;        dLimu=0.01;         nLimu=round(abs(Limu_max-Limu_min)./dLimu);
phic_max=0.55;       phic_min=0.35;        dphic=0.01;         nphic=round(abs(phic_max-phic_min)./dphic);
ks_max=40.00;         ks_min=25.00;        dks=0.1;            nks=round(abs(ks_max-ks_min)./dks);
mus_max=40.00;         mus_min=25.00;      dmus=0.1;           nmus=round(abs(mus_max-mus_min)./dmus);
rhos_max=2.75;         rhos_min=2.5;       drhos=0.01;        nrhos=round(abs(rhos_max-rhos_min)./drhos);

%%%%%%%%%%%%%% Minimum range of global parameters only for QC
% % % lamdai_max=1.65;     lamdai_min=1.55;     dlamdai=0.005;       nlamdai=round(abs(lamdai_max-lamdai_min)./dlamdai);
% % % Cx_max=122.0;       Cx_min=118.0;          dCx=1;              nCx=round(abs(Cx_max-Cx_min)./dCx);
% % % Lik_max=0.105;        Lik_min=0.095;        dLik=0.01;          nLik=round(abs(Lik_max-Lik_min)./dLik);
% % % Limu_max=0.085;       Limu_min=0.075;       dLimu=0.01;         nLimu=round(abs(Limu_max-Limu_min)./dLimu);
% % % phic_max=0.46;       phic_min=0.44;        dphic=0.01;         nphic=round(abs(phic_max-phic_min)./dphic);
% % % ks_max=78.00;         ks_min=76.00;        dks=0.1;            nks=round(abs(ks_max-ks_min)./dks);
% % % mus_max=33.00;         mus_min=31.00;      dmus=0.1;           nmus=round(abs(mus_max-mus_min)./dmus);
% % % rhos_max=2.72;         rhos_min=2.70;       drhos=0.01;        nrhos=round(abs(rhos_max-rhos_min)./drhos);



%%%%%%%%%%%%%%plotting 
Cx_lim=[Cx_min Cx_max];
lamdai_lim=[lamdai_min lamdai_max];
phic_lim=[phic_min phic_max];
Lik_lim=[Lik_min Lik_max];
Limu_lim=[Limu_min Limu_max];
rhos_lim=[rhos_min rhos_max];
ks_lim=[ks_min ks_max];
mus_lim=[mus_min mus_max];

%%%%%%%%%%% pre-defined matrices
phim_mod=zeros(nz,1);    phii_mod=zeros(nz,1);   phit_mod=zeros(nz,1);  
swm_mod=zeros(nz,1);     swi_mod=zeros(nz,1);    swt_mod=zeros(nz,1);
phim_trial=zeros(nz,1);  phii_trial=zeros(nz,1); phit_trial=zeros(nz,1);
swm_trial=zeros(nz,1);   swi_trial=zeros(nz,1);  swt_trial=zeros(nz,1);

npar=3; %%%%% number of continious variables to be estimated....
%%%%%%%%%%%%%%(here intergranular porosity, & total porosity and total water saturations)

%%%%% the following 4 matrices has three columns,...
%%%%% 1st column is intergranular porosity,...
%%%%% 2nd column is total porosity...
%%%%% 3rd column is total water saturation
emod_local=zeros(nz,npar); 
etrial_local=zeros(nz,npar); 
model_mod=zeros(nz,npar); 
model_trial=zeros(nz,npar);
model_min=zeros(nz,npar); 
model_max=zeros(nz,npar); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%% plotting of constant model parameters(constant value for the entire log interval) 
lamdai_plt=zeros(niter,nruns);
Cx_plt=zeros(niter,nruns);
Lik_plt=zeros(niter,nruns);
Limu_plt=zeros(niter,nruns);
phic_plt=zeros(niter,nruns);
rhos_plt=zeros(niter,nruns);
ks_plt=zeros(niter,nruns);
mus_plt=zeros(niter,nruns);


err_plt=zeros(niter,nruns); %%% tracer for error at different iteration per run
emod_local_plt=zeros(niter,nruns,nz,npar);


temp_min=0.00001;%%%common min temp for global and local 
decay_local=ones(nz,npar).*decay; %%%% Same decay for all parameters, 
tmpi_local=ones(nz,npar).*temp0; %%% initial temp for all layers and all model parameters (4 here)


%%%%%%%%%%%%%%%%%%%%%% known and useful bounds on model parameters 
%%%%%%%%%%% how to obtain maximum phit for sampling 
%%%%%%%%%%%% compute density porosity assumiung water-filled and known heavy matrix for phit upper bound 
rho_ma=rhos_max; %%%% come from sdem_sonic_sandstone.m
phitw_den=(den-rho_ma)./(rhow-rho_ma);

ubphit=0.01; %%%% upper bound 
phit_max=phitw_den+ubphit;%%%% upper bound of total porosity

%%%%%%%%%%%% compute density porosity assumiung oil-filled and known light matrix for phit lower bound 
rho_ma=rhos_min; %%%% come from sdem_sonic_carbonate.m
phitg_den=(den-rho_ma)./(rho_oil-rho_ma);

lbphit=0.01; %%%% upper bound 
phit_min=phitg_den-lbphit;%%%% lower bound of total porosity


phimax= max(phit_max); %%%% comes from figure 7000
phimax_line= phimax.*ones(size(DEPTH)); %%%% comes from figure 7000
phic_line=phic.*ones(size(DEPTH));

%%%%%%%%%%%
figure(13)
plot(phit,DEPTH,'k',phitw_den,DEPTH,'b',phitg_den,DEPTH,'r','LineWidth',lw)
hold on
plot(phic_line,DEPTH,'m','LineWidth',lw)
hold on
plot(phit_min,DEPTH,'--k',phit_max,DEPTH,':k','LineWidth',lw)
title('Porosity(v/v)','FontSize',fst)
ylabel('Relative Depth(m)','FontSize',fsy)
legend('True Total Porosity','Water-filled Density Porosity(Heavy Matrix)',...
    'Oil-filled Density Porosity(Light Matrix)','Critical Porosity','Priori Lower Bound','Priori Upper Bound')
xlim([0 phic+0.05])
ylim(depth_lim_all)
grid MINOR
grid on
set(gca,'YDir','reverse','FontSize',fset)       

%%%%%%%%%%%%%%%%%%%%%% other known and useful bounds on model parameters 
phii_min=0.0001; 
phii_max= phit_max; %%%% comes from figure 7000

dphi=0.001; 
nphii=round(abs(phii_max-phii_min)/dphi);

nphit=round(abs(phit_max-phit_min)/dphi);

sw_min=0.0; 
sw_max=1.0;
dsw=0.001; 
nsw=round(abs(sw_max-sw_min)/dsw);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Main VFSA loop
for jrun=1:nruns    
         
% % %         lamdai_mod=monte_carlo_sample_500try(lamdai_min,dlamdai,nlamdai,lamdai_max);
% % %         Cx_mod=monte_carlo_sample_500try(Cx_min,dCx,nCx,Cx_max);
% % %         Lik_mod=monte_carlo_sample_500try(Lik_min,dLik,nLik,Lik_max);
% % %         Limu_mod=monte_carlo_sample_500try(Limu_min,dLimu,nLimu,Limu_max);
% % %         phic_mod=monte_carlo_sample_500try(phic_min,dphic,nphic,phic_max);    
% % %         rhos_mod=monte_carlo_sample_500try(rhos_min,drhos,nrhos,rhos_max); 
% % %         ks_mod=monte_carlo_sample_500try(ks_min,dks,nks,ks_max); 
% % %         mus_mod=monte_carlo_sample_500try(mus_min,dmus,nmus,mus_max);
        
% % % % % %         %%%%% rough initial guess for global parameters only for better
% % % % % %         %%%%% plotting of SEG 2019
        lamdai_mod=lamdai_min;
        Cx_mod=Cx_min;
        Lik_mod=Lik_max;
        Limu_mod=Limu_min;
        phic_mod=phic_max;    
        rhos_mod=rhos_min; 
        ks_mod=ks_max; 
        mus_mod=mus_max; 

                       
        for jz=1:nz
            
            logm=-1;
            while logm==-1
                phit_mod(jz)=monte_carlo_sample_500try(phit_min(jz),dphi,nphit(jz),phit_max(jz));%%% Initial guess for total porosity
                phii_mod(jz)=monte_carlo_sample_500try(phii_min,dphi,nphii(jz),phii_max(jz));%%% Initial guess  
                phim_mod(jz)=phit_mod(jz)-phii_mod(jz);
                if phim_mod(jz)< 0.0
                    logm=-1;
                else
                    logm=1;
                end
            end 
            
            swt_mod(jz)=monte_carlo_sample_500try(sw_min,dsw,nsw,sw_max);%%% Initial guess
        end
        swm_mod=swt_mod;
        swi_mod=swt_mod;
        sot_mod=1-swt_mod;
        sgt_mod=1-swt_mod-sot_mod;    
        
        
        
        %%%%%% 1st column is intergranular porosity
        %%%%%% 2nd column is total porosity
        %%%%%% 3rd column is total saturation which is equal to intergranular and vuyyg saturation 
        model_min(:,1)=phii_min; model_min(:,2)=phit_min; model_min(:,3)=sw_min;  
        model_max(:,1)=phii_max; model_max(:,2)=phit_max; model_max(:,3)=sw_max;
        
              
               
        %%%%%% compute gloabl error
        [emod,logr,logk,logmu]=objfun_sdem_density_sonic_resist_carbonate_v01(resist,dtco,dtsm,den,phit_mod,...
                               phic_mod,phii_mod,phim_mod,swi_mod,swm_mod,lamdai_mod,lamdam,tmpc,Cx_mod,wetcase,...
                               Lik_mod,Limu_mod,Lmk,Lmmu,swt_mod,sot_mod,sgt_mod,satype,kw,k_oil,k_gas,rhow,rho_oil,rho_gas,ks_mod,mus_mod,rhos_mod);
        
        %%%%%% compute local error  and Initialize local temperatures 
        if wloc~=0
            for jz=1:nz
                
                 [emod_local(jz,1:npar),logr,logk,logmu]=objfun_sdem_density_sonic_resist_carbonate_v01(resist(jz),dtco(jz),dtsm(jz),den(jz),phit_mod(jz),...
                 phic_mod,phii_mod(jz),phim_mod(jz),swi_mod(jz),swm_mod(jz),lamdai_mod,lamdam,tmpc,Cx_mod,wetcase,...
                               Lik_mod,Limu_mod,Lmk,Lmmu,swt_mod(jz),sot_mod(jz),sgt_mod(jz),satype,kw,k_oil,k_gas,rhow,rho_oil,rho_gas,ks_mod,mus_mod,rhos_mod);                
                
            end
            emod_local=wloc.*emod_local+wglob.*emod;
            model_mod=[phii_mod,phit_mod,swt_mod]; 
            ktemp=1;
            tmp_local=tmpi_local.*exp(-decay_local.*(ktemp-1).^0.5); %% local temp(matrix[nlayers,nparam=3])
            local_temp(:,:,ktemp)=tmp_local;
            jjtemp=1; %%%% for Reannealing
        end
        
        %%%%%%%%%%%%%%%%%%% start global tempreture loop
        
        jtemp=1;
        
        lamdai_plt(jtemp,jrun)=lamdai_mod;
        Cx_plt(jtemp,jrun)=Cx_mod;
        Lik_plt(jtemp,jrun)=Lik_mod;
        Limu_plt(jtemp,jrun)=Limu_mod; 
        phic_plt(jtemp,jrun)=phic_mod; 
        rhos_plt(jtemp,jrun)=rhos_mod; 
        ks_plt(jtemp,jrun)=ks_mod; 
        mus_plt(jtemp,jrun)=mus_mod; 
        err_plt(jtemp,jrun)=emod;     
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% VFSA tempreture loop
        while(jtemp<=niter-1)
            temp(jtemp)=temp0.*exp(-decay.*(jtemp-1).^0.5);
            tmp=temp0.*exp(-decay.*(jtemp-1).^0.5); 
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%% VFSA move loop
            for jmov=1:nmov
                
                %%%%%%%%%%% walk with global tempreture
                
                lamdai_trial=walk(lamdai_mod,lamdai_min,lamdai_max,tmp); lamdai_trial=round((lamdai_trial-lamdai_min)/dlamdai)*dlamdai+lamdai_min;
                Cx_trial=walk(Cx_mod,Cx_min,Cx_max,tmp); Cx_trial=round((Cx_trial-Cx_min)/dCx)*dCx+Cx_min;
                phic_trial=walk(phic_mod,phic_min,phic_max,tmp); phic_trial=round((phic_trial-phic_min)/dphic)*dphic+phic_min;
                Lik_trial=walk(Lik_mod,Lik_min,Lik_max,tmp); Lik_trial=round((Lik_trial-Lik_min)/dLik)*dLik+Lik_min;
                Limu_trial=walk(Limu_mod,Limu_min,Limu_max,tmp); Limu_trial=round((Limu_trial-Limu_min)/dLimu)*dLimu+Limu_min;    
                rhos_trial=walk(rhos_mod,rhos_min,rhos_max,tmp); rhos_trial=round((rhos_trial-rhos_min)/drhos)*drhos+rhos_min;
                ks_trial=walk(ks_mod,ks_min,ks_max,tmp); ks_trial=round((ks_trial-ks_min)/dks)*dks+ks_min;
                mus_trial=walk(mus_mod,mus_min,mus_max,tmp); mus_trial=round((mus_trial-mus_min)/dmus)*dmus+mus_min;
                     
                                             
                for jz=1:nz
                       logw=-1;
                       while logw==-1
                           phit_trial(jz)=walk(phit_mod(jz),phit_min(jz),phit_max(jz),tmp); phit_trial(jz)=round(( phit_trial(jz)-phit_min(jz))/dphi)*dphi+ phit_min(jz);
                           phii_trial(jz)=walk(phii_mod(jz),phii_min,phii_max(jz),tmp); phii_trial(jz)=round(( phii_trial(jz)-phii_min)/dphi)*dphi+ phii_min;
                           phim_trial(jz)=phit_trial(jz)-phii_trial(jz);
                           if phim_trial(jz) < 0.0
                               logw=-1;
                           else
                               logw=1;
                           end
                       end                           
                       
                       swt_trial(jz)=walk(swt_mod(jz),sw_min,sw_max,tmp); swt_trial(jz)=round(( swt_trial(jz)-sw_min)/dsw)*dsw+ sw_min;                       
                end        

                swi_trial=swt_trial;
                swm_trial=swt_trial;
                sot_trial=1-swt_trial;
                sgt_trial=1-swt_trial-sot_trial;
                
                

                [etrial,logr,logk,logmu]=objfun_sdem_density_sonic_resist_carbonate_v01(resist,dtco,dtsm,den,phit_trial,...
                               phic_trial,phii_trial,phim_trial,swi_trial,swm_trial,lamdai_trial,lamdam,tmpc,Cx_trial,wetcase,...
                               Lik_trial,Limu_trial,Lmk,Lmmu,swt_trial,sot_trial,sgt_trial,satype,kw,k_oil,k_gas,rhow,rho_oil,rho_gas,ks_trial,mus_trial,rhos_trial);

                                 
                 %%%%%%%%%%% walk with local tempreture
                  if (wloc~=0 && jtemp >= iter_loc && mod(jtemp,freq_loc)==0)

% % %                  if (wloc~=0 && jtemp >= iter_loc)
                     for jz=1:nz
                        logw=-1;
                        while logw==-1 
                            phii_trial(jz)=walk(phii_mod(jz),phii_min,phii_max(jz),tmp_local(jz,1));  phii_trial(jz)=round(( phii_trial(jz)-phii_min)/dphi)*dphi+ phii_min;
                            phit_trial(jz)=walk(phit_mod(jz),phit_min(jz),phit_max(jz),tmp_local(jz,2));  phit_trial(jz)=round(( phit_trial(jz)-phit_min(jz))/dphi)*dphi+ phit_min(jz);
                            phim_trial(jz)=phit_trial(jz)-phii_trial(jz);
                            if phim_trial(jz) < 0.0
                               logw=-1;
                            else
                               logw=1;
                            end
                        end           
                        
                        swt_trial(jz)=walk(swt_mod(jz),sw_min,sw_max,tmp_local(jz,3));  swt_trial(jz)=round(( swt_trial(jz)-sw_min)/dsw)*dsw+ sw_min;
                     
                        swi_trial(jz)=swt_trial(jz);
                        swm_trial(jz)=swt_trial(jz);
                        sot_trial(jz)=1-swt_trial(jz);
                        sgt_trial(jz)=1-swt_trial(jz)-sot_trial(jz);
                        
                        [etrial_local(jz,1:npar),logr,logk,logmu]=objfun_sdem_density_sonic_resist_carbonate_v01(resist(jz),dtco(jz),dtsm(jz),den(jz),phit_trial(jz),...
                         phic_trial,phii_trial(jz),phim_trial(jz),swi_trial(jz),swm_trial(jz),lamdai_trial,lamdam,tmpc,Cx_trial,wetcase,...
                         Lik_trial,Limu_trial,Lmk,Lmmu,swt_trial(jz),sot_trial(jz),sgt_trial(jz),satype,kw,k_oil,k_gas,rhow,rho_oil,rho_gas,ks_trial,mus_trial,rhos_trial);   
               
                     end
                     %%%%% recompute global error using models locally updated
                     [etrial,logr,logk,logmu]=objfun_sdem_density_sonic_resist_carbonate_v01(resist,dtco,dtsm,den,phit_trial,...
                               phic_trial,phii_trial,phim_trial,swi_trial,swm_trial,lamdai_trial,lamdam,tmpc,Cx_trial,wetcase,...
                               Lik_trial,Limu_trial,Lmk,Lmmu,swt_trial,sot_trial,sgt_trial,satype,kw,k_oil,k_gas,rhow,rho_oil,rho_gas,ks_trial,mus_trial,rhos_trial);                                   
                     
                     model_trial=[phii_trial,phit_trial,swt_trial]; 
                     etrial_local=wloc.*etrial_local+wglob.*etrial;
                 end                    
                             
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Global error update                
                if etrial< emod
                    %%%% hist_updat
                    emod=etrial;
                    lamdai_mod=lamdai_trial;
                    Cx_mod=Cx_trial;
                    phic_mod=phic_trial;
                    rhos_mod=rhos_trial;
                    ks_mod=ks_trial;
                    mus_mod=mus_trial;
                    Lik_mod=Lik_trial;
                    Limu_mod=Limu_trial;                    
                    for jz=1:nz
                        phii_mod(jz)=phii_trial(jz);
                        phit_mod(jz)=phit_trial(jz);
                        phim_mod(jz)=phit_trial(jz)-phii_trial(jz);
                        swt_mod(jz)=swt_trial(jz);
                        swi_mod(jz)=swt_trial(jz);
                        swm_mod(jz)=swt_trial(jz);                        
                    end                    
                else
                    arg=(etrial-emod)./temp(jtemp);
                    if arg>1.e6
                        pde=0.001;
                    else
                        pde=exp(-arg);
                    end
                    if pde>rand
                        %%%% hist_updat
                        emod=etrial;
                        lamdai_mod=lamdai_trial;
                        Cx_mod=Cx_trial;
                        phic_mod=phic_trial;
                        rhos_mod=rhos_trial;
                        ks_mod=ks_trial;
                        mus_mod=mus_trial;
                        Lik_mod=Lik_trial;
                        Limu_mod=Limu_trial; 
                        for jz=1:nz
                            phii_mod(jz)=phii_trial(jz);
                            phit_mod(jz)=phit_trial(jz);
                            phim_mod(jz)=phit_trial(jz)-phii_trial(jz); 
                            swt_mod(jz)=swt_trial(jz);
                            swi_mod(jz)=swt_trial(jz);
                            swm_mod(jz)=swt_trial(jz);      
                        end
                    end
                end
                
                %%%%%%%%%%%%%%%%%%% local error update
                if (wloc~=0 && jtemp >= iter_loc && mod(jtemp,freq_loc)==0)
% % %                 if (wloc~=0 && jtemp >= iter_loc)
                    for jr=1:nz %%%% loop over layers
                        for jm=1:npar %%% loop over Micro  and matrix porosity & saturations 
                            if (etrial_local(jr,jm)<=emod_local(jr,jm))%%% local error to update model here and to update local_tmp later if reannealing is on
                                emod_local(jr,jm)=etrial_local(jr,jm);
                                model_mod(jr,jm)=model_trial(jr,jm);
                            else
                                arg=(etrial_local(jr,jm)-emod_local(jr,jm))./tmp_local(jr,jm);
                                if arg > 1.e6
                                    pde=0.001;
                                else
                                    pde=exp(-arg);
                                end

                                    if pde>rand %%% hist_updat
                                        emod_local(jr,jm)=etrial_local(jr,jm);
                                        model_mod(jr,jm)=model_trial(jr,jm);
                                    end                
                            end                        
                        end
                    end
                    phii_mod=model_mod(:,1);
                    phit_mod=model_mod(:,2);
                    phim_mod=phit_mod-phii_mod;
                    swt_mod=model_mod(:,3);  
                    swi_mod=swt_mod;
                    swm_mod=swt_mod;
                end  
                                
% % % %                 %%%%%%%%%%%%% grab global statistics for all jrun,jtemp,jmov
% % % %                 glamdai(jrun,jtemp,jmov)=lamdai_mod; %%%%% Global lamdai for PPD,etc. 
% % % %                 gCx(jrun,jtemp,jmov)=Cx_mod;
% % % %                 gphic(jrun,jtemp,jmov)=phic_mod;
% % % %                 grhos(jrun,jtemp,jmov)=rhos_mod;
% % % %                 gks(jrun,jtemp,jmov)=ks_mod;
% % % %                 gmus(jrun,jtemp,jmov)=mus_mod;
% % % %                 gLik(jrun,jtemp,jmov)=Lik_mod;
% % % %                 gLimu(jrun,jtemp,jmov)=Limu_mod;
% % % %                 
% % % %                 
% % % %                 gphii(jrun,jtemp,jmov,:)=phii_mod; %%%%% Global intergranular porosity for PPD,etc.
% % % %                 gphit(jrun,jtemp,jmov,:)=phit_mod; %%%%% Global total porosity for PPD,etc. 
% % % %                 gphim(jrun,jtemp,jmov,:)=phim_mod; %%%%% Global Micro-porosity for PPD,etc.
% % % %                 gswt(jrun,jtemp,jmov,:)=swt_mod; %%%%% Global total water saturation for PPD,etc. 
% % % %                 gswi(jrun,jtemp,jmov,:)=swi_mod; %%%%% Global intergranular water saturation for PPD,etc. 
% % % %                 gswm(jrun,jtemp,jmov,:)=swm_mod; %%%%% Global Micro  water saturation for PPD,etc.
% % % %                 gphiiw(jrun,jtemp,jmov,:)=phii_mod.*swi_mod;
% % % %                 gphimw(jrun,jtemp,jmov,:)=phim_mod.*swm_mod;
       
                
            end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end move
            
            %%%%%%%%%%%%%%%% for plottig             
            err_plt(jtemp+1,jrun)=emod;               
            lamdai_plt(jtemp+1,jrun)=lamdai_mod;
            Cx_plt(jtemp+1,jrun)=Cx_mod;
            Lik_plt(jtemp+1,jrun)=Lik_mod;
            Limu_plt(jtemp+1,jrun)=Limu_mod;
            phic_plt(jtemp+1,jrun)=phic_mod;
            rhos_plt(jtemp+1,jrun)=rhos_mod;
            ks_plt(jtemp+1,jrun)=ks_mod;
            mus_plt(jtemp+1,jrun)=mus_mod;
            emod_local_plt(jtemp+1,jrun,:,:)=emod_local(:,:);
            
            
            %%%%%%%%%%%%% grab global statistics for all jrun,jtemp
                gemod(jrun,jtemp)=emod;   
                glamdai(jrun,jtemp)=lamdai_mod; %%%%% 
                gCx(jrun,jtemp)=Cx_mod;
                gphic(jrun,jtemp)=phic_mod;
                grhos(jrun,jtemp)=rhos_mod;
                gks(jrun,jtemp)=ks_mod;
                gmus(jrun,jtemp)=mus_mod;
                gLik(jrun,jtemp)=Lik_mod;
                gLimu(jrun,jtemp)=Limu_mod;
                
                
                gphii(jrun,jtemp,:)=phii_mod; %%%%% Global intergranular porosity for PPD,etc.
                gphit(jrun,jtemp,:)=phit_mod; %%%%% Global total porosity for PPD,etc. 
                gphim(jrun,jtemp,:)=phim_mod; %%%%% Global Micro-porosity for PPD,etc.
                gswt(jrun,jtemp,:)=swt_mod; %%%%% Global total water saturation for PPD,etc. 
                gswi(jrun,jtemp,:)=swi_mod; %%%%% Global intergranular water saturation for PPD,etc. 
                gswm(jrun,jtemp,:)=swm_mod; %%%%% Global Micro  water saturation for PPD,etc.
                gphiiw(jrun,jtemp,:)=phii_mod.*swi_mod;
                gphimw(jrun,jtemp,:)=phim_mod.*swm_mod;   
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Local temp update
        if (wloc~=0 && jtemp >= iter_loc && mod(jtemp,freq_loc)==0)
            if reannealing~=0
                deplot1(jtemp)=0.0; dmplot1(jtemp)=0.0;splot1(jtemp)=0.0;tplot1(jtemp)=0.0;rplot1(jtemp)=0.0;
                if (jtemp >= iter_reann && mod(jtemp,freq_reann)==0)
                   [tmp_local,smax,R,deplot,dmplot,splot,tplot,rplot]=temp_scale(tmp_local,emod_local,etrial_local,model_mod,model_trial,...
                                                               model_min,model_max,nz,npar,gama,reannealing);%%% Reanealing
                    deplot1(jtemp)=deplot; dmplot1(jtemp)=dmplot;splot1(jtemp)=splot;tplot1(jtemp)=tplot;rplot1(jtemp)=rplot;
                    jjtemp=1;
                    tmpi_local=tmp_local;
                end
            end
            tmp_local=tmpi_local.*exp(-decay_local.*(jjtemp).^0.5); %%%Annealing
            jjtemp=jjtemp+1;               
                for jr=1:nz %%%% loop over layers
                   for jm=1:npar %%% loop over Micro  and matrix porosity and saturations 
                       if tmp_local(jr,jm)<=temp_min
                           tmp_local(jr,jm)=temp_min;
                       end
                   end
                end             
               local_temp(:,:,jtemp+1)=tmp_local; %%%%%%plotting
        end               
            

            %%%%%%%%%%%%% Exit from temp loop if emod is so small
            if emod<=error_tresh
                jtemp=niter;                
            end
            sprintf('Run#      Iter#     Globerror           Sal(ppk)        Critical Porosity        Lik   Limu   intergranular-lambda  Matrix Density  Bulk Modulus  Shear Modulus')
            sprintf('%5.0f     %5.0f        %12.10f          %3.0f              %5.2f              %5.2f     %5.2f      %5.2f     %5.2f      %5.2f     %5.2f',...
                jrun,jtemp,emod,Cx_mod,phic_mod,Lik_mod,Limu_mod,lamdai_mod,rhos_mod,ks_mod,mus_mod)
            
                      
            jtemp=jtemp+1;  
        end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Compute mean models
phim_mean=zeros(nz,1);
phii_mean=zeros(nz,1);
phit_mean=zeros(nz,1);
swt_mean=zeros(nz,1);
swi_mean=zeros(nz,1);
swm_mean=zeros(nz,1);


for jrun=1:nruns
        vec=gemod(jrun,:);
        gemod_min(jrun)=min(vec);
        inerr1=find(vec==min(vec));
        inerr2=length(inerr1);
        if inerr2 > 1
            inerr=inerr1(inerr2);
        else
            inerr=inerr1;
        end
        lamdai_opt(jrun)=glamdai(jrun,inerr);
        phic_opt(jrun)=gphic(jrun,inerr);
        Cx_opt(jrun)=gCx(jrun,inerr);
        Lik_opt(jrun)=gLik(jrun,inerr);  
        Limu_opt(jrun)=gLimu(jrun,inerr); 
        rhos_opt(jrun)=grhos(jrun,inerr);       
        ks_opt(jrun)=gks(jrun,inerr); 
        mus_opt(jrun)=gmus(jrun,inerr);
        phii_opt(jrun,:)=reshape(gphii(jrun,inerr,:),[nz 1]);
        phit_opt(jrun,:)=reshape(gphit(jrun,inerr,:),[nz 1]);     
        swt_opt(jrun,:)=reshape(gswt(jrun,inerr,:),[nz 1]);
end

phim_opt=phit_opt-phii_opt;
swi_opt=swt_opt;
swm_opt=swt_opt;
phiiw_opt=phii_opt.*swi_opt;
phimw_opt=phim_opt.*swm_opt;
sot_opt=1-swt_opt;
sgt_opt=1-swt_opt-sot_opt;  


vec1=(1-gemod_min);
sumvec=sum(vec1);
lamdai_mean=dot(vec1,lamdai_opt)/sumvec;  
phic_mean=dot(vec1,phic_opt)/sumvec;
Cx_mean=dot(vec1,Cx_opt)/sumvec; 
Lik_mean=dot(vec1,Lik_opt)/sumvec;  
Limu_mean=dot(vec1,Limu_opt)/sumvec;
rhos_mean=dot(vec1,rhos_opt)/sumvec; 
ks_mean=dot(vec1,ks_opt)/sumvec;
mus_mean=dot(vec1,mus_opt)/sumvec;  

for jz=1:nz
    phii_mean(jz)=dot(vec1(:),phii_opt(:,jz))./sumvec;
    phit_mean(jz)=dot(vec1(:),phit_opt(:,jz))./sumvec;
    swt_mean(jz)=dot(vec1(:),swt_opt(:,jz))./sumvec;
end    

phim_mean=phit_mean-phii_mean;
swi_mean=swt_mean;
swm_mean=swt_mean;
phiiw_mean=phii_mean.*swi_mean;
phimw_mean=phim_mean.*swm_mean;
sot_mean=1-swt_mean;
sgt_mean=1-swt_mean-sot_mean;   

        
    
    
% % % if nruns==1
% % %     
% % %     gemod_min=min(gemod);
% % %     inerr1=find(gemod==gemod_min);
% % %     inerr2=length(inerr1);
% % %     if inerr2 > 1
% % %         inerr=inerr1(inerr2);
% % %     else
% % %         inerr=inerr1;
% % %     end
% % %     lamdai_mean=glamdai(inerr-1);
% % %     Cx_mean=gCx(inerr-1);
% % %     phic_mean=gphic(inerr-1);
% % %     Lik_mean=gLik(inerr-1);
% % %     Limu_mean=gLimu(inerr-1);
% % %     ks_mean=gks(inerr-1);
% % %     mus_mean=gmus(inerr-1);
% % %     rhos_mean=grhos(inerr-1);  
% % %     phii_mean=reshape(gphii(1,inerr-1,:),[np 1]);
% % %     phit_mean=reshape(gphit(1,inerr-1,:),[np 1]);
% % %     phim_mean=phit_mean-phii_mean;
% % %     swt_mean=reshape(gswt(1,inerr-1,:),[np 1]);
% % %     swi_mean=swt_mean;
% % %     swm_mean=swt_mean;
% % %     phiiw_mean=phii_mean.*swi_mean;
% % %     phimw_mean=phim_mean.*swm_mean;
% % %     sot_mean=1-swt_mean;
% % %     sgt_mean=1-swt_mean-sot_mean;   
% % % 
% % % end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Compute smooth mean models
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% for continous variables
phii_smth=zeros(size(phii_mean)); phit_smth=zeros(size(phit_mean)); swt_smth=zeros(size(swi_mean)); 

phii_smth(1)=(phii_mean(1)+phii_mean(2))/2;
phii_smth(2)=(phii_mean(2)+phii_mean(3))/2;
for jz=3:nz-2
    phii_smth(jz)=(phii_mean(jz)+phii_mean(jz-1)+phii_mean(jz+1))/3;
end
phii_smth(nz-1)=(phii_mean(nz-1)+phii_mean(nz-2))/2;
phii_smth(nz)=(phii_mean(nz)+phii_mean(nz-1))/2;

phit_smth(1)=(phit_mean(1)+phit_mean(2))/2;
phit_smth(2)=(phit_mean(2)+phit_mean(3))/2;
for jz=3:nz-2
    phit_smth(jz)=(phit_mean(jz)+phit_mean(jz-1)+phit_mean(jz+1))/3;
end
phit_smth(nz-1)=(phit_mean(nz-1)+phit_mean(nz-2))/2;
phit_smth(nz)=(phit_mean(nz)+phit_mean(nz-1))/2;


phim_smth=phit_smth-phii_smth;



swt_smth(1)=(swt_mean(1)+swt_mean(2))/2;
swt_smth(2)=(swt_mean(2)+swt_mean(3))/2;
for jz=3:nz-2
    swt_smth(jz)=(swt_mean(jz)+swt_mean(jz-1)+swt_mean(jz+1))/3;
end
swt_smth(nz-1)=(swt_mean(nz-1)+swt_mean(nz-2))/2;
swt_smth(nz)=(swt_mean(nz)+swt_mean(nz-1))/2;

swi_smth=swt_smth;
swm_smth=swt_smth;



phiiw_smth=phii_smth.*swi_smth;
phimw_smth=phim_smth.*swm_smth;

sot_smth=1-swt_smth;
sgt_smth=1-swt_smth-sot_smth;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% compute simulated data with PPD
%%%%%%%%%%%%%%% Velocities  

if nruns>1
    
    for jrun=1:nruns
        [dtp_sim_ppd(:,jrun),dts_sim_ppd(:,jrun),vp_sim_ppd(:,jrun),vs_sim_ppd(:,jrun),ai_sim_ppd(:,jrun),si_sim_ppd(:,jrun)...
            ,vr_sim_ppd(:,jrun),pr_sim_ppd(:,jrun),ksat_sim_ppd(:,jrun),musat_sim_ppd(:,jrun),den_sim_ppd(:,jrun),logk,logmu]=...
             sdem_sonic_carbonate_v01...
             (phit_opt(jrun,:),phii_opt(jrun,:),phim_opt(jrun,:),phic_opt(jrun),Lik_opt(jrun),Limu_opt(jrun),Lmk,Lmmu,...
             swt_opt(jrun,:),sot_opt(jrun,:),sgt_opt(jrun,:),satype,kw,k_oil,k_gas,rhow,rho_oil,rho_gas,ks_opt(jrun),mus_opt(jrun),rhos_opt(jrun));


    %%%%%%%%%%%%%% Resistivity
        [resist_sim_ppd(:,jrun),ffac_sim_ppd(:,jrun),rw_sim_ppd(:,jrun),logr]...
            =sdem_resist_carbonate(phit_opt(jrun,:),phii_opt(jrun,:),phim_opt(jrun,:),swi_opt(jrun,:),swm_opt(jrun,:),lamdai_opt(jrun),lamdam,tmpc,Cx_opt(jrun),wetcase);

    end
    
    for jz=1:nz
        dtp_sim_mean(jz,1)=mean(dtp_sim_ppd(jz,:));
        dts_sim_mean(jz,1)=mean(dts_sim_ppd(jz,:));
        den_sim_mean(jz,1)=mean(den_sim_ppd(jz,:));
        resist_sim_mean(jz,1)=mean(resist_sim_ppd(jz,:));
    end
    
else
        
    [dtp_sim_mean,dts_sim_mean,vp_sim_mean,vs_sim_mean,ai_sim_mean,si_sim_mean,vr_sim_mean,pr_sim_mean,ksat_sim_mean,musat_sim_mean,den_sim_mean,logk,logmu]=...
             sdem_sonic_carbonate_v01...
             (phit_mean,phii_mean,phim_mean,phic_mean,Lik_mean,Limu_mean,Lmk,Lmmu,...
             swt_mean,sot_mean,sgt_mean,satype,kw,k_oil,k_gas,rhow,rho_oil,rho_gas,ks_mean,mus_mean,rhos_mean);

    %%%%%%%%%%%%%% Resistivity
    [resist_sim_mean,ffac_sim_mean,rw_sim_mean,logr]=sdem_resist_carbonate(phit_mean,phii_mean,phim_mean,swi_mean,swm_mean,lamdai_mean,lamdam,tmpc,Cx_mean,wetcase);
   
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% compute simulated data with
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% smooth mean models for continous variables 
%%%%%%%%%%%%%%% Velocities  
[dtp_sim_smth,dts_sim_smth,vp_sim_smth,vs_sim_smth,ai_sim_smth,si_sim_smth,vr_sim_smth,pr_sim_smth,ksat_sim_smth,musat_sim_smth,den_sim_smth,logk,logmu]=...
         sdem_sonic_carbonate_v01...
         (phit_smth,phii_smth,phim_smth,phic_mean,Lik_mean,Limu_mean,Lmk,Lmmu,...
         swt_smth,sot_smth,sgt_smth,satype,kw,k_oil,k_gas,rhow,rho_oil,rho_gas,ks_mean,mus_mean,rhos_mean);

%%%%%%%%%%%%% Resistivity
[resist_sim_smth,ffac_sim_smth,rw_sim_smth,logr]=sdem_resist_carbonate(phit_smth,phii_smth,phim_smth,swi_smth,swm_smth,lamdai_mean,lamdam,tmpc,Cx_mean,wetcase);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%plotting
nitery=length(err_plt);
lamdai_true(1:nitery,1)=lamdai;
Cx_true(1:nitery,1)=Cx;
phic_true(1:nitery,1)=phic;
rhos_true(1:nitery,1)=rho_ma;
ks_true(1:nitery,1)=k_ma;
mus_true(1:nitery,1)=mu_ma;
Lik_true(1:nitery,1)=Lik;
Limu_true(1:nitery,1)=Limu;
lamdai_mean_line(1:nitery,1)=lamdai_mean;
Cx_mean_line(1:nitery,1)=Cx_mean;
phic_mean_line(1:nitery,1)=phic_mean;
rhos_mean_line(1:nitery,1)=rhos_mean;
ks_mean_line(1:nitery,1)=ks_mean;
mus_mean_line(1:nitery,1)=mus_mean;
Lik_mean_line(1:nitery,1)=Lik_mean;
Limu_mean_line(1:nitery,1)=Limu_mean;


krun=1:nruns;
for jtemp=1:niter
    err_mean_line(jtemp)=mean(err_plt(jtemp,:));
    lamdai_plt_mean(jtemp)=mean(lamdai_plt(jtemp,:));
    phic_plt_mean(jtemp)=mean(phic_plt(jtemp,:));
    Lik_plt_mean(jtemp)=mean(Lik_plt(jtemp,:));
    Limu_plt_mean(jtemp)=mean(Limu_plt(jtemp,:));
    Cx_plt_mean(jtemp)=mean(Cx_plt(jtemp,:));
    ks_plt_mean(jtemp)=mean(ks_plt(jtemp,:));
    mus_plt_mean(jtemp)=mean(mus_plt(jtemp,:));
    rhos_plt_mean(jtemp)=mean(rhos_plt(jtemp,:));
end

    err_true=zeros(size(err_mean_line));

    figure(14)
    runy=1:nruns;
    subplot(3,3,1)  
    semilogx(err_plt(:,runy),'c','LineWidth',lw)
    hold on
    semilogx(err_true,'b','LineWidth',lw)
    hold on
    semilogx(err_mean_line,'r','LineWidth',lw)
    xlim([0 jtemp-1]); grid;
    ylabel('Objective function','FontSize',fs);    title('Objective function','FontSize',fs);
    set(gca,'FontSize',fs)
    subplot(3,3,2)
    semilogx(lamdai_plt(:,runy),'c','LineWidth',lw); 
    hold on
    semilogx(lamdai_true,'b','LineWidth',lw)
    hold on
    semilogx(lamdai_plt_mean,'r','LineWidth',lw)
    xlim([0 jtemp-1]);
    ylabel('Model parameters','FontSize',fs);    
    title('Intergranular Lithology Exponent','FontSize',fs);
    ylim(lamdai_lim);grid; set(gca,'FontSize',fs) 
    subplot(3,3,3)
    semilogx(phic_plt(:,runy),'c','LineWidth',lw); 
    hold on
    semilogx(phic_true,'b','LineWidth',lw)
    hold on
    semilogx(phic_plt_mean,'r','LineWidth',lw)
    xlim([0 jtemp-1]);
    ylabel('Model parameters','FontSize',fs);   
    title('Critical Porosity(v/v)','FontSize',fs);
    ylim(phic_lim); grid; set(gca,'FontSize',fs) 
    subplot(3,3,4)
    semilogx(Cx_plt(:,runy),'c','LineWidth',lw); 
    hold on
    semilogx(Cx_true,'b','LineWidth',lw)
    hold on
    semilogx(Cx_plt_mean,'r','LineWidth',lw)
    xlim([0 jtemp-1]);
    ylabel('Model parameters','FontSize',fs);    title('Salinity (Kppm)','FontSize',fs);
    ylim(Cx_lim); grid; set(gca,'FontSize',fs) 
    subplot(3,3,5)
    semilogx(Lik_plt(:,runy),'c','LineWidth',lw);
    hold on
    semilogx(Lik_true,'b','LineWidth',lw)
    hold on
    semilogx(Lik_plt_mean,'r','LineWidth',lw)
    xlim([0 jtemp-1]);
    title('Intergranular Length Scale for Bulk Modulus ','FontSize',fs);
    ylim(Lik_lim);grid; set(gca,'FontSize',fs) 
    subplot(3,3,6)
    semilogx(Limu_plt(:,runy),'c','LineWidth',lw);
    hold on
    semilogx(Limu_true,'b','LineWidth',lw)
     hold on
    semilogx(Limu_plt_mean,'r','LineWidth',lw)
    xlim([0 jtemp-1]);
    title('Intergranular Length Scale for Shear Modulus ','FontSize',fs);
    ylim(Limu_lim); grid;  set(gca,'FontSize',fs) 
    subplot(3,3,7)
    semilogx(rhos_plt(:,runy),'c','LineWidth',lw);
    hold on
    semilogx(rhos_true,'b','LineWidth',lw)
    hold on
    semilogx(rhos_plt_mean,'r','LineWidth',lw)
    xlim([0 jtemp-1]);
    xlabel('Iteration number','FontSize',fs);    
    ylabel('Model parameters','FontSize',fs);    title('Matrix Density(gr/cm3) ','FontSize',fs);
    ylim(rhos_lim); grid;  set(gca,'FontSize',fs) 
    subplot(3,3,8)
    semilogx(ks_plt(:,runy),'c','LineWidth',lw);
    hold on
    semilogx(ks_true,'b','LineWidth',lw)
    hold on
    semilogx(ks_plt_mean,'r','LineWidth',lw)
    xlim([0 jtemp-1]);
    xlabel('Iteration number','FontSize',fs);    
    title('Matrix Bulk Modulus(GPa) ','FontSize',fs);
    ylim(ks_lim); grid;  set(gca,'FontSize',fs) 
    subplot(3,3,9)
    semilogx(mus_plt(:,runy),'c','LineWidth',lw);
    hold on
    semilogx(mus_true,'b','LineWidth',lw)
    hold on
    semilogx(mus_plt_mean,'r','LineWidth',lw)
    xlim([0 jtemp-1]);
    xlabel('Iteration number','FontSize',fs);    
    title('Matrix Shear Modulus(GPa) ','FontSize',fs);
    ylim(mus_lim); grid;  set(gca,'FontSize',fs) 


    
if (nruns>1)
    nhist=15;
    figure(15)
    subplot(3,3,1) 
    hist(err_plt(niter,:),nhist)
    xlim([-max(err_plt(niter,:)) max(err_plt(niter,:))]);
    grid; xlabel('Global error','FontSize',fs);ylabel('Frequency','FontSize',fs); set(gca,'FontSize',fs)
    subplot(3,3,2) 
    hist(lamdai_opt,nhist)
    xlim(lamdai_lim);grid; xlabel('Intergranular lithology exponent','FontSize',fs); 
    ylabel('Frequency','FontSize',fs);set(gca,'FontSize',fs)
    subplot(3,3,3) 
    hist(phic_opt,nhist)
    xlim(phic_lim);grid; xlabel('Critical porosity(v/v)','FontSize',fs); ylabel('Frequency','FontSize',fs);
    set(gca,'FontSize',fs)
    subplot(3,3,4) 
    hist(Cx_opt,nhist)
    xlim(Cx_lim);grid; xlabel('Salinity(Kppm)','FontSize',fs); ylabel('Frequency','FontSize',fs);
    set(gca,'FontSize',fs)
    subplot(3,3,5) 
    hist(Lik_opt,nhist)
    xlim(Lik_lim); grid; xlabel('Intergranular L parameter for bulk modulus','FontSize',fs); 
    ylabel('Frequency','FontSize',fs);  set(gca,'FontSize',fs)
    subplot(3,3,6) 
    hist(Limu_opt,nhist)
    xlim(Limu_lim); grid; xlabel('Intergranular L parameter for shear modulus','FontSize',fs); 
    ylabel('Frequency','FontSize',fs); set(gca,'FontSize',fs)
    subplot(3,3,7) 
    hist(rhos_opt,nhist)
    xlim(rhos_lim); grid; xlabel('Matrix Density(gr/cm3)','FontSize',fs); 
    ylabel('Frequency','FontSize',fs); set(gca,'FontSize',fs)
    subplot(3,3,8) 
    hist(ks_opt,nhist)
    xlim(ks_lim); grid; xlabel('Matrix Bulk Modulus(GPa)','FontSize',fs); 
    ylabel('Frequency','FontSize',fs); set(gca,'FontSize',fs)
    subplot(3,3,9) 
    hist(mus_opt,nhist)
    xlim(mus_lim); grid; xlabel('Matrix Shear Modulus(GPa)','FontSize',fs); 
    ylabel('Frequency','FontSize',fs); set(gca,'FontSize',fs)
end
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% compare data and model (mean model no smoothing)
figure(16)
subplot(1,4,1)
semilogx(resist,DEPTH,'b',resist_sim_mean,DEPTH,'r','LineWidth',lw)
ylabel('Relative Depth(m)','FontSize',fsy);
legend('True','Inverted')
ylim(depth_lim_all);
xlim(resist_xlim);
title('Resistivity(ohm-m)','FontSize',fst);grid MINOR;grid on;set(gca,'YDir','reverse','FontSize',fset)
subplot(1,4,3)
plot(dtco,DEPTH,'b',dtp_sim_mean,DEPTH,'r','LineWidth',lw)
title('DTCO(microsec/m)','FontSize',fst);legend('True','Inverted');xlim(dtp_xlim);ylim(depth_lim_all);
grid MINOR;grid on;set(gca,'YDir','reverse','FontSize',fset)
subplot(1,4,4)
plot(dtsm,DEPTH,'b',dts_sim_mean,DEPTH,'r','LineWidth',lw)
title('DTSM(microsec/m)','FontSize',fst);legend('True','Inverted');xlim(dts_xlim);ylim(depth_lim_all);
grid MINOR;grid on;set(gca,'YDir','reverse','FontSize',fset)
subplot(1,4,2)
plot(den,DEPTH,'b',den_sim_mean,DEPTH,'r','LineWidth',lw)
title('Density(gr/cm3)','FontSize',fst);legend('True','Inverted');xlim(den_xlim);ylim(depth_lim_all);
grid MINOR;grid on;set(gca,'YDir','reverse','FontSize',fset)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% compare data and model (PPD)
if nruns>1
    nrunny=1:nruns;
    figure(17)
    subplot(1,4,1)
    semilogx(resist_sim_ppd(:,nrunny),DEPTH,'c','LineWidth',lw)
    hold on
    semilogx(resist,DEPTH,'b','LineWidth',lw)
    hold on
    semilogx(resist_sim_mean,DEPTH,'r','LineWidth',lw)
    ylabel('Relative Depth(m)','FontSize',fsy);
    ylim(depth_lim_all);
    xlim(resist_xlim);
    title('Resistivity(ohm-m)','FontSize',fst);grid MINOR;grid on;set(gca,'YDir','reverse','FontSize',fset)
    subplot(1,4,3)
    plot(dtp_sim_ppd(:,nrunny),DEPTH,'c','LineWidth',lw)
    hold on
    plot(dtco,DEPTH,'b','LineWidth',lw)
    hold on
    plot(dtp_sim_mean,DEPTH,'r','LineWidth',lw)
    title('DTCO(microsec/m)','FontSize',fst);
    xlim(dtp_xlim);ylim(depth_lim_all);
    grid MINOR;grid on;set(gca,'YDir','reverse','FontSize',fset)
    subplot(1,4,4)
    plot(dts_sim_ppd(:,nrunny),DEPTH,'c','LineWidth',lw)
    hold on
    plot(dtsm,DEPTH,'b','LineWidth',lw)
    hold on
    plot(dts_sim_mean,DEPTH,'r','LineWidth',lw)
    title('DTSM(microsec/m)','FontSize',fst);
    xlim(dts_xlim);ylim(depth_lim_all);
    grid MINOR;grid on;set(gca,'YDir','reverse','FontSize',fset)
    subplot(1,4,2)
    plot(den_sim_ppd(:,nrunny),DEPTH,'c','LineWidth',lw)
    hold on
    plot(den,DEPTH,'b','LineWidth',lw)
    hold on
    plot(den_sim_mean,DEPTH,'r','LineWidth',lw)
    title('Density(gr/cm3)','FontSize',fst);
    xlim(den_xlim);ylim(depth_lim_all);
    grid MINOR;grid on;set(gca,'YDir','reverse','FontSize',fset)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% compare data and model (smoothing)
figure(18)
subplot(1,4,1)
semilogx(resist,DEPTH,'b',resist_sim_smth,DEPTH,'r','LineWidth',lw)
ylabel('Relative Depth(m)','FontSize',fsy);
legend('True','Smoothed inverted')
ylim(depth_lim_all);
xlim(resist_xlim);
title('Resistivity(ohm-m)','FontSize',fst);grid on;set(gca,'YDir','reverse','FontSize',fset)
subplot(1,4,3)
plot(dtco,DEPTH,'b',dtp_sim_smth,DEPTH,'r','LineWidth',lw)
title('DTCO(microsec/m)','FontSize',fst);legend('True','Smoothed inverted');xlim(dtp_xlim);ylim(depth_lim_all);
grid MINOR;grid on;set(gca,'YDir','reverse','FontSize',fset)
subplot(1,4,4)
plot(dtsm,DEPTH,'b',dts_sim_smth,DEPTH,'r','LineWidth',lw)
title('DTSM(microsec/m)','FontSize',fst);legend('True','Smoothed inverted');xlim(dts_xlim);ylim(depth_lim_all);
grid MINOR;grid on;set(gca,'YDir','reverse','FontSize',fset)
subplot(1,4,2)
plot(den,DEPTH,'b',den_sim_smth,DEPTH,'r','LineWidth',lw)
title('Density(gr/cm3)','FontSize',fst);legend('True','Smoothed inverted');
xlim(den_xlim);ylim(depth_lim_all);
grid MINOR;grid on;set(gca,'YDir','reverse','FontSize',fset)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plotting models
figure(19)
subplot(1,4,1)
plot(phii_mean,DEPTH,'r',phii,DEPTH,'b','LineWidth',lw)
legend('Inverted','True');ylabel('Relative Depth(m)','FontSize',fsy);title('Intergranular porosity(v/v)','FontSize',fst)
xlim(phi_xlim);ylim(depth_lim_all);grid MINOR;grid on;set(gca,'YDir','reverse','FontSize',fset)       
subplot(1,4,2)
plot(phim_mean,DEPTH,'r',phim,DEPTH,'b','LineWidth',lw)
legend('Inverted','True');title('Micro-porosity(v/v)','FontSize',fst)
xlim(phi_xlim);ylim(depth_lim_all);grid MINOR;grid on;set(gca,'YDir','reverse','FontSize',fset)       
subplot(1,4,3)
plot(phit,DEPTH,'k',phit_mean,DEPTH,':k',phiiw,DEPTH,'b',phiiw_mean,DEPTH,'r','LineWidth',lw)
title('Porosity(v/v)','FontSize',fst);legend('True total','Inverted total','True Intg. water-filled','Inverted Intg. water-filled')
xlim(phi_xlim);ylim(depth_lim_all);grid MINOR;grid on;set(gca,'YDir','reverse','FontSize',fset)       
subplot(1,4,4)
plot(swt_mean,DEPTH,'r',swt,DEPTH,'b','LineWidth',lw)
title('Total water saturation(v/v)','FontSize',fst);legend('Inverted','True')
xlim(sw_xlim);ylim(depth_lim_all);grid MINOR;grid on;set(gca,'YDir','reverse','FontSize',fset)



figure(20)
subplot(1,7,1)
plot(phii_mean,DEPTH,'r',phii,DEPTH,'b','LineWidth',lw)
legend('Inverted','True');ylabel('Relative Depth(m)','FontSize',fsy);title('Intergranular porosity(v/v)','FontSize',fst)
xlim(phi_xlim);ylim(depth_lim_all);grid MINOR;grid on;set(gca,'YDir','reverse','FontSize',fset)       
subplot(1,7,2)
plot(phim_mean,DEPTH,'r',phim,DEPTH,'b','LineWidth',lw)
legend('Inverted','True');title('Micro-porosity(v/v)','FontSize',fst)
xlim(phi_xlim);ylim(depth_lim_all);grid MINOR;grid on;set(gca,'YDir','reverse','FontSize',fset)       
subplot(1,7,3)
plot(phit,DEPTH,'k',phit_mean,DEPTH,':k',phiiw,DEPTH,'b',phiiw_mean,DEPTH,'r','LineWidth',lw)
title('Porosity(v/v)','FontSize',fst);legend('True total','Inverted total','True Intg. water-filled','Inverted Intg. water-filled')
xlim(phi_xlim);ylim(depth_lim_all);grid MINOR;grid on;set(gca,'YDir','reverse','FontSize',fset)       
subplot(1,7,4)
semilogx(resist,DEPTH,'b',resist_sim_mean,DEPTH,'r','LineWidth',lw)
ylabel('Relative Depth(m)','FontSize',fsy);
legend('True','Inverted')
ylim(depth_lim_all);
xlim(resist_xlim);
title('Resistivity(ohm-m)','FontSize',fst);grid MINOR;grid on;set(gca,'YDir','reverse','FontSize',fset)
subplot(1,7,5)
plot(den,DEPTH,'b',den_sim_mean,DEPTH,'r','LineWidth',lw)
title('Density(gr/cm3)','FontSize',fst);legend('True','Inverted');xlim(den_xlim);ylim(depth_lim_all);
grid MINOR;grid on;set(gca,'YDir','reverse','FontSize',fset)
subplot(1,7,6)
plot(dtco,DEPTH,'b',dtp_sim_mean,DEPTH,'r','LineWidth',lw)
title('DTCO(microsec/m)','FontSize',fst);legend('True','Inverted');xlim(dtp_xlim);ylim(depth_lim_all);
grid MINOR;grid on;set(gca,'YDir','reverse','FontSize',fset)
subplot(1,7,7)
plot(dtsm,DEPTH,'b',dts_sim_mean,DEPTH,'r','LineWidth',lw)
title('DTSM(microsec/m)','FontSize',fst);legend('True','Inverted');xlim(dts_xlim);ylim(depth_lim_all);
grid MINOR;grid on;set(gca,'YDir','reverse','FontSize',fset)





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(21)
subplot(1,4,1)
plot(phii_smth,DEPTH,'r',phii,DEPTH,'b','LineWidth',lw)
legend('Smoothed inverted','True');ylabel('Relative Depth(m)','FontSize',fsy);title('Intergranular porosity(v/v)','FontSize',fst)
xlim(phi_xlim);ylim(depth_lim_all);grid MINOR;grid on;set(gca,'YDir','reverse','FontSize',fset)       
subplot(1,4,2)
plot(phim_smth,DEPTH,'r',phim,DEPTH,'b','LineWidth',lw)
legend('Smoothed inverted','True');title('Micro-porosity(v/v)','FontSize',fst)
xlim(phi_xlim);ylim(depth_lim_all);grid MINOR;grid on;set(gca,'YDir','reverse','FontSize',fset)       
subplot(1,4,3)
plot(phit,DEPTH,'k',phit_smth,DEPTH,':k',phiiw,DEPTH,'b',phiiw_smth,DEPTH,'r','LineWidth',lw)
title('Porosity(v/v)','FontSize',fst);legend('True total','Smoothed inverted total','True Intg. water-filled','Smoothed inverted Intg. water-filled')
xlim(phi_xlim);ylim(depth_lim_all);grid MINOR;grid on;set(gca,'YDir','reverse','FontSize',fset)       
subplot(1,4,4)
plot(swt_smth,DEPTH,'r',swt,DEPTH,'b','LineWidth',lw)
title('Water saturation(v/v)','FontSize',fst);legend('Smoothed inverted','True')
xlim(sw_xlim);ylim(depth_lim_all);grid MINOR;grid on;set(gca,'YDir','reverse','FontSize',fset)



if nruns>1
        runy=1:nruns;
        figure(22)
        subplot(1,4,1)
        plot(phii_opt(runy,:),DEPTH,'c','LineWidth',lw)
        hold on
        plot(phii_mean,DEPTH,'r',phii,DEPTH,'b','LineWidth',lw)
        title('Intergranular porosity(v/v)','FontSize',fst)
        ylabel('Relative Depth(m)','FontSize',fsy)
        xlim(phi_xlim); ylim(depth_lim_all);grid MINOR; grid on; set(gca,'YDir','reverse','FontSize',fset)       
        subplot(1,4,2)
        plot(phim_opt(runy,:),DEPTH,'c','LineWidth',lw)
        hold on
        plot(phim_mean,DEPTH,'r',phim,DEPTH,'b','LineWidth',lw)
        title('Micro-porosity(v/v)','FontSize',fst)
        xlim(phi_xlim);ylim(depth_lim_all);grid MINOR ;grid on;set(gca,'YDir','reverse','FontSize',fset)       
        subplot(1,4,3)
        plot(phiiw_opt(runy,:),DEPTH,'c',phit_opt(runy,:),DEPTH,'y','LineWidth',lw)
        hold on
        plot(phiiw_mean,DEPTH,'r',phiiw,DEPTH,'b','LineWidth',lw)
        hold on
        plot(phit_mean,DEPTH,':k',phit,DEPTH,'k','LineWidth',lw)
        title('Porosity(v/v)','FontSize',fst)
        xlim(phi_xlim);ylim(depth_lim_all);grid MINOR;grid on;set(gca,'YDir','reverse','FontSize',fset)       
        subplot(1,4,4)
        plot(swt_opt(runy,:),DEPTH,'c','LineWidth',lw)
        hold on
        plot(swt_mean,DEPTH,'r',swt,DEPTH,'b','LineWidth',lw)
        title('Total water saturation(v/v)','FontSize',fst)
        xlim(sw_xlim);ylim(depth_lim_all);grid MINOR;grid on;set(gca,'YDir','reverse','FontSize',fset)

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting PPD for all runs 
if nruns>1
        nd=30;
        dz=1; %%% one feet 
        dd=(1:1:nz).*dz;
        aphiimax=max((max(phii_opt))); 
        aphiimin=min(min(phii_opt));  
        aphiilog=linspace(aphiimin,aphiimax,nd);

        bphimmax=max((max(phim_opt))); 
        bphimmin=min(min(phim_opt));
        bphimlog=linspace(bphimmin,bphimmax,nd);

        phitw_opt=phiiw_opt+phimw_opt;
        cphitwmax=max((max(phitw_opt))); 
        cphiimin=min(min(phitw_opt));  
        cphitwlog=linspace(cphiimin,cphitwmax,nd);

        dswtmax=max((max(swt_opt))); 
        dswtmin=min(min(swt_opt));  
        dswtlog=linspace(dswtmin,dswtmax,nd);

        for jz=1:nz
            phii_ppd(jz,:)=hist(phii_opt(:,jz),aphiilog)./nruns;
            phim_ppd(jz,:)=hist(phim_opt(:,jz),bphimlog)./nruns;
            phitw_ppd(jz,:)=hist(phitw_opt(:,jz),cphitwlog)./nruns;
            swt_ppd(jz,:)=hist(swt_opt(:,jz),dswtlog)./nruns;
        end

        phii_ppd=phii_ppd./max(max(phii_ppd));
        phim_ppd=phim_ppd./max(max(phim_ppd));
        phitw_ppd=phitw_ppd./max(max(phitw_ppd));
        swt_ppd=swt_ppd./max(max(swt_ppd));


        
        
        figure(23)
        subplot(1,4,1);
        imagesc(aphiilog,dd,phii_ppd);
        colorbar;colormap cool;hold on 
        plot(phii,DEPTH,'b','LineWidth',lw)
        hold on
        plot(phii_mean,DEPTH,'r','LineWidth',lw)
        title('Intergranular porosity(v/v)','FontSize',fst);ylabel('Relative Depth(m)','FontSize',fsy)
        xlim(phi_xlim);
        ylim(depth_lim_all);grid MINOR; grid on; set(gca,'YDir','reverse','FontSize',fset)    
        subplot(1,4,2);
        imagesc(bphimlog,dd,phim_ppd);
        colorbar;colormap cool;hold on 
        plot(phim,DEPTH,'b','LineWidth',lw)
        hold on
        plot(phim_mean,DEPTH,'r','LineWidth',lw)
        title('Micro-porosity(v/v)','FontSize',fst)
        xlim(phi_xlim);
        ylim(depth_lim_all);grid MINOR; grid on; set(gca,'YDir','reverse','FontSize',fset)    
        subplot(1,4,3);
        imagesc(cphitwlog,dd,phitw_ppd);
        colorbar;colormap cool;hold on 
        plot(phiiw,DEPTH,'b',phit,DEPTH,'k','LineWidth',lw)
        hold on
        plot(phiiw_mean,DEPTH,'r',phit_mean,DEPTH,':k','LineWidth',lw)
        legend('True Intg. water-filled','True total');title('Porosity(v/v)','FontSize',fst)
        xlim(phi_xlim);
        ylim(depth_lim_all);grid MINOR; grid on; set(gca,'YDir','reverse','FontSize',fset)       
        subplot(1,4,4);
        imagesc(dswtlog,dd,swt_ppd);
        colorbar;colormap cool;hold on 
        plot(swt,DEPTH,'b','LineWidth',lw)
        hold on
        plot(swt_mean,DEPTH,'r','LineWidth',lw)
        title('Total water saturation(v/v)','FontSize',fst)
        xlim(sw_xlim);
        ylim(depth_lim_all);grid MINOR; grid on; set(gca,'YDir','reverse','FontSize',fset) 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting PPD for all runs, all temp, all moves  
% % % nd=30;
% % % dz=1; %%% one feet 
% % % dd=(1:1:np).*dz;
% % % aphiimax=max(max(max((max(gphii))))); 
% % % aphiimin=min(min(min(min(gphii))));  
% % % aphiilog=linspace(aphiimin,aphiimax,nd);
% % % 
% % % kn=nruns*(niter-1)*nmov;
% % % 
% % % bphimmax=max(max(max((max(gphim))))); 
% % % bphimmin=min(min(min(min(gphim))));
% % % bphimlog=linspace(bphimmin,bphimmax,nd);
% % % 
% % % gphitw=gphiiw+gphimw;
% % % cphitwmax=max(max(max((max(gphitw))))); 
% % % cphiimin=min(min(min(min(gphitw))));  
% % % cphitwlog=linspace(cphiimin,cphitwmax,nd);
% % % 
% % % dswtmax=max(max(max((max(gswt))))); 
% % % dswtmin=min(min(min(min(gswt))));  
% % % dswtlog=linspace(dswtmin,dswtmax,nd);
% % % 
% % % for jp=1:np
% % %     gphii_reshape=reshape(gphii(:,:,:,jp),[kn,1]); gphii_ppd(jp,:)=hist(gphii_reshape,aphiilog)./kn;
% % %     gphim_reshape=reshape(gphim(:,:,:,jp),[kn,1]); gphim_ppd(jp,:)=hist(gphim_reshape,bphimlog)./kn;
% % %     gphitw_reshape=reshape(gphitw(:,:,:,jp),[kn,1]); gphitw_ppd(jp,:)=hist(gphitw_reshape,cphitwlog)./kn;
% % %     gswt_reshape=reshape(gswt(:,:,:,jp),[kn,1]); gswt_ppd(jp,:)=hist(gswt_reshape,dswtlog)./kn;
% % % end
% % % 
% % % gphii_ppd=gphii_ppd./max(max(gphii_ppd));
% % % gphim_ppd=gphim_ppd./max(max(gphim_ppd));
% % % gphitw_ppd=gphitw_ppd./max(max(gphitw_ppd));
% % % gswt_ppd=gswt_ppd./max(max(gswt_ppd));
% % % 
% % % 
% % % 
% % % 
% % % figure(6)
% % % subplot(1,4,1);
% % % imagesc(aphiilog,dd,gphii_ppd);
% % % colorbar;colormap cool;hold on 
% % % plot(phii,DEPTH,'b','LineWidth',lw)
% % % hold on
% % % plot(phii_mean,DEPTH,'r','LineWidth',lw)
% % % title('Intergranular porosity(v/v)','FontSize',fst);ylabel('Relative Depth(m)','FontSize',fsy)
% % % xlim(phi_xlim);
% % % ylim(depth_lim_all);grid MINOR; grid on; set(gca,'YDir','reverse','FontSize',fset)    
% % % subplot(1,4,2);
% % % imagesc(bphimlog,dd,gphim_ppd);
% % % colorbar;colormap cool;hold on 
% % % plot(phim,DEPTH,'b','LineWidth',lw)
% % % hold on
% % % plot(phim_mean,DEPTH,'r','LineWidth',lw)
% % % xlim(phi_xlim);
% % % title('Micro-porosity(v/v)','FontSize',fst); ylim(depth_lim_all);grid MINOR; grid on; set(gca,'YDir','reverse','FontSize',fset)    
% % % subplot(1,4,3);
% % % imagesc(cphitwlog,dd,gphitw_ppd);
% % % colorbar;colormap cool;hold on 
% % % plot(phimw+phiiw,DEPTH,'b',phit,DEPTH,'k','LineWidth',lw)
% % % hold on
% % % plot(phimw_mean+phiiw_mean,DEPTH,'r',phit_mean,DEPTH,':k','LineWidth',lw)
% % % legend('True water-filled','True total');title('Porosity(v/v)','FontSize',fst)
% % % xlim(phi_xlim);
% % % ylim(depth_lim_all);grid MINOR; grid on; set(gca,'YDir','reverse','FontSize',fset)       
% % % subplot(1,4,4);
% % % imagesc(dswtlog,dd,gswt_ppd);
% % % colorbar;colormap cool;hold on 
% % % plot(swt,DEPTH,'b','LineWidth',lw)
% % % hold on
% % % plot(swt_mean,DEPTH,'r','LineWidth',lw)
% % % xlim(sw_xlim);
% % % title('Total water saturation(v/v)','FontSize',fst); ylim(depth_lim_all);grid MINOR; grid on; set(gca,'YDir','reverse','FontSize',fset) 

%%%%%%%%%%%%%%%%% behavior of local errors, local models, local temp, and
%%%%%%%%%%%%%%%%% ratio of error (etrial/emod) and derivatives of errors
%%%%%%%%%%%%%%%%% wrt model;all for one layers (jp=24) and for one model
%%%%%%%%%%%%%%%%% parameters (here intergranular porosity)
% % % if reannealing~=0
% % %     figure(24)
% % %     subplot(2,3,1)
% % %     plot(deplot1,'b','LineWidth',lw)
% % %     xlim([0 jtemp-1]);xlabel('Iteration number','FontSize',fs);ylabel('Plot','FontSize',fs);title('Error(new)-Error(old)','FontSize',fst); grid MINOR; grid on ; set(gca,'FontSize',fset)    
% % %     subplot(2,3,2)
% % %     plot(dmplot1,'b','LineWidth',lw)
% % %     xlim([0 jtemp-1]);xlabel('Iteration number','FontSize',fs);ylabel('Plot','FontSize',fs);title('Model(new)-Model(old)','FontSize',fst); grid MINOR ;grid on ;set(gca,'FontSize',fset) 
% % %     subplot(2,3,3)
% % %     plot(splot1,'b','LineWidth',lw)
% % %     xlim([0 jtemp-1]);xlabel('Iteration number','FontSize',fs);ylabel('Plot','FontSize',fs);title('(Error(new)-Error(old))/(Model(new)-Model(old))','FontSize',fst);grid MINOR ;grid on; set(gca,'FontSize',fset) 
% % %     subplot(2,3,4)
% % %     plot(tplot1,'b','LineWidth',lw)
% % %     xlim([0 jtemp-1]);xlabel('Iteration number','FontSize',fs);ylabel('Plot','FontSize',fs);title('Local Temperature','FontSize',fst);grid MINOR;grid on;set(gca,'FontSize',fset)     
% % %     subplot(2,3,5)
% % %     plot(rplot1,'b','LineWidth',lw)
% % %     xlim([0 jtemp-1]);xlabel('Iteration number','FontSize',fs);ylabel('Plot','FontSize',fs);title('Error(new)/Error(old)','FontSize',fst);grid MINOR;grid on;set(gca,'FontSize',fset)    
% % %     subplot(2,3,6)
% % %     plot(emod_local_plt(:,:,24,1),'b','LineWidth',lw)
% % %     xlim([0 jtemp-1]);xlabel('Iteration number','FontSize',fs);ylabel('Plot','FontSize',fs);title('Local Error','FontSize',fst);grid MINOR;grid on;set(gca,'FontSize',fset)     
% % % end


%%%%%%%%%%%%%%%%%%%%%%%QC 
nhist=40;
figure(25)
subplot(2,4,1)
hist(phii,nhist)
xlim([0 0.22]);
ylabel('Frequency','FontSize',fst)
title('Intergranular porosity(v/v)','FontSize',fst)
grid MINOR
grid on
set(gca,'FontSize',fset)
subplot(2,4,2)
hist(phim,nhist)
xlim([0 0.22]);
ylim([0 20]);
% % % ylabel('Frequency','FontSize',fst)
title('Micro-porosity(v/v)','FontSize',fst)
grid MINOR
grid on
set(gca,'FontSize',fset)
subplot(2,4,3)
hist(phit,nhist)
xlim([0 0.22]);
ylim([0 15]);
% % % ylabel('Frequency','FontSize',fst)
title('Total porosity(v/v)','FontSize',fst)
grid MINOR
grid on
set(gca,'FontSize',fset)
subplot(2,4,4)
hist(swt,nhist)
% % % ylabel('Frequency','FontSize',fst)
title('Water saturation(v/v)','FontSize',fst)
grid MINOR
grid on
set(gca,'FontSize',fset)

subplot(2,4,5)
hist(phii_mean,nhist)
xlim([0 0.22]);
ylim([0 12]);
ylabel('Frequency','FontSize',fst)
title('Inverted intergranular porosity(v/v)','FontSize',fst)
grid MINOR
grid on
set(gca,'FontSize',fset)
subplot(2,4,6)
hist(phim_mean,nhist)
xlim([0 0.22]);
% % % ylabel('Frequency','FontSize',fst)
title('Inverted Micro-porosity(v/v)','FontSize',fst)
grid MINOR
grid on
set(gca,'FontSize',fset)
subplot(2,4,7)
hist(phit_mean,nhist)
xlim([0 0.22]);
% % % ylabel('Frequency','FontSize',fst)
title('Inverted total porosity(v/v)','FontSize',fst)
grid MINOR
grid on
set(gca,'FontSize',fset)
subplot(2,4,8)
hist(swt_mean,nhist)
ylim([0 40]);
% % % ylabel('Frequency','FontSize',fst)
title('Inverted water saturation(v/v)','FontSize',fst)
grid MINOR
grid on
set(gca,'FontSize',fset)

%%%%%%%%%%%%%%%%%%%%%%%QC 
nhist=40;
figure(26)
subplot(2,4,1)
hist(resist,nhist)
ylabel('Frequency','FontSize',fst)
title('Resistivity (ohm-m)','FontSize',fst)
grid MINOR
grid on
set(gca,'FontSize',fset)
subplot(2,4,2)
hist(den,nhist)
xlim([2.2 2.6]);
% % ylabel('Frequency','FontSize',fst)
title('Density (gr/cm3)','FontSize',fst)
grid MINOR
grid on
set(gca,'FontSize',fset)
subplot(2,4,3)
hist(dtco,nhist)
% % ylabel('Frequency','FontSize',fst)
title('DTCO(microsec/m)','FontSize',fst)
grid MINOR
grid on
set(gca,'FontSize',fset)
subplot(2,4,4)
hist(dtsm,nhist)
% % % ylabel('Frequency','FontSize',fst)
title('DTSM (microsec/m)','FontSize',fst)
grid MINOR
grid on
set(gca,'FontSize',fset)
subplot(2,4,5)
hist(resist_sim_mean,nhist)
ylabel('Frequency','FontSize',fst)
title('Inverted Resistivity (ohm-m)','FontSize',fst)
grid MINOR
grid on
set(gca,'FontSize',fset)
subplot(2,4,6)
hist(den_sim_mean,nhist)
ylim([0 15]);
% % % ylabel('Frequency','FontSize',fst)
title('Inverted Density (gr/cm3)','FontSize',fst)
grid MINOR
grid on
set(gca,'FontSize',fset)
subplot(2,4,7)
hist(dtp_sim_mean,nhist)
ylim([0 15]);
% % % ylabel('Frequency','FontSize',fst)
title('Inverted DTCO(microsec/m)','FontSize',fst)
grid MINOR
grid on
set(gca,'FontSize',fset)
subplot(2,4,8)
hist(dts_sim_mean,nhist)
ylim([0 12]);
% % % ylabel('Frequency','FontSize',fst)
title('Inverted DTSM (microsec/m)','FontSize',fst)
grid MINOR
grid on
set(gca,'FontSize',fset)


if nruns>1
    nrunny=1:nruns;
    figure(27)
    
        subplot(1,7,1)
        plot(phii_opt(runy,:),DEPTH,'c','LineWidth',lw)
        hold on
        plot(phii_mean,DEPTH,'r','LineWidth',lw)
        hold on
        plot(phii,DEPTH,'b','LineWidth',lw)
% %         legend('Inverted Realizations','Inverted Mean','True')
        title('Intergranular porosity(v/v)','FontSize',fst)
        ylabel('Relative Depth(m)','FontSize',fsy)
        xlim(phi_xlim); ylim(depth_lim_all);grid MINOR; grid on; set(gca,'YDir','reverse','FontSize',fset)       
        subplot(1,7,2)
        plot(phim_opt(runy,:),DEPTH,'c','LineWidth',lw)
        hold on
        plot(phim_mean,DEPTH,'r',phim,DEPTH,'b','LineWidth',lw)
        title('Micro-porosity(v/v)','FontSize',fst)
        xlim(phi_xlim);ylim(depth_lim_all);grid MINOR ;grid on;set(gca,'YDir','reverse','FontSize',fset)       
        subplot(1,7,3)
        plot(phiiw_opt(runy,:)+phimw_opt(runy,:),DEPTH,'c',phit_opt(runy,:),DEPTH,'y','LineWidth',lw)
        hold on
        plot(phimw_mean+phiiw_mean,DEPTH,'r',phimw+phiiw,DEPTH,'b','LineWidth',lw)
        hold on
        plot(phit_mean,DEPTH,':k',phit,DEPTH,'k','LineWidth',lw)
        title('Porosity(v/v)','FontSize',fst)
        xlim(phi_xlim);ylim(depth_lim_all);grid MINOR;grid on;set(gca,'YDir','reverse','FontSize',fset)         
    subplot(1,7,4)
    semilogx(resist_sim_ppd(:,nrunny),DEPTH,'c','LineWidth',lw)
    hold on
    semilogx(resist,DEPTH,'b','LineWidth',lw)
    hold on
    semilogx(resist_sim_mean,DEPTH,'r','LineWidth',lw)
    ylim(depth_lim_all);
    xlim(resist_xlim);
    title('Resistivity(ohm-m)','FontSize',fst);grid MINOR;grid on;set(gca,'YDir','reverse','FontSize',fset)
    subplot(1,7,6)
    plot(dtp_sim_ppd(:,nrunny),DEPTH,'c','LineWidth',lw)
    hold on
    plot(dtco,DEPTH,'b','LineWidth',lw)
    hold on
    plot(dtp_sim_mean,DEPTH,'r','LineWidth',lw)
    title('DTCO(microsec/m)','FontSize',fst);
    xlim(dtp_xlim);ylim(depth_lim_all);
    grid MINOR;grid on;set(gca,'YDir','reverse','FontSize',fset)
    subplot(1,7,7)
    plot(dts_sim_ppd(:,nrunny),DEPTH,'c','LineWidth',lw)
    hold on
    plot(dtsm,DEPTH,'b','LineWidth',lw)
    hold on
    plot(dts_sim_mean,DEPTH,'r','LineWidth',lw)
    title('DTSM(microsec/m)','FontSize',fst);
    xlim(dts_xlim);ylim(depth_lim_all);
    grid MINOR;grid on;set(gca,'YDir','reverse','FontSize',fset)
    subplot(1,7,5)
    plot(den_sim_ppd(:,nrunny),DEPTH,'c','LineWidth',lw)
    hold on
    plot(den,DEPTH,'b','LineWidth',lw)
    hold on
    plot(den_sim_mean,DEPTH,'r','LineWidth',lw)
    title('Density(gr/cm3)','FontSize',fst);
    xlim(den_xlim);ylim(depth_lim_all);
    grid MINOR;grid on;set(gca,'YDir','reverse','FontSize',fset)
end


