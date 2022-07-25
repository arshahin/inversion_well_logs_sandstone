function [resist,ffac,rw,logr]=sdem_resist_carbonate(phit,phii,phiv,swi,swv,lamdai,lamdav,tmpc,salppk,wetcase)

%%%%%%%%%%%%%%%%%%%% version v01
%%%%%% code to compute formation factor and resistivity response of one(or several) core
%%%%%% plugs or reservoir depth intervals with dual porosity and dual
%%%%%% saturation 


%%%%% wetcase 1 is water wet
%%%%% wetcase 2 is oil wet 


%%%%%%%%%%%%%%%%%%%% Input
%%%%%% phii,phiv: Intergranular and Vuggyy porosities where phit=phii+phiv
%%%%%% swi and swv are intergranular and vuggy water saturations where
%%%%%% swt.phit=swi.phii+swv.phiv

% % % Compute Rw using Arp's equation
%%%%% ppm_br  ppm of NaCl in water, 
%%%%% tmpc tempreture in celcious(oC)
ppm_br=salppk.*1000;
tmpf=1.8*tmpc+32;%%% Farenhite
rw=(0.0123+(3647.5)./(ppm_br.^0.955)).*((81.77)/(tmpf+6.77));


%%%%%%% make them all non zero 
phii=phii+eps;

phiv=phiv+eps;

phit=phit+eps;

swi=swi+eps;

swv=swv+eps;

%%% water-filled porosities 
phiiw=phii.*swi; 
phivw=phiv.*swv;
phitw=phiiw+phivw;
swt=(phiiw+phivw)./phit;


if wetcase==1  %%%% water wet 
    ffac=(1./phiiw).^lamdai.*(phiiw./phitw).^lamdav;
    resist=rw.*ffac;
    logr=1;
elseif wetcase==2  %%%%% oil wet 
      nv=lamdav;
      ni=lamdai+(lamdai-lamdav).*(phiv./phii);
      ffac1=(1./phii).^lamdai.*(phii./phit).^lamdav;
      ffac2=(1./swi).^ni.*(swi./swt).^nv;  
      ffac=ffac1.*ffac2;  
      resist=rw.*ffac;      
      logr=1;  
    else
        sprintf('The wetcase is not in library')
    logr=-1;
    rw=-999;
    ffac=-999;
    resist=-999;        
    return;
       
end
        
    







