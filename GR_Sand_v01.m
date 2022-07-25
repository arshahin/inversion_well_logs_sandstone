function [GR,logr]=GR_Sand_v01(phi_sh,rhow,den,GR_clay,csh,rho_clay,rho_ma,GRma,phii)

cma=1-csh-phii;  %%%% concentration of matrix mineral constituents 

rho_sh=phi_sh.*rhow+(1-phi_sh).*rho_clay; %%% shale density


GR=(csh.*rho_sh.*GR_clay+cma.*rho_ma.*GRma)./den;
logr=1;
end