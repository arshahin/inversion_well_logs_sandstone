function [veclogout,loge]=edge_log_intepol_noise(veclogin,Shale_log_value,nz,nsh,ksh,pwn_add)

veclogin(1:nsh)=Shale_log_value;
veclogin((nz-nsh):nz)=Shale_log_value;
%%%%%%%%% Linear interpolation at edge 
slop_veclogin1=(veclogin(nsh+1)-veclogin(ksh))/(nsh+1-ksh);
intersect_veclogin1=veclogin(ksh)-slop_veclogin1*ksh;
for jz=ksh:nsh    
    veclogin(jz)=jz*slop_veclogin1+intersect_veclogin1;   
end

slop_veclogin2=(veclogin(nz-nsh-1)-veclogin(nz-ksh))/(nz-nsh-1-nz+ksh);
intersect_veclogin2=veclogin(nz-ksh)-slop_veclogin2*(nz-ksh);

for jz=(nz-nsh):(nz-ksh)    
    veclogin(jz)=jz*slop_veclogin2+intersect_veclogin2;   
end
%%%%%%%%%%%%% adding noise to top and bottom parts of log 
veclogin1=veclogin(1:nsh);
veclogin2=veclogin((nsh+1):(nz-nsh));
veclogin3=veclogin((nz-nsh+1):nz);

noisy1=mean(veclogin1).*pwn_add*randn(length(veclogin1),1);  
noisy3=mean(veclogin3).*pwn_add*randn(length(veclogin3),1);

veclogin1_noisy=veclogin1+noisy1;      
veclogin2_noisy=veclogin2; 
veclogin3_noisy=veclogin3+noisy3; 

veclogout=[veclogin1_noisy;veclogin2_noisy;veclogin3_noisy];

loge=1;

end 