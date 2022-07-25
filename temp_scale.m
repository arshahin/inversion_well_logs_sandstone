function [tmp_local,smax,R,deplot,dmplot,splot,tplot,rplot]=temp_scale(tmp_local,emod_local,etrial_local,model_mod,model_trial,model_min,model_max,nlayers,nparam,gama,reannealing)

% Rescale tmp_local according to Ingnber 1993 to include sensitivity of different model parameters to 
% update(reannealing) local tempreture based on the derivative of local errors wrt local model 
        

%% emod_local(nlayers,nparam) 
%% tmp_local (nlayers,nparam) 
%% model_mod(nlayers,nparam)
%% model_trial(nlayers,nparam)
%% model_max(nlayers,nparam)
%% model_min(nlayers,nparam)

 %% Initialize plots if not assinged later
R=0; 
deplot=0;
dmplot=0;
splot=0;
tplot=0;
rplot=0;

if reannealing==1 %%%% ubdapte temp based on the Ingber sensitivity 
    smax=-1.e10; %%% big number to initialize the Smax
    for jlay=1:nlayers
        for jpar=1:nparam 
            de=etrial_local(jlay,jpar)-emod_local(jlay,jpar);
            dm=model_trial(jlay,jpar)-model_mod(jlay,jpar);
            if(abs(dm)>0.005) %%% exclude dm==0
                intv=model_max(jlay,jpar)-model_min(jlay,jpar);
                s=abs(intv*de/dm); 
                smax=max(smax,s);
            end
         end
    end  

    for jlay=1:nlayers
        for jpar=1:nparam  
            de=etrial_local(jlay,jpar)-emod_local(jlay,jpar);
            dm=model_trial(jlay,jpar)-model_mod(jlay,jpar);
            if(abs(dm)>0.005)%%% exclude dm==0
                intv=model_max(jlay,jpar)-model_min(jlay,jpar);
                s=abs(intv*de/dm); 
                R=(smax/s)^gama;
                ER=(etrial_local(jlay,jpar)/emod_local(jlay,jpar))^gama; %%% error ratio
                tmp_local(jlay,jpar)=tmp_local(jlay,jpar)*R;
                if tmp_local(jlay,jpar)>=10
                    tmp_local(jlay,jpar)=10;
                end 
                    
                if(jpar==1 && jlay==24)%%% for plotting the fitting parametres of one of the layers (QC on fly), here intergranular porosity of layer 24th
                    deplot=de;
                    dmplot=dm;
                    splot=s;
                    rplot=ER;
                    tplot=tmp_local(jlay,jpar);
                end
            end
        end
    end
elseif reannealing==2 %%%%% Udpate local temps based on the LTR(Armando) 
    smax=0;
     for jlay=1:nlayers
        for jpar=1:nparam
            if(emod_local(jlay,jpar)>0.005) %%% exclude emod_local(jlay,jpar)==0
                R=etrial_local(jlay,jpar)/emod_local(jlay,jpar); %% R=new_error/old_error, when new_error>old_error==> R>1==> increase temp
                R=R^gama;
                tmp_local(jlay,jpar)=tmp_local(jlay,jpar)*R;
                if(jpar==1 && jlay==24)%%% for plotting the fitting parametres of one of the layers (QC on fly), here intergranular porosity of layer 24th
                    de=etrial_local(jlay,jpar)-emod_local(jlay,jpar);
                    dm=model_trial(jlay,jpar)-model_mod(jlay,jpar);
                    intv=model_max(jlay,jpar)-model_min(jlay,jpar);
                    s=abs(intv*de/dm); 
                    deplot=de;
                    dmplot=dm;
                    splot=s;
                    tplot=tmp_local(jlay,jpar);
                    rplot=R;  %%%error ratio 
                end
            end
        end
     end   
else
    sprintf('reannealing must be either 1 or 2')
    return
end

    