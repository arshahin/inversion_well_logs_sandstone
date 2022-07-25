
function [sample,count]=mc_sampling_reproducable(pdf,xmin,xmax,nx,randf,count,methodmc)

%%% methodmc==2 is recommended

%%%% see example in mc_sampling_test.m
%%%%% Monte Carlo sampling from a general pdf

% set the grid and discritizing the interval
dx=(xmax-xmin)./(nx-1); %% ======> xmax=xmin+(nx-1).*dx  , so x(k)=xmin+(k-1).*dx    while x(1)=xmin and x(nx)=xmax
aa=0:1:nx-1;
x=xmin+aa.*dx;

cdf=cumsum(pdf);
cp=max(cdf);
cdf=cdf./cp;

if methodmc==1; %%%Alireza method    
    cdf_rand=randf(count); %%% get random number 
    count=count+1;
    if (cdf_rand < cdf(1))
        ind1=1;    
    elseif(cdf_rand > cdf(nx))
        ind1=nx;
    else    
        pp=abs(cdf-cdf_rand);
        qq=min(pp);    
        ind1=find(pp==qq);
        if max(size(ind1))>1; ind1=round(mean(ind1)); end;
    end
    sample=x(ind1); % mc random sample
    
elseif methodmc==2; %%% matlab interpolation
    cdf_rand=randf(count); %%% get random number 
    count=count+1;
    sample= interp1(cdf,x,cdf_rand,'spline');
    while (sample < xmin || sample >xmax) ;
        cdf_rand=randf(count); %%% get random number 
        count=count+1;
        sample= interp1(cdf,x,cdf_rand,'spline'); 
    end
else
    sprinf('methodmc must be either 1 or 2')
    return
end
    

if (sample <xmin||sample >xmax); sprintf('Error: sample is not in interval')
    return
end
