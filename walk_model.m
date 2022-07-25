function [x_new,pdf,cdf]=walk_model(x_old,dx,xmin,xmax,tmp)

%example
% clc
% close all
% clear all
% xmin=-7.0;
% xmax=3.0;
% x_old=-10;
% dx=0.05; 
% tmp=0.001;

% set the grid and discritizing the interval
% dx=(xmax-xmin)./(nbins-1); %% ======> xmax=xmin+(nbins-1).*dx  , so x(k)=xmin+(k-1).*dx    while x(1)=xmin and x(nbins)=xmax
% scale=1.0;
if (x_old<xmin||x_old>xmax); sprintf('Error: x_old is outside of interval')
    return;
end; 

nbins=(xmax-xmin)./dx+1;
aa=0:1:nbins-1;
x=xmin+aa.*dx;

%%%%%%%%%%%%%%%% Calculate Cauchy like CDF 
 y=(x-x_old)./(xmax-xmin);

% for i=1:nbins;
%     if x(i)<= x_old; y(i)=(x(i)-x_old)./(x_old-xmin); %%% DISTANCE VECTOR equivalent with y=(x-x_old)./(xmax-xmin); if don't do this way, cdf will not be in [0 1] interval
%     else y(i)=(x(i)-x_old)./(xmax-x_old);
%     end
% end   

c=abs(y)+tmp;
d=log(1+1./(tmp+eps));
cauch=1./(2.*c.*d); %%%% PDF from Ingber 1993 pg. 107 from book Sen/Stoffa
cdf=cumsum(cauch);


ss=max(cdf);

if ss~=1; cdf=cdf./ss; %%% Normalized Cauchy like CDF
    pdf(2:nbins)=cdf(2:nbins)-cdf(1:nbins-1); %%%%%Final Cauchy like PDF
    pdf(1)=cdf(1);
else
    pdf(2:nbins)=cdf(2:nbins)-cdf(1:nbins-1); %%%%%Final Cauchy like PDF
    pdf(1)=cdf(1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% get the new sample
x_rand=rand;
if (x_rand < cdf(1)),ind1=1;
elseif(x_rand > cdf(nbins));ind1=nbins;
else    
    pp=abs(cdf-x_rand);
    qq=min(pp);    
    ind1=find(pp==qq);
    if max(size(ind1))>1; ind1=round(mean(ind1)); end;
end
x_new=x(ind1); % fianl new model parameter	

if (x_new <xmin||x_new >xmax); sprintf('X-new is not in interval'); 
    return; end; 
 
% % 
% figure(1)
% plot(x,cdf);
% 
% figure(2)
% plot(x,pdf);

% figure(3);hist(x_new,50)
