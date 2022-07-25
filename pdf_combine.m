function [x_new,pdf,cdf]=pdf_combine(prior,nx,xmin,xmax,x_old,tmp,w,methodmc)
%%% methodmc==2 is recommended

%%% function to combine prior information with vfsa 
%%% prior is a prior pdf e.g. p-wave velocity or porosity prior data from well logs
%%% Note that sum(prior)==1; but function is written in the way that normalized prior 
%%% nx: number of bins to discritize interval xmin and xmax
%%% x_old is the initial solution or guess
%%% tmp is the tempreture
%%% w is the weight that we assign to prior,so 1-w is the weight assigned to vfsa sampling ,e.g. if w=0.2, it means that 
%%% we calculate sampling probability based on the follwoing equation p=(prior.^0.2).*(pvfsa.^0.8)

%%%% check input data
if (x_old<xmin||x_old>xmax); sprintf('Error: x_old is outside of interval')
    return;
end
if (w<0||w>1); sprintf('Error: Weights shouldn,t be less than zero or greater than one')
    return;
end

if w~=0;
%%%% Prior CDF to be normalized
    cdf_prior=cumsum(prior);
    ss=max(cdf_prior);
    if ss~=1; cdf_prior=cdf_prior./ss; prior(1)=cdf_prior(1); prior(2:nx)=cdf_prior(2:nx)-cdf_prior(1:nx-1); end;
end

if w~=1;
    % set the grid and discritizing the interval
    dx=(xmax-xmin)./(nx-1); %% ======> xmax=xmin+(nx-1).*dx  , so x(k)=xmin+(k-1).*dx    while x(1)=xmin and x(nx)=xmax
    aa=0:1:nx-1;
    x=xmin+aa.*dx;
    %%%%%%%%%%%%%%%% Calculate Cauchy like PDF 
    y=(x-x_old)./(xmax-xmin); %%% DISTANCE VECTOR
    c=abs(y)+tmp;
    d=log(1+1./(tmp+eps));
    cauch=1./(2.*c.*d); %%%% PDF from Ingber 1993 pg. 107 from book Sen/Stoffa
    cvfsa=cumsum(cauch);
    ss=max(cvfsa);

    if ss~=1; cvfsa=cvfsa./ss; %%% Normalized Cauchy like CDF
        pvfsa(1)=cvfsa(1);
        pvfsa(2:nx)=cvfsa(2:nx)-cvfsa(1:nx-1); %%%%%Final Cauchy like PDF    
    else
        pvfsa(1)=cvfsa(1);
        pvfsa(2:nx)=cvfsa(2:nx)-cvfsa(1:nx-1); %%%%%Final Cauchy like PDF    
    end
end

if w==0; pdf=pvfsa;
elseif w==1; pdf=prior;
else pdf=(prior.^w).*(pvfsa.^(1-w)); %% Combined PDF 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% compute the cdf  for the combined ppd and then normalize 
cdf=cumsum(pdf);
ss=max(cdf);
if ss~=1; cdf=cdf./ss; pdf(1)=cdf(1);pdf(2:nx)=cdf(2:nx)-cdf(1:nx-1); end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% get the new sample
x_new=mc_sampling(pdf,xmin,xmax,nx,methodmc);


