
function [pdf1d,meshx]=pdf1d(d,dmin,dmax,nbins)
%%% to speed up the code, put nx at 50 or less
nx=nbins;
% % xmin=min(d);
% % xmax=max(d);
xmin=dmin;
xmax=dmax;

binx=abs(xmax-xmin)/nx;

N=max(size(d));
A=zeros(nx,1);

for k=1:nx
    for i=1:N
        if (d(i)<=(xmin+k*binx)) && ...
           d(i)>(xmin+(k-1)*binx) 
            A(k,1)=A(k,1)+1;
        end
    end
end
pdf1d=A./N;%%  density function
% % % Nr=nsmth;%%% number of samples for smoothing. Put it one for no smoothing 
%%% Nr=3 is usually a good choice 
% % % pdf1d=smooth(pdf1d,Nr);
meshx=linspace(xmin,xmax,nx);
% % plot(meshx,pdf1d)

end


                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                