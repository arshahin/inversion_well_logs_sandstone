function [hist2d_stick,meshx,meshy]=hist2d(array,array_min,array_max,nbins,nsmth,msmth)
    d=array;  
    [n,m]=size(d);
    
    ymin=array_min;
    ymax=array_max;
    ny=nbins;
    
    xmin=0;
    xmax=m;
    nx=m;    
    
    hist2d_stick=zeros(nbins,m);    
    for jk=1:m
        [hist2d_stick(:,jk),~]=pdf1d(d(:,jk),ymin,ymax,ny);
% %         hist2d_stick(:,jk)=histogram(d(:,jk),nbins);
% %         [hist2d_stick(:,jk),~] = histcounts(d(:,jk),nbins);
    end  
    Nr=nsmth;%%% smoothing in row direction 
    Nc=msmth; %%% smmothing in column direction 
    hist2d_stick = smooth2a(hist2d_stick,Nr,Nc);
% %     hist2d_stick=fliplr(hist2d_stick);
    [meshx,meshy]=meshgrid(linspace(xmin,xmax,nx),linspace(ymin,ymax,ny));
end
    
    