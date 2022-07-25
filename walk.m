function walk=walk(xmod,xmin,xmax,tmp)

ntry=1;
arand=rand;
ayy=0.0;
dif=arand-0.5;
if(dif<0.0) 
    ayy=-1.0;
end
if(dif>=0.0) 
    ayy=1.0;
end
pwr=abs(2.*arand-1.0);
if tmp<=eps
    yy=0.0;
else
    yy=ayy.*tmp.*((1+1./tmp).^pwr-1.0);
end
xmod1=xmod+yy.*(xmax-xmin);

while (xmod1<xmin||xmod1>xmax)
   
    ntry=ntry+1;
    arand=rand;
    ayy=0.0;
    dif=arand-0.5;
    if(dif<0.0) 
        ayy=-1.0;
    end
    if(dif>=0.0) 
        ayy=1.0;
    end
    pwr=abs(2.*arand-1);
    if tmp<=eps
        yy=0.0;
    else
        yy=ayy.*tmp.*((1+1./tmp).^pwr-1.0);
    end
    xmod1=xmod+yy.*(xmax-xmin);

    if(ntry>100) 
       break
    end

end

walk=xmod1;

