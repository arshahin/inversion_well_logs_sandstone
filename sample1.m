function sample1=sample1(xmod,dx,nx,xmin,xmax)
ix=2.*(rand-0.5).*nx;
sample1=xmod+ix.*dx;
while (sample1 < xmin || sample1 >xmax) ;
    ix=2.*(rand-0.5).*nx;
    sample1=xmod+ix.*dx;
end
    







