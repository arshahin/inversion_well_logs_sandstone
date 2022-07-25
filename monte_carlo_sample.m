function sample=monte_carlo_sample(xmin,dx,nx,xmax)
ix=2.*(rand-0.5).*nx;
sample=xmin+ix.*dx;
while (sample < xmin || sample >xmax) 
    ix=2.*(rand-0.5).*nx;
    sample=xmin+ix.*dx;
end
    







