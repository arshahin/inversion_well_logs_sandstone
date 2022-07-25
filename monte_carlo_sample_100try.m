function sample=monte_carlo_sample_500try(xmin,dx,nx,xmax)
ntry=1;
ix=2.*(rand-0.5).*nx;
sample=xmin+ix.*dx;
while (sample < xmin || sample >xmax) 
    ntry=ntry+1;
    ix=2.*(rand-0.5).*nx;
    sample=xmin+ix.*dx;
    if(ntry>500) 
    break
    end
end
end
    







