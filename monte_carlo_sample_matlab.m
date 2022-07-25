function sample=monte_carlo_sample_matlab(xmin,xmax)
% % % % %  Use the rand function to draw the values from a uniform distribution in the open interval, (xmin,xmax)
ntry=1;
sample=(xmax-xmin)*rand+xmin;
while (sample < xmin || sample >xmax) 
    ntry=ntry+1;
    sample=(xmax-xmin)*rand+xmin;
    if(ntry>100) 
    break
    end
end
end
    





