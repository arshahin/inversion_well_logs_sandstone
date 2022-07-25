
function vec_smth=smth_nth_handy(vec,optl)


%%%% optl has to be odd 
if mod(optl,2)==0.0
    sprintf('Operator lenght must be an odd number')
    vec_smth=-999;
    return
end 

nt=length(vec);

if optl > nt
    sprintf('Operator lenght must be smaller than input vector')
    vec_smth=-999;
    return
end 


nc=(optl+1)/2;
vec_smth=zeros(nt,1);

for jp=1:nc-1
    vec_smth(jp)=sum(vec(1:optl))/optl;
end
    
for jp=nc:(nt-nc+1)
    vec_smth(jp)=sum(vec((jp-nc+1):(jp+nc-1)))/optl;
end

for jp=(nt-nc+2):nt
    vec_smth(jp)=sum(vec((nt-optl+1):nt))/optl;
end

end


