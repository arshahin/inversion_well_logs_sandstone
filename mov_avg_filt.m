
function vec_movavg=mov_avg_filt(vec,optlen)
% % % Find the moving average of the data 


vs1=size(vec,1);%%%% no. of row 
vs2=size(vec,2);%%%% no of columns 
if vs2>vs1
    vec=vec';
end

ext=ones(optlen+1,1).*vec(1);
vec1=[ext;vec];
nv1=max(size(vec1));

b=(1/optlen)*ones(1,optlen);
    a=1;
    vec2=filter(b,a,vec1);
    vec_movavg=vec2(optlen+2:nv1);    
    
end


