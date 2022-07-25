
function vec_smth=smth3_handy(vec)

nt=max(length(vec));

vec_smth(1)=(vec(1)+vec(2))/2;
vec_smth(2)=(vec(2)+vec(3))/2;
for jp=3:nt-2
    vec_smth(jp)=(vec(jp)+vec(jp-1)+vec(jp+1))/3;
end
vec_smth(nt-1)=(vec(nt-1)+vec(nt-2))/2;
vec_smth(nt)=(vec(nt)+vec(nt-1))/2;

end