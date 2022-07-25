function nv=norm_vfsa_apply(A,B,scale,normd)
a=norm_vfsa((A-scale.*B),normd);
b=norm_vfsa((A+scale.*B),normd);
nv=2.*a./(a+b+eps);