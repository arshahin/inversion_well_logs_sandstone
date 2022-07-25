function B=subgrid_1d(A,s)

%%%% Convert the course grid vector A to fine grid vector B

%%%% s is the scaling factor to convert

%%%% A vector of size (n,1) 
%%%% B vector of size (s*n,1) 

s=round(s);
n=length(A);
for i=1:n;
    B(((i-1)*s+1):(i*s))=A(i);
end




