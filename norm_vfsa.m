function norma=norm_vfsa(A,c)
%%% function to calculate the norm of matrices elementwise, so it is
%%% totally different from the normal definition of matrix norm
%% A matrix of form (n,m)
%% c is scalar or norm degree
%% Sometimes for L3, L5, and L(2k+1) and when A(i,j) are negative, the norm
%% ouput will be complex, so be carful by QC, to make it work for now, I
%% did absolute value for odd norms greater than one 

if c==0
    sprintf('Norm zero is undefined >>>>> divided by zero')
    return;
end
   
if c==1 %%% Norm L1
    B=abs(A);
    norma=(sum(sum(B))).^(1./c);       
else
%     B=A.^c; %%Norm Lc
    B=abs(A.^c); %%Norm Lc   absolute value for odd norms greater than one 
    norma=(sum(sum(B))).^(1./c);  
%     norma=(sum(sum(B)));   
end

% norma=real(norma);
        