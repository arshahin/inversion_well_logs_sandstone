function mw=meanw(A,W)

%%% Function to calculate the weigthed mean of a vector,matrix, and cube
%%% A is a vector,matrix, or cube containig the values to be avaraged
%%% W is a vector,matrix, or cube containig the weighting coefficients
%%% mw is a scalar accounting for weighted mean
%%% Note that it won't work for A that has higher dimension than cube like A(N,M,K), but this code can be generalized to handle higher dimensions 

s=size(A);
l=length(s);

if l==3;
    WW=sum(sum(sum(W)));
    if WW==0; sprintf('Summation of weigths is zero >>>>> divided by zero') 
        return;
    else
    mw=sum(sum(sum(A.*W)))./WW; %% Three dimensional cubes like A(N,M,K)
    end
elseif l==2;
    n=min(s);
    if n==1;
        WW=sum(W);
        if WW==0; sprintf('Summation of weigths is zero >>>>> divided by zero')
            return;
        else
            mw=sum(A.*W)./WW; %% vector like A(N,1) or A(1,N)
        end
    else
        WW=sum(sum(W));
        if WW==0; sprintf('Summation of weigths is zero >>>>> divided by zero')
            return;
        else
            mw=sum(sum(A.*W))./WW; %% matrix like A(N,M)
        end
    end
end
        
    
 