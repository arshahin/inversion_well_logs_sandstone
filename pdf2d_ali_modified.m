function [pdf2d_ali,cp1,cp2]=pdf2d_ali_modified(x,xmin,xmax,nx,y,ymin,ymax,ny)

binx=abs(xmax-xmin)./nx;
biny=abs(ymax-ymin)./ny;

n=size(x,1);
m=size(y,2);

if n<=m; x=x';
end

n=size(y,1);
m=size(y,2);

if n<=m; y=y';
end

d=[x,y];
N=size(d,1);
A=zeros(ny,nx);
for k=1:nx;
    for m=1:ny;
        for i=1:N;
            if ((d(i,1)<=(xmin+k*binx)) & (d(i,1)>(xmin+(k-1)*binx)) & (d(i,2)<=(ymin+m*biny)) & (d(i,2)>(ymin+(m-1)*biny)))
                A(m,k)=A(m,k)+1;
            end
        end
    end
end

cp1=sum(sum(A));%%% Check point cp1 should be N_ctotal number of pairs
pdf2d_ali=A./N;%% Joint density function
cp2=sum(sum(pdf2d_ali));%%% Check point cp2 should be 1_cumulative probability
%%%% smoothing
%optbin=floor(sqrt(length(d)/0.03));
%smfilt=fspecial('gaussian',floor(optbin/3),floor(sqrt(optbin)/2));
smfilt=fspecial('average');
pdf2d_ali=filter2(smfilt,pdf2d_ali);
% figure(5)
% p_label=linspace(pmin,pmax,biny);
% k_label=linspace(kmin,kmax,binx);
% surfc(p_label,k_label,Jpdf);
% xlabel('Porosity')
% ylabel('Permeability(mD)')
% title('Joint Probability Distribution of Permeability and Porosity')

% x_label=linspace(xmin,xmax,binx);
% y_label=linspace(ymax,ymin,biny);
% imagesc(x_label,y_label,pdf2d_ali);
% colorbar
% xlabel('X')
% ylabel('Y')
% title('Joint Probability Distribution')
% 

                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                