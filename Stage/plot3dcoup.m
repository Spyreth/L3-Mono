function plot3dcoup(para,NN,flux,T,seed,cov,direction,position)

NNd = NN*NN ; 
% plot3dtest(80,600,0.00278,300,1,'0.100','i',200) % où para = 80, N=600,
% '0.100' = coverage, une coupe de i=200
% coverage à mettre avec 2 chiffres après la virgule

stra=append('./Donnees_G/Para',num2str(para),'/a_N_',num2str(NN),'_flux_',num2str(flux),'_T_',num2str(T),'_seed_',num2str(seed),'_COV_',num2str(cov),'.dat');
strb=append('./Donnees_G/Para',num2str(para),'/b_N_',num2str(NN),'_flux_',num2str(flux),'_T_',num2str(T),'_seed_',num2str(seed),'_COV_',num2str(cov),'.dat');

fidua=fopen(stra,'r');
A=fread(fidua,'int32');
tail=size(A);
%Tmax=floor(tail(1)/NN^2)
%uua=reshape(A,NN,NN,Tmax);
uuA=reshape(A,NN,NN);
uua=transpose(uuA) ;
whos uua

fidub=fopen(strb,'r');
B=fread(fidub,'int32');
%uub=reshape(B,NN,NN,Tmax);
uuB=reshape(B,NN,NN);
uub=transpose(uuB) ;
%max(max(uua))
%max(max(uub)) 

a1x = 1; 
a1y = 0; 
a2x = 1/2; 
a2y = sqrt(3)/2;
dax = -1/4;
day = -1/(4*sqrt(3)); 
dbx = 1/4;
dby = 1/(4*sqrt(3)); 

vrax = zeros(NN,NN);
vray = zeros(NN,NN);
vraz = zeros(NN,NN);
vrbx = zeros(NN,NN);
vrby = zeros(NN,NN);
vrbz = zeros(NN,NN);
Vrx = 1:2*NN*NN; 
Vry = 1:2*NN*NN; 
Vrz = 1:2*NN*NN; 


if direction == 'i' %car 'reshape' de Matlab donne une matrice transposée 
    k=1;
    for j=1:1:NN          
        coup(k)= uua(position,j);
        coup(k+1)= uub(position,j);
        k = k + 2;
    end
else
    k=1; 
    for i=1:1:NN                    
        coup(k)= uua(i,position);
        coup(k+1)= uub(i,position);
        k = k + 2;
    end
 end

 %
 for i = 1:1:NN;
    for j=1:1:NN;
       vrax(i,j) = (i-1)*a1x + (j-1)*a2x + dax ; 
       vray(i,j) = (i-1)*a1y + (j-1)*a2y + day ;
       vraz(i,j) = uua(i,j); 
       %
       vrbx(i,j) = (i-1)*a1x + (j-1)*a2x + dbx ; 
       vrby(i,j) = (i-1)*a1y + (j-1)*a2y + dby ;
       vrbz(i,j) = uub(i,j);
    end
  end

zmax = max(max(max(uua))); 

%S = repmat([2],numel(Vrx),1);
%C = repmat([1 0 1],numel(Vrx),1); % C marker color 


for k = 1:1:NN*NN
    j = floor((k-1)/NN) + 1 ; 
    i = k - (j-1)*NN; 
    Vrx(k) = vrax(i,j);
    Vry(k) = vray(i,j);
    Vrz(k) = vraz(i,j);
    Vrx(k+NNd) = vrbx(i,j);
    Vry(k+NNd) = vrby(i,j);
    Vrz(k+NNd) = vrbz(i,j);
end

%tiledlayout(2,1)
%nexttile
title(['Para:',num2str(para),' , \theta =',num2str(str2double(cov)),'ML'])
%pl = scatter3(Vrx,Vry,Vrz,S,C,'filled');
pl = scatter3(Vrx,Vry,Vrz,[2],Vrz,'filled');
axis([0 3*NN/2 0 3*NN/2 0 3*zmax])
view([0  0 1])
map1 = [1 0 0
    1 0.5 0
    1 0.8 0
    1 1 0];
%colormap(map1)
%colorbar('Ticks',[0.5,1.5,2.5,3.5], ...
%    'TickLabels',{'Sub','h=1','h=2','h>=3'})                        %// this is optional, just to make sure the colorbar does not vary

%clim([0 4])  ;

colormap(autumn(5))

plstruct = get(pl);
legend(plstruct.Children, append('Para',num2str(para),', \theta = ',num2str(cov),'ML'));
%w = waitforbuttonpress;
hold on
x = linspace(0,3/2*NN, 10);
if direction=='i' %car 'reshape' de Matlab donne une matrice transposée 
    y=sqrt(3)*(x-position) ;
else 
    y=0*x+sqrt(3)*position/2;
end
z=0*x+2*zmax ;
plot3(x,y,z,'LineWidth',1,'color','black')
legend(append(direction,' = ', num2str(position)));

figure
coupx = 0.5:0.5:NN;
plot(coupx,coup)
set(gca,'YLim',[-1 10]) %to be fitted
title('Topographic')

%nexttile
%C=[A,B];
%%histogram(A,Normalization="percentage");
%hc = histcounts(C);
%b = bar(hc);
%% percent of total for each bar
%s = compose('%.3f%%', hc / sum(hc) * 100);
%yOffset = 100; % tweat, as necessary
%text(b.XData, b.YEndPoints + yOffset,s);

para 

fclose(fidua);
fclose(fidub);






