function plot3dvideo(para,NN,flux,Temperature,seed,cov_i,cov_f)

NNd = NN*NN ; 
% plot3dvideo(100,600,0.000278,300,1,'0.000','0.100') % où para = 100,
% N=200,T=300 2.0-3.0 ML
% coverage,Seed 1

stra=append('./Donnees_G/Para',num2str(para),'/Movie_a_N_',num2str(NN),'_flux_',num2str(flux),'_T_',num2str(Temperature),'_seed_',num2str(seed),'_COV_',num2str(cov_i),'-',num2str(cov_f),'.dat');
strb=append('./Donnees_G/Para',num2str(para),'/Movie_b_N_',num2str(NN),'_flux_',num2str(flux),'_T_',num2str(Temperature),'_seed_',num2str(seed),'_COV_',num2str(cov_i),'-',num2str(cov_f),'.dat');

fidua=fopen(stra,'r');
A=fread(fidua,'int32');
tail=size(A);
Tmax=floor(tail(1)/NN^2);
uua=reshape(A,NN,NN,Tmax);

whos uua

fidub=fopen(strb,'r');
B=fread(fidub,'int32');
uub=reshape(B,NN,NN,Tmax);


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

strvideo = append('./Movie/Para',num2str(para),'/growthfilePara',num2str(para),' N = ',num2str(NN),' T=',num2str(Temperature),' seed',num2str(seed),' from ', num2str(cov_i),'ML to ',num2str(cov_f),'ML' );
v = VideoWriter(strvideo,'MPEG-4');
open(v);

fig1 = figure(1);

for tps = 1:1:Tmax
for i = 1:1:NN
    for j=1:1:NN
    vrax(i,j) = (i-1)*a1x + (j-1)*a2x + dax ; 
    vray(i,j) = (i-1)*a1y + (j-1)*a2y + day ;
    vraz(i,j) = uua(j,i,tps);
    %
    vrbx(i,j) = (i-1)*a1x + (j-1)*a2x + dbx ; 
    vrby(i,j) = (i-1)*a1y + (j-1)*a2y + dby ;
    vrbz(i,j) = uub(j,i,tps);
    
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
   %
    %[r,g,b] = color2(Vrz(k)) ;
    % C(k,1) = r ;
    % C(k,2) = g ;
    % C(k,3) = b ;
    %
    % [r,g,b] = color2(Vrz(k+NNd)) ;
    % C(k+NNd,1) = r ; 
    % C(k+NNd,2) = g ; 
    % C(k+NNd,3) = b ; 

end


pl = scatter3(Vrx,Vry,Vrz,[2],Vrz,'filled');

hold off
grid on

videoframe=(str2double(cov_f)-str2double(cov_i))/(Tmax-1); 
title(['Para ',num2str(para)]);
subtitle(['T=',num2str(Temperature),'k, \theta =',num2str(str2double(cov_i)+videoframe*(tps-1)),'ML'])

axis([0 3*NN/2 0 3*NN/2 0 zmax])
view([0 0 3*zmax])
map = [0 0 0
    0 0 1
    0 1 1
    1 1 0];
map1 = [1 0 0
    1 0.5 0
    1 0.8 0
    1 1 0
    0 1 0];
colormap(map1)
%colormap(autumn(5))
colorbar('Ticks',[0.5,1.5,2.5,3.5,4.5], ...
    'TickLabels',{'Sub','h=1','h=2','h=3','h>=4'})                        %// this is optional, just to make sure the colorbar does not vary

clim([0 5])  ;

F=getframe(fig1);
writeVideo(v,F)

end 
para

fclose(fidua);
fclose(fidub);

close(v);
end
%fclose(fiduSi0);
%fclose(fiduAg0);
