function listvideo(para,NN,flux,Temperature,seed,target,num)
span = target/num;
formatSpec = '%.3f';

for i = 1:num
    cov_i = span*(i-1);
    cov_f = span*i;
    Cov_i = num2str(cov_i,formatSpec);
    Cov_f = num2str(cov_f,formatSpec);    
    plot3dvideo(para,NN,flux,Temperature,seed,Cov_i,Cov_f);
end

end