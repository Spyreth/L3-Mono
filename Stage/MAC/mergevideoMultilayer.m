function mergevideoMultilayer(para,NN,Temperature,seed,target,num)
%The span of each video
span = target/num;
formatSpec = '%.3f';

%create new video fichier
vfile = append('./Movie/Para',num2str(para),'/MergedVideoPara',num2str(para),' N = ',num2str(NN),' T=',num2str(Temperature),' seed',num2str(seed),' from 0.000ML to ',num2str(target,formatSpec),'ML'); 
v=VideoWriter(vfile,'MPEG-4');
open(v)

% Set up parameters
for i = 1:num
    cov_i = span*(i-1);
    cov_f = span*i;
    Cov_i = num2str(cov_i,formatSpec);
    Cov_f = num2str(cov_f,formatSpec);
    Vidfile = append('./Movie/Para',num2str(para),'/growthfilePara',num2str(para),' N = ',num2str(NN),' T=',num2str(Temperature),' seed',num2str(seed),' from ',Cov_i,'ML to ',Cov_f,'ML.mp4'); % Video i
    Vid=VideoReader(Vidfile);
    % Iterate on all frames in video i and write one frame at a time
    while hasFrame(Vid) 
        Video = readFrame(Vid); % read each frame
        writeVideo(v,Video) % write each frame
    end
end

close(v)