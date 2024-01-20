function videoslowdown(para,NN,Temperature,seed,target,time)

formatSpec = '%.3f';

videoFile = append('./Movie/Para',num2str(para),'/MergedVideoPara',num2str(para),' N = ',num2str(NN),' T=',num2str(Temperature),' seed',num2str(seed),' from 0.000ML to ',num2str(target,formatSpec),'ML.mp4'); % ex : ".\Movie\Para1000\growthfilePara1000 N = 300 T=300 seed1 from 0.000ML to 0.100ML.mp4"
videoObj = VideoReader(videoFile);
outputFile =append( './Movie/Para',num2str(para),'/Video(',num2str(time),' times slow) Para',num2str(para), ' NN=',num2str(NN), ' T=',num2str(Temperature), ' seed',num2str(seed),' from ', num2str(0.000),'ML to ',num2str(target,formatSpec),'ML.mp4') ;
outputVideo = VideoWriter(outputFile, 'MPEG-4');
outputVideo.FrameRate = videoObj.FrameRate ; % Change this factor to adjust the speed
open(outputVideo);
while hasFrame(videoObj)
    frame = readFrame(videoObj);
    
    % Repeat each frame to slow down the video
    for i = 1:time  % This repeats each frame (time) times, making it (time)x slower
        writeVideo(outputVideo, frame);
    end
end

close(outputVideo);