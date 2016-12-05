clear all; close all;

workingDir = 'D:/Plot';
imageNames = dir(fullfile(workingDir,'*.png'));
imageNames = {imageNames.name}';
vect = zeros(1,length(imageNames));
for ii = 1:length(imageNames)
   name = imageNames{ii};
   new_name = strrep(name, 'plot', '');
   new_name = strrep(new_name, '.png', '');
   vect(ii) = str2double(new_name);
end
[vect, idx] = sort(vect);
sorted_imageNames = imageNames;
for ii = 1:length(imageNames)
   sorted_imageNames{ii} = imageNames{idx(ii)};
end
imageNames = sorted_imageNames';

outputVideo = VideoWriter(fullfile('shuttle_out.avi'));
outputVideo.FrameRate = 60;
open(outputVideo)

Params = genConf();

for ii = 1:length(imageNames)
   img = imread(fullfile(workingDir,imageNames{ii}));
   p = Params(ii,:);
   text = ['Right Mirror Angle: ',num2str(p(1)), ' Camera position: ', num2str(p(2)), ' Head Position: (',num2str(p(3)),',',num2str(p(4)),')'];
   img = insertText(img,[10,10],text);
   writeVideo(outputVideo,img)
end

close(outputVideo)