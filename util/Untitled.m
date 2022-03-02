clear
folder = 'H:\My Drive\NMRC_umn\Projects\FCAnalysis\exp\data\Animals\Kitty\KittyMinMaxPSD';
inVideoFilename = 'Kitty20150408_2-1394 Desktop Video Camera.avi';
inVideoFile = fullfile(folder, inVideoFilename);

trialTimeXlsFile = fullfile(folder, 'KittyTrialTime.xlsx');

videoObject = VideoReader(inVideoFile);
frames = read(videoObject, [1 1000]);
fig = figure;
ax = axes(fig);
for fi = 1 : size(frames, 4)
    frame = frames(:, :, :, fi);
    image(frame,'Parent', ax);
    drawnow;
    if mod(fi, 15) == 1
        tic;
    end
    if mod(fi, 15) == 0
        tdur = toc;
        disp(['fi = ' num2str(fi) 'fs =' num2str(15/tdur)])
    end
end

