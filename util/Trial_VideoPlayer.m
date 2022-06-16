classdef Trial_VideoPlayer < matlab.apps.AppBase
    properties (Access = public)
        VideoPlayerUIFigure             matlab.ui.Figure
        GridLayout                      matlab.ui.container.GridLayout
        Panel                           matlab.ui.container.Panel       
        videoObjects
        VideoUIAxesArray
        PlayPauseButtonsArray
        trialsFrame_strEnd
        FramesCell         % store frames of each trial
    end
    
    properties (Access = private)
        VideoNum_Row = 4;
        VideoUIAxes_width = 300;
        VideoUIAxes_height = 250;
        videoUIAxes_x0 = 30;
        videoUIAxes_y0 = 30;
        videoBtn_width = 100;
        videoBtn_height = 20;
    end
    
    methods (Access = private)
        
         % Create UIFigure and components
        function createComponents(app, inVideoFile, trialTimeXlsFile)
            
            % Create VideoPlayerUIFigure and hide until all components are created
            app.VideoPlayerUIFigure = uifigure('Visible', 'on');
            app.VideoPlayerUIFigure.Position = [100 100 800 350];
            app.VideoPlayerUIFigure.Name = 'trial Video Player';
            
            
            % Create GridLayout
            app.GridLayout = uigridlayout(app.VideoPlayerUIFigure);
            app.GridLayout.ColumnWidth = {'1x'};
            app.GridLayout.RowHeight = {300, '1x'};
            app.GridLayout.ColumnSpacing = 0;
            app.GridLayout.RowSpacing = 0;
            app.GridLayout.Padding = [0 0 0 0];
            app.GridLayout.Scrollable = 'on';

            
            % Create UpperPanel
            app.Panel = uipanel(app.GridLayout);
            app.Panel.Layout.Row = 1;
            app.Panel.Layout.Column = 1;
            app.Panel.Scrollable = 'on';
            
            
            % extract date and bk
            tmp = regexp(inVideoFile, '[0-9]{8}', 'match');
            dateofexp_str = tmp{1};
            tmp = regexp(inVideoFile, '_[0-9]*-', 'match');
            tmp = tmp{1};
            bk = str2num(tmp(2:end-1));
            trialT_strEnd = ExtractTrial_tStrEnd(trialTimeXlsFile, dateofexp_str, bk);
            videoObject = VideoReader(inVideoFile);
            app.trialsFrame_strEnd = trialT_strEnd * videoObject.FrameRate;
            clear tmp dateofexp_str bk
            
            
            % Create VideoUIAxes, Button for each trial   
            ntrials = size(app.trialsFrame_strEnd, 1);
            if ntrials > 2
                nshowTrials = 2;
            else
                nshowTrials = ntrials;
            end
            for tri = 1: nshowTrials
                coli = mod(tri, app.VideoNum_Row);
                rowi =  ceil(tri / app.VideoNum_Row);
                x_axes = app.videoUIAxes_x0 + (coli -1)*app.VideoUIAxes_width;
                y_axes = app.videoUIAxes_y0 + (rowi -1)*app.VideoUIAxes_height;
                x_btn = x_axes + app.VideoUIAxes_width /2 - app.videoBtn_width/2;
                y_btn = y_axes - 20;
                
                
                % create Axes
                videoUIAxes = uiaxes(app.Panel);
                title(videoUIAxes, ['trial='  num2str(tri)]);
                videoUIAxes.Position = [x_axes y_axes app.VideoUIAxes_width app.VideoUIAxes_height];
                videoUIAxes.Tag = ['hAxes-' num2str(tri)];
                app.VideoUIAxesArray(tri) = videoUIAxes;
                
                % create button
                videoButton = uibutton(app.Panel, 'push');
                videoButton.ButtonPushedFcn = createCallbackFcn(app, @VideoPlayPauseButtonPushed, true);
                videoButton.Position = [x_btn y_btn app.videoBtn_width app.videoBtn_height];
                videoButton.Text = 'Play';
                videoButton.Tag = ['hBtn-' num2str(tri)];
                app.PlayPauseButtonsArray(tri) = videoButton;
                
                
                % read frames 
                frames = read(videoObject,[app.trialsFrame_strEnd(tri,1) app.trialsFrame_strEnd(tri,2)]);
                app.FramesCell{tri} = frames;
                
                % show the first frame
                frame1 = frames(:,:,:,1);
                imshow(frame1, 'Parent', videoUIAxes);
                drawnow;
                clear frame1
                
                
                clear coli rowi x_axes y_axes x_btn y_btn
                clear videoUIAxes videoButton frames
            end
            
         
             % Show the figure after all components are created
             app.VideoPlayerUIFigure.Visible = 'on';

            
        end
    end
    
    
    methods (Access = private)
        function PlayingFrames(app, tri)
            frames = app.FramesCell{tri};
            videoUIAxes = findobj(app.Panel, 'Tag', ['hAxes-' num2str(tri)]);
            for framei = 1 : size(frames, 4)
                frame = frames(:, :, :, framei);
                image(frame, 'Parent', videoUIAxes);
                drawnow;
            end
            clear frames
        end
        
        
    end
    
    
    % Callbacks that handle component events
    methods (Access = private)
        
        % Button pushed function: VideoPlayButton
        function VideoPlayPauseButtonPushed(app, event)
            tag = event.Source.Tag;
            tmp = regexp(tag, '[0-9]*', 'match');
            tri = str2num(tmp{1});
            switch event.Source.Text
                case 'Play'
                    event.Source.Text = 'Pause';
                    PlayingFrames(app, tri);
                case 'Pause'
                     event.Source.Text = 'Play';
                otherwise
                     event.Source.Text = 'Wrong';
            end
        end
    end
    
    
    methods (Access = public)
        
        % Construct app
        function app = Trial_VideoPlayer(inVideoFile, trialTimeXlsFile)
            
            % Create UIFigure and components
            createComponents(app, inVideoFile, trialTimeXlsFile);
            
            if nargout == 0
                clear app
            end
        end
        
        % Code that executes before app deletion
        function delete(app)
            
            delete(app.VideoPlayerUIFigure);
        end
    end
end