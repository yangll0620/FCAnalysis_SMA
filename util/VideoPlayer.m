classdef VideoPlayer < matlab.apps.AppBase
    properties (Access = public)
        VideoPlayerUIFigure             matlab.ui.Figure
        GridLayout                     matlab.ui.container.GridLayout
        UpperPanel                      matlab.ui.container.Panel
        VideoUIAxes                    matlab.ui.control.UIAxes % showing video
        BottomPanel                     matlab.ui.container.Panel
        VideoPlayPauseButton            matlab.ui.control.Button
        videoObject
    end
    
    methods (Access = private)
        
         % Create UIFigure and components
        function createComponents(app)
            
            % Create VideoPlayerUIFigure and hide until all components are created
            app.VideoPlayerUIFigure = uifigure('Visible', 'off');
            app.VideoPlayerUIFigure.Position = [100 100 400 350];
            app.VideoPlayerUIFigure.Name = 'Video Player';
            
            
            % Create GridLayout
            app.GridLayout = uigridlayout(app.VideoPlayerUIFigure);
            app.GridLayout.ColumnWidth = {'1x'};
            app.GridLayout.RowHeight = {300, '1x'};
            app.GridLayout.ColumnSpacing = 0;
            app.GridLayout.RowSpacing = 0;
            app.GridLayout.Padding = [0 0 0 0];
            app.GridLayout.Scrollable = 'on';

            
            % Create UpperPanel
            app.UpperPanel = uipanel(app.GridLayout);
            app.UpperPanel.Layout.Row = 1;
            app.UpperPanel.Layout.Column = 1;
            app.UpperPanel.Scrollable = 'on';
            
            
            % Create VideoUIAxes 
            app.VideoUIAxes = uiaxes(app.UpperPanel);
            title(app.VideoUIAxes, 'Video')
            app.VideoUIAxes.Position = [30 36 326 250];
            
            
            % Create RightPanel
            app.BottomPanel = uipanel(app.GridLayout);
            app.BottomPanel.Layout.Row = 2;
            app.BottomPanel.Layout.Column = 1;
            app.BottomPanel.Scrollable = 'on';
            
            
            % Create VideoPlayButton
            app.VideoPlayPauseButton = uibutton(app.BottomPanel, 'push');
            app.VideoPlayPauseButton.ButtonPushedFcn = createCallbackFcn(app, @VideoPlayPauseButtonPushed, true);
            app.VideoPlayPauseButton.Position = [19 71 108 22];
            app.VideoPlayPauseButton.Text = 'Play';
            
            
            % Show the figure after all components are created
            app.VideoPlayerUIFigure.Visible = 'on';
            
            
            input_video_file = 'H:\My Drive\NMRC_umn\Projects\FCAnalysis\exp\data\Animals\Kitty\KittyMinMaxPSD\Kitty20150408_2-1394 Desktop Video Camera.avi';
            app.videoObject = VideoReader(input_video_file);
            
            % Display first frame
            frame_1 = read(app.videoObject,1);
            imshow(frame_1, 'Parent', app.VideoUIAxes);
            drawnow;
            
        end
    end
    
    
    methods (Access = private)
        function PlayingVideo(app)
            for framei = 1 : app.videoObject.NumberOfFrames
                frame = read(app.videoObject,framei);
                imshow(frame, 'Parent', app.VideoUIAxes);
                drawnow;
            end
        end
    end
    
    
    % Callbacks that handle component events
    methods (Access = private)
        
        % Button pushed function: VideoPlayButton
        function VideoPlayPauseButtonPushed(app, event)
            switch app.VideoPlayPauseButton.Text
                case 'Play'
                    app.VideoPlayPauseButton.Text = 'Pause';
                    PlayingVideo(app);
                case 'Pause'
                    app.VideoPlayPauseButton.Text = 'Play';
            end  
        end
    end
    
    
    methods (Access = public)
        
        % Construct app
        function app = VideoPlayer
            
            % Create UIFigure and components
            createComponents(app);
            
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