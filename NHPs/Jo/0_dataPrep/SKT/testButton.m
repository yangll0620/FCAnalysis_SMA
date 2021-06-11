function testButton()

% auto save without click
reply = input('Auto save mode [y] or no [n] ', 's');
if isempty(reply) || ~strcmpi(reply, 'y')
    autoSaveMode =  'n';
else
    autoSaveMode = 'y';
end
clear reply


% subplot/Figure parameters
fig_left = 1;
fig_bottom = 41;
fig_width = 1920;
fig_height = 960;


subp_endLeft = 0.97;


% position of Next and previous button
btn_width = 50;
btn_height = 30;
pos_btn_next_left = (subp_endLeft) * fig_width;
pos_btn_next_bottom = (0.05) * fig_height;
pos_btn_prev_left = pos_btn_next_left;
pos_btn_prev_bottom = pos_btn_next_bottom + btn_height + 10;

% Finish Button
pos_btn_finish_left = (subp_endLeft) * fig_width;
pos_btn_finish_bottom = 0.5 * fig_height;


% create a new figure
tagnum = 1;
tag = ['fig' num2str(tagnum)];
fig = figure('Tag', tag, 'Name', tag, 'NumberTitle','off');
annotation(fig,'textbox',...
    [0.35 0.45 0.27 0.2],...
    'String',{tag});


% add next button
c_next = uicontrol(fig, 'Style','pushbutton', 'String', 'Next', 'Position', [450 20 60 30]);
c_next.Callback = @btn_next;

max_figN = 5;
if strcmpi(autoSaveMode, 'y')
    
    % extract current fig tag and save it
    fig_current = gcf;
    tag = fig_current.Tag;
    saveas(gcf, tag, 'tif');
    
    nextFig_isLast = func_next();
    while(~nextFig_isLast)
        
        % extract current fig tag and save it
        fig_current = gcf;
        tag = fig_current.Tag;
        saveas(gcf, tag, 'tif');
        
        % show next fig
        nextFig_isLast = func_next();        
    end
    
    % extract last fig tag and save it
    fig_current = gcf;
    tag = fig_current.Tag;
    saveas(gcf, tag, 'tif');
    func_finish();
end


end


function nextFig_isLast = func_next()
% return 0 if the last figure, else 1

fig_current = gcf;
tag = fig_current.Tag;
tagnum = str2num(tag(length('fig')+1 : end));
close(fig_current);

%%% show next figure
tagnum = tagnum + 1;
tag = ['fig' num2str(tagnum)];
fig_next = figure('Tag', tag, 'Name', tag, 'NumberTitle','off');
annotation(fig_next,'textbox',...
    [0.35 0.45 0.27 0.2],...
    'String',{tag},...
    'FitBoxToText','off');

nextFig_isLast = false;
max_figN = 5;
if tagnum < max_figN
    % add next button
    c_next = uicontrol(fig_next, 'Style','pushbutton', 'String', 'Next', 'Position', [450 20 60 30]);
    c_next.Callback = @btn_next;
    
    % add previous button
    c_prev = uicontrol(fig_next, 'Style','pushbutton', 'String', 'Previous', 'Position', [50 20 60 30]);
    c_prev.Callback = @btn_prev;
    
    nextFig_isLast = false;
end

if tagnum == max_figN
    
    % add previous button
    c_prev = uicontrol(fig_next, 'Style','pushbutton', 'String', 'Previous', 'Position', [50 20 60 30]);
    c_prev.Callback = @btn_prev;
    
    % add Finish button
    c_Finish = uicontrol(fig_next, 'Style','pushbutton', 'String', 'Finish', 'Position', [250 120 60 30]);
    c_Finish.Callback = @btn_finish;
    
    nextFig_isLast = true;
end
end

function btn_next(src,event)
func_next();
end %


function func_prev()
fig_current = gcf;
tag = fig_current.Tag;
tagnum = str2num(tag(length('fig')+1 : end));
close(fig_current);

%%% show previous figure
tagnum = tagnum - 1;
tag = ['fig' num2str(tagnum)];
fig_prev = figure('Tag', tag, 'Name', tag, 'NumberTitle','off');
annotation(fig_prev,'textbox',...
    [0.35 0.45 0.27 0.2],...
    'String',{tag},...
    'FitBoxToText','off');

if tagnum >= 2
    % add next button
    c_next = uicontrol(fig_prev, 'Style','pushbutton', 'String', 'Next', 'Position', [450 20 60 30]);
    c_next.Callback = @btn_next;
    
    % add previous button
    c_prev = uicontrol(fig_prev, 'Style','pushbutton', 'String', 'Previous', 'Position', [50 20 60 30]);
    c_prev.Callback = @btn_prev;
end

if tagnum == 1
    % add next button
    c_next = uicontrol(fig_prev, 'Style','pushbutton', 'String', 'Next', 'Position', [450 20 60 30]);
    c_next.Callback = @btn_next;
end

end

function btn_prev(src,event)
func_prev();
end

function func_finish()
fig_current = gcf;
tag = fig_current.Tag;
tagnum = str2num(tag(length('fig')+1 : end));
close(fig_current);
end
function btn_finish(src,event)
func_finish();
end


