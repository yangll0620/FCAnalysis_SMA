function areas_vis = testGUI()
areas = {'M1', 'lSMA', 'lVA', 'lVPLo','lCd', 'lVLo', 'rSMA', 'rMC','rVLo','rVPLo','rVA','STN','GP'};
fig = uifigure('Name','Select Visual Location','Position', [500 500 200 300]);
lbx = uilistbox(fig, ...
    'Position', [25 25 100 250], ...
    'Items',areas, ...
    'Multiselect', 'on' );

btn_OK = uibutton(fig,'Text','OK', ...
    'Position', [145, 20, 50, 25], ...
    'ButtonPushedFcn', @(btn,event) plotButtonPushed(btn,1));

btn_Cancle = uibutton(fig,'Text','Cancle', ...
    'Position', [145, 70, 50, 25]);

% Create the function for the ButtonPushedFcn callback
function plotButtonPushed(btn,a)
    areas_vis = lbx.Value;
    close(fig);
end

end