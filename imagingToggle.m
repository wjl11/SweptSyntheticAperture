function imagingToggle(hObject, eventdata)

IM_STATE = evalin('base','IM_STATE');
pw_image = evalin('base','pw_image');
bmode_image = evalin('base','bmode_image');
Control = evalin('base','Control');

if strcmpi(IM_STATE,'pa')
    Control(1).Command = 'set&Run';
    Control(1).Parameters = {'Parameters',1,'startEvent',pw_image};
    evalin('base','Resource.Parameters.startEvent = pw_image;');
    
    Control(2).Command = 'set&Run';
    Control(2).Parameters = {'DisplayWindow',1,'clrWindow',1};
    evalin('base','Resource.DisplayWindow(1).clrWindow = 1;');
    
    assignin('base','Control',Control);

    disp(['Jump to event: ' num2str(pw_image)])
    disp('Switch to PW imaging mode.')
    IM_STATE = 'pw';
elseif strcmpi(IM_STATE,'pw')
    Control(1).Command = 'set&Run';
    Control(1).Parameters = {'Parameters',1,'startEvent',bmode_image};
    evalin('base','Resource.Parameters.startEvent = bmode_image;');
    
    Control(2).Command = 'set&Run';
    Control(2).Parameters = {'DisplayWindow',1,'clrWindow',1};
    evalin('base','Resource.DisplayWindow(1).clrWindow = 1;');
    
    assignin('base','Control',Control);

    disp(['Jump to event: ' num2str(bmode_image)])
    disp('Switch to PA B-mode imaging mode.')
    IM_STATE = 'pa';
else
    error('Imaging state invalid.')
end
assignin('base','IM_STATE',IM_STATE);
