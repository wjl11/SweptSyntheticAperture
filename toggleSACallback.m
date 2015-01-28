function toggleSACallback(hObject, eventdata)
state = get(hObject,'Value');
switch state
    case 0 
        bmode_start = evalin('base','bmode_start'); 
        Control = evalin('base','Control');
        Control(1).Command = 'set&Run';
        Control(1).Parameters = {'Parameters',1,'startEvent',bmode_start};
        evalin('base','Resource.Parameters.startEvent = bmode_start;');
        assignin('base','Control',Control);
        
        disp(['Jump to event: ' num2str(bmode_start)])
        disp('Return to b-mode imaging.')
        
        
    case 1
        bmode_end = evalin('base','bmode_end'); 
        Control = evalin('base','Control');
        Control(1).Command = 'set&Run';
        Control(1).Parameters = {'Parameters',1,'startEvent',bmode_end};
        evalin('base','Resource.Parameters.startEvent = bmode_end;');
        assignin('base','Control',Control);
        
        disp(['Jump to event: ' num2str(bmode_end)])
        disp('Save B-mode & initiate SA imaging.')
end
