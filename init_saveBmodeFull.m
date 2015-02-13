function init_saveBmodeFull(hObject, eventdata)

bmode_acq = evalin('base','bmode_acq'); 
Control = evalin('base','Control');

Control(1).Command = 'set&Run';
Control(1).Parameters = {'Parameters',1,'startEvent',bmode_acq};
evalin('base','Resource.Parameters.startEvent = bmode_acq;');
assignin('base','Control',Control);

disp(['Jump to event: ' num2str(bmode_acq)])
