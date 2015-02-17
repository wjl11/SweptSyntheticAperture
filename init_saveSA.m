function init_saveSA(hObject,eventdata)

SA_acq = evalin('base','SA_acq'); 
Control = evalin('base','Control');

Control(1).Command = 'set&Run';
Control(1).Parameters = {'Parameters',1,'startEvent',SA_acq};
evalin('base','Resource.Parameters.startEvent = SA_acq;');
assignin('base','Control',Control);

disp(['Jump to event: ' num2str(SA_acq)])