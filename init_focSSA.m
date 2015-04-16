function init_focSSA(hObject, eventdata)

focSSA_acq = evalin('base','focSSA_acq'); 
Control = evalin('base','Control');
Control(1).Command = 'set&Run';
Control(1).Parameters = {'Parameters',1,'startEvent',focSSA_acq};
evalin('base','Resource.Parameters.startEvent = focSSA_acq;');
assignin('base','Control',Control);

disp(['Jump to event: ' num2str(focSSA_acq)])
disp('Initiate focused SSA imaging.')

SSA_TYPE = 'foc';
assignin('base','SSA_TYPE',SSA_TYPE);