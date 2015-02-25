function init_pwSSA(hObject, eventdata)

SSA_acq = evalin('base','SSA_acq'); 
Control = evalin('base','Control');
Control(1).Command = 'set&Run';
Control(1).Parameters = {'Parameters',1,'startEvent',SSA_acq};
evalin('base','Resource.Parameters.startEvent = SSA_acq;');
assignin('base','Control',Control);

disp(['Jump to event: ' num2str(SSA_acq)])
disp('Initiate planewave SSA imaging.')

SSA_TYPE = 'pw';
assignin('base','SSA_TYPE',SSA_TYPE);
