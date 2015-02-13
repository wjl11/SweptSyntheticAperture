function init_steerSSA(hObject, eventdata)

steerSSA_acq = evalin('base','steerSSA_acq'); 
Control = evalin('base','Control');
Control(1).Command = 'set&Run';
Control(1).Parameters = {'Parameters',1,'startEvent',steerSSA_acq};
evalin('base','Resource.Parameters.startEvent = steerSSA_acq;');
assignin('base','Control',Control);

disp(['Jump to event: ' num2str(steerSSA_acq)])
disp('Initiate steered SA imaging.')

SSA_TYPE = 'steer_pw';
assignin('base','SSA_TYPE',SSA_TYPE);
