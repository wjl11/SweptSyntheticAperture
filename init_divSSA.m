function init_divSSA(hObject, eventdata)

divSSA_acq = evalin('base','divSSA_acq'); 
Control = evalin('base','Control');
Control(1).Command = 'set&Run';
Control(1).Parameters = {'Parameters',1,'startEvent',divSSA_acq};
evalin('base','Resource.Parameters.startEvent = divSSA_acq;');
assignin('base','Control',Control);

disp(['Jump to event: ' num2str(divSSA_acq)])
disp('Initiate diverging SSA imaging.')

SSA_TYPE = 'div_pw';
assignin('base','SSA_TYPE',SSA_TYPE);
