function init_div2SSA(hObject, eventdata)

div2SSA_acq = evalin('base','div2SSA_acq'); 
Control = evalin('base','Control');
Control(1).Command = 'set&Run';
Control(1).Parameters = {'Parameters',1,'startEvent',div2SSA_acq};
evalin('base','Resource.Parameters.startEvent = div2SSA_acq;');
assignin('base','Control',Control);

disp(['Jump to event: ' num2str(div2SSA_acq)])
disp('Initiate single element diverging SSA imaging.')

SSA_TYPE = 'div_single';
assignin('base','SSA_TYPE',SSA_TYPE);
