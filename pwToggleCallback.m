function pwToggleCallback(hObject, eventdata)

pw_image = evalin('base','pw_image'); 
Control = evalin('base','Control');
Control(1).Command = 'set&Run';
Control(1).Parameters = {'Parameters',1,'startEvent',pw_image};
evalin('base','Resource.Parameters.startEvent = pw_image;');
assignin('base','Control',Control);

disp(['Jump to event: ' num2str(pw_image)])
disp('Switch to PW imaging mode.')
