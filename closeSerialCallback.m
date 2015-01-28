function closeSerialCallback(hObject, eventdata)

try s_port = evalin('base','s_port');
    evalin('base','fclose(s_port)');
    s_port = evalin('base','s_port');
    if strcmpi(s_port.Status,'closed')
        disp(['Serial port connection closed: ' get(s_port,'Name')])
    else
        disp('Failed to close connection.');
    end
    evalin('base','delete(s_port)');
    evalin('base','clear s_port');
catch
    disp('Performing instrument reset...')
    instrreset;
end

disp('Available serial ports: ')
disp(instrhwinfo('serial'))