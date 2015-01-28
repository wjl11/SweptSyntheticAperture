function abortSweepCallback(hObject, eventdata)

%******************** SERIAL SETUP START **********************************
try s_port = evalin('base','s_port');
catch
    devices = instrhwinfo('serial');
    while isempty(devices.SerialPorts)
        disp('Searching for serial port...')
        devices = instrhwinfo('serial');
        pause(0.5)
    end
    s_port = serial(devices.SerialPorts,...
     'BaudRate',9600,...
     'DataBits',8,...
     'Parity','none',...
     'Terminator','NUL',...
     'StopBits',2,...
     'Timeout',0.5,...
     'InputBufferSize', 4096);
    
    fopen(s_port);
    if strcmpi(s_port.Status,'open')
        disp(['Serial port connected: ' get(s_port,'Name')])
        s_port
    else
        error('Serial connection failed. Aborting operation.');
    end
    s_port.ReadAsyncMode = 'continuous';
    assignin('base','s_port',s_port);
    
    % set to bipolar mode
    fprintf(s_port,'Set DisplayPolarity BIPOLAR');
    if strcmpi(fscanf(s_port),'Ok')
    else error('Failed to configure polarity. Aborting operation.')
    end
    
    % set acc function
    fprintf(s_port,['Set AccelFunc ' num2str(acc_fnc)]);
    if strcmpi(fscanf(s_port),'Ok')
    else error('Failed to configure acceleration function. Aborting operation.')
    end
end

if strcmpi(s_port.Status,'open')
else error('Serial connection failed. Aborting operation.');
end
%********************** SERIAL SETUP END **********************************

fprintf(s_port,'Set MoveAbort');
if strcmpi(fscanf(s_port),'Ok')
    disp('Action aborted.')
else
    error('Command not sent.')
end