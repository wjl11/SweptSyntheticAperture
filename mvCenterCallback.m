function mvCenterCallback(hObject, eventdata)
sweep_range = evalin('base','SERIAL.sweep_range');
sweep_limits = evalin('base','SERIAL.sweep_limits');
scan_velocity = evalin('base','SERIAL.scan_velocity');
norm_velocity = evalin('base','SERIAL.norm_velocity');
step = evalin('base','SERIAL.step');
acc_fnc = evalin('base','SERIAL.acc_fnc');
rs232Toggle = evalin('base','rs232Toggle');

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
mid_range = round((sweep_range(2)-sweep_range(1))/2)+sweep_range(1);
if rs232Toggle == 1
    % set positioning velocity
    fprintf(s_port,['Set Velocity ' num2str(norm_velocity)]);
    if strcmpi(fscanf(s_port),'Ok')
    else error('Failed to configure velocity. Aborting operation.')
    end

    % get current position of tdr
    fprintf(s_port, 'Get Position');
    tt_pos = str2double(fscanf(s_port));
    disp(['Current position: ' num2str(tt_pos) ' deg'])

    if tt_pos > mid_range
        tt_dir = 'CW';
        dest = mid_range;
        cmd = ['GoTo ' tt_dir ' ' num2str(dest)];
    %     disp(cmd)
    %     disp('*** Press any key to send command ***')
    %     pause()
        fprintf(s_port,cmd);
        if strcmpi(fscanf(s_port),'Ok')
        else error('Failed to send command. Aborting operation.')
        end
    elseif tt_pos < mid_range
        tt_dir = 'CCW';
        dest = mid_range;
        cmd = ['GoTo ' tt_dir ' ' num2str(dest)];
    %     disp(cmd)
    %     disp('*** Press any key to send command ***')
    %     pause()
        fprintf(s_port,cmd);
        if strcmpi(fscanf(s_port),'Ok')
        else error('Failed to send command. Aborting operation.')
        end
    else
        disp('Position in center. Not moving.')
    end

    assignin('base','tt_pos',tt_pos);
    assignin('base','tt_dir',tt_dir);
else
    warning('RS232 debug mode ON.')
end

