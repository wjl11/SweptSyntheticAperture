function testSweepCallback(hObject, eventdata)

sweep_range = evalin('base','SERIAL.sweep_range');
scan_range = evalin('base','SERIAL.scan_range');
sweep_limits = evalin('base','SERIAL.sweep_limits');
scan_velocity = evalin('base','SERIAL.scan_velocity');
norm_velocity = evalin('base','SERIAL.norm_velocity');
step = evalin('base','SERIAL.step');
acc_fnc = evalin('base','SERIAL.acc_fnc');

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

% get current position of tdr
fprintf(s_port, 'Get Position');
tt_pos = str2double(fscanf(s_port));
disp(['Current position: ' num2str(tt_pos) ' deg'])
if tt_pos == sweep_range(1)
    tt_dir = 'CCW';
    dest = sweep_range(2);
elseif tt_pos == sweep_range(2)
    tt_dir = 'CW';
    dest = sweep_range(1);
else
    warning('Current position invalid. Move TDR to home.');
    % move tdr to home position
    if tt_pos > sweep_range(1)
        tt_dir = 'CW';
    elseif tt_pos < sweep_range(1)
        tt_dir = 'CCW';
    else
        disp('TDR position at home.')
        tt_dir = 'CCW';
    end
    fprintf(s_port,['Set Velocity ' num2str(norm_velocity)]);
    if strcmpi(fscanf(s_port),'Ok')
    else error('Failed to configure velocity. Aborting operation.')
    end

    cmd = ['GoTo ' tt_dir ' ' num2str(sweep_range(1))];
    disp(cmd)
    disp('*** Press any key to send command ***')
    pause()
    fprintf(s_port,cmd);
    if strcmpi(fscanf(s_port),'Ok')
    else error('Failed to return home. Aborting operation.')
    end

    while tt_pos ~= sweep_range(1)
        fprintf(s_port,'Get Position');
        tt_pos = str2double(fscanf(s_port));
        disp(['TDR position: ' num2str(tt_pos) ' deg'])
    end
    
    tt_dir = 'CCW';
    dest = sweep_range(2);
end

% set sweep velocity
fprintf(s_port,['Set Velocity ' num2str(scan_velocity)]);
if strcmpi(fscanf(s_port),'Ok')
else error('Failed to configure velocity. Aborting operation.')
end

% move tdr to opposite side
cmd = ['GoTo ' tt_dir ' ' num2str(dest)];
disp(cmd)
disp('*** Press any key to send command ***')
pause()
fprintf(s_port,cmd);
if strcmpi(fscanf(s_port),'Ok')
else error('Failed to return home. Aborting operation.')
end