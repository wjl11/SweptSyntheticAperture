function startAngledSACallback(hObject, eventdata)
sweep_range = evalin('base','SERIAL.sweep_range');
sweep_limits = evalin('base','SERIAL.sweep_limits');
scan_velocity = evalin('base','SERIAL.scan_velocity');
norm_velocity = evalin('base','SERIAL.norm_velocity');
step = evalin('base','SERIAL.step');
acc_fnc = evalin('base','SERIAL.acc_fnc');
rs232Toggle = evalin('base','rs232Toggle');

%******************** SERIAL SETUP START **********************************
if rs232Toggle == 1
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
end
%********************** SERIAL SETUP END **********************************

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
