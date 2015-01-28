function toggleSACallback(hObject, eventdata)
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

state = get(hObject,'Value');
switch state
    case 0 % return to b-mode
        
        % if not at home position skip section below and continue scanning
        
        % get current position of tdr
        fprintf(s_port, 'Get Position');
        tt_pos = str2double(fscanf(s_port));
        disp(['Current position: ' num2str(tt_pos) ' deg'])
        try tt_dir = evalin('base','tt_dir');
        catch tt_dir = '';
        end
        if tt_pos == sweep_range(1) || tt_pos == sweep_range(2)
            bmode_start = evalin('base','bmode_start'); 
            Control = evalin('base','Control');
            Control(1).Command = 'set&Run';
            Control(1).Parameters = {'Parameters',1,'startEvent',bmode_start};
            evalin('base','Resource.Parameters.startEvent = bmode_start;');
            assignin('base','Control',Control);

            disp(['Jump to event: ' num2str(bmode_start)])
            disp('Return to b-mode imaging.')
        else
            disp('Finishing SA scan cycle...')
        end
%        count = evalin('base','count');
%        if count < 10
%             bmode_start = evalin('base','bmode_start'); 
%             Control = evalin('base','Control');
%             Control(1).Command = 'set&Run';
%             Control(1).Parameters = {'Parameters',1,'startEvent',bmode_start};
%             evalin('base','Resource.Parameters.startEvent = bmode_start;');
%             assignin('base','Control',Control);
% 
%             disp(['Jump to event: ' num2str(bmode_start)])
%             disp('Return to b-mode imaging.')
%        end
       
    case 1 % jump to SA imaging
        bmode_end = evalin('base','bmode_end'); 
        Control = evalin('base','Control');
        Control(1).Command = 'set&Run';
        Control(1).Parameters = {'Parameters',1,'startEvent',bmode_end};
        evalin('base','Resource.Parameters.startEvent = bmode_end;');
        assignin('base','Control',Control);
        
        disp(['Jump to event: ' num2str(bmode_end)])
        disp('Save B-mode & initiate SA imaging.')
        
        % move transducer to home and start turn table scanning movement
        
        % set move to home velocity
        fprintf(s_port,['Set Velocity ' num2str(norm_velocity)]);
        if strcmpi(fscanf(s_port),'Ok')
        else error('Failed to configure velocity. Aborting operation.')
        end
        
        % get current position of tdr
        fprintf(s_port, 'Get Position');
        tt_pos = str2double(fscanf(s_port));
        disp(['Current position: ' num2str(tt_pos) ' deg'])
        if tt_pos > sweep_range(1)
            tt_dir = 'CW';
        elseif tt_pos < sweep_range(1)
            tt_dir = 'CCW';
        else
            disp('TDR position at home.')
            tt_dir = 'CCW';
        end
        
        % move tdr to home position
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
        
        % begin sweep for scan
        if tt_pos ~= sweep_range(1)
            error('Position mismatch.')
        end
        
         % set scan velocity
        fprintf(s_port,['Set Velocity ' num2str(scan_velocity)]);
        if strcmpi(fscanf(s_port),'Ok')
        else error('Failed to configure velocity. Aborting operation.')
        end
        
        tt_dir = 'CCW';
        cmd = ['GoTo ' tt_dir ' ' num2str(sweep_range(2))];
        disp(cmd)
        disp('*** Press any key to send command ***')
        pause()
        fprintf(s_port,cmd);    
end

assignin('base','tt_pos',tt_pos);
assignin('base','tt_dir',tt_dir);
