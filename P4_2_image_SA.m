function P4_2_image_SA(RData)
persistent tt_pos_prev
sweep_range = evalin('base','SERIAL.sweep_range');
scan_range = evalin('base','SERIAL.scan_range');
sweep_limits = evalin('base','SERIAL.sweep_limits');
scan_velocity = evalin('base','SERIAL.scan_velocity');
norm_velocity = evalin('base','SERIAL.norm_velocity');
step = evalin('base','SERIAL.step');
acc_fnc = evalin('base','SERIAL.acc_fnc');

SAdata = evalin('base','SAdata');
SAVE_STATE = evalin('base','SAVE_STATE');

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
if isempty(tt_pos_prev)
    tt_pos_prev = sweep_range(1);
end
tic
fprintf(s_port, 'Get Position');
tt_pos = str2double(fscanf(s_port));
disp(['Current position: ' num2str(tt_pos) ' deg'])
try tt_dir = evalin('base','tt_dir');
catch; error('Direction not set.');
end
toc
% add b-mode scanning if within certain range return to b-mode event [edit]

% if TDR within scan range, save rf and position into SAdata
if tt_pos >= scan_range(1) && tt_pos <= scan_range(2)
    disp('Process SA frames')
    if isempty(SAdata.rf)
        SAdata.rf(:,:,1) = RData(:,[1:32 97:128]);
        SAdata.th = tt_pos;
        SAdata.t(:,1) = clock;
    else
        i = length(SAdata.th)+1;
        SAdata.rf(:,:,i) = RData(:,[1:32 97:128]);
        SAdata.th(i) = tt_pos;
        SAdata.t(:,i) = clock;
    end
end

% check of end of scan range has been reached
Control = evalin('base','Control');
Control(1).Command = 'set&Run';
if tt_pos == sweep_range(1) || tt_pos == sweep_range(2) 
    
    if (SAVE_STATE(1) == 0 || SAVE_STATE(2) == 0) && tt_pos ~= tt_pos_prev && ~isnan(tt_pos_prev) 
        % save b-mode and SA data
        save_start = evalin('base','save_start'); 
        Control(1).Parameters = {'Parameters',1,'startEvent',save_start};
        evalin('base','Resource.Parameters.startEvent = save_start;');
        assignin('base','Control',Control);
        
    else
%         if SAVE_STATE(1) == 1 && SAVE_STATE(2) == 1
        
        if evalin('base','get(UI(1).handle,''Value'')') == 0  
            % return to guidance b-mode ultrasound if toggle turned off
            bmode_start = evalin('base','bmode_start'); 
            Control(1).Parameters = {'Parameters',1,'startEvent',bmode_start};
            evalin('base','Resource.Parameters.startEvent = bmode_start;');
            assignin('base','Control',Control);

            disp(['Jump to event: ' num2str(bmode_start)])
            disp('Return to b-mode imaging.')
        elseif evalin('base','get(UI(1).handle,''Value'')') == 1
            % continue SA sequence if toggle remains on
            bmode_end = evalin('base','bmode_end'); 
            Control(1).Parameters = {'Parameters',1,'startEvent',bmode_end};
            evalin('base','Resource.Parameters.startEvent = bmode_end;');
            assignin('base','Control',Control);

            disp(['Jump to event: ' num2str(bmode_end)])
            disp('Save B-mode & initiate SA imaging.')

            % switch direction based on current position
            if strcmpi(tt_dir,'CW') && tt_pos == sweep_range(1)
                tt_dir = 'CCW';
                cmd = ['GoTo ' tt_dir ' ' num2str(sweep_range(2))];
                % move to next position
                disp(cmd)
                fprintf(s_port,cmd);
                if strcmpi(fscanf(s_port),'Ok')
                else error('Failed to return home. Aborting operation.')
                end
            elseif strcmpi(tt_dir,'CCW') && tt_pos == sweep_range(2)
                tt_dir = 'CW';
                cmd = ['GoTo ' tt_dir ' ' num2str(sweep_range(1))];
                % move to next position
                disp(cmd)
                fprintf(s_port,cmd);
                if strcmpi(fscanf(s_port),'Ok')
                else error('Failed to return home. Aborting operation.')
                end
            end
            
%             disp('*** Press any key to send command ***')
%             pause()
        end
        SAVE_STATE = [0,0];
    end
else % continue SA sequence if scan range not reached
    SA_start = evalin('base','SA_start');
    Control(1).Parameters = {'Parameters',1,'startEvent',SA_start};
    evalin('base','Resource.Parameters.startEvent = SA_start;');
    assignin('base','Control',Control);
%     disp('Range not reached. Continue SA imaging.')
end
tt_pos_prev = tt_pos;
assignin('base','tt_pos',tt_pos);
assignin('base','tt_dir',tt_dir);
assignin('base','SAdata',SAdata);
assignin('base','SAVE_STATE',SAVE_STATE);

% try count = evalin('base','count');
% catch
%     count = 1;
% end
% 
% count = count+1;
% if count < 10
%     disp(count)
%     SA_start = evalin('base','SA_start');
%     Control = evalin('base','Control');
%     Control(1).Command = 'set&Run';
%     Control(1).Parameters = {'Parameters',1,'startEvent',SA_start};
%     evalin('base','Resource.Parameters.startEvent = SA_start;');
%     assignin('base','Control',Control);
%     imagesc(RData)
%     
% else
%     state = evalin('base','get(UI(1).handle,''Value'')')
%     if state == 0
%         bmode_start = evalin('base','bmode_start'); 
%         Control = evalin('base','Control');
%         Control(1).Command = 'set&Run';
%         Control(1).Parameters = {'Parameters',1,'startEvent',bmode_start};
%         evalin('base','Resource.Parameters.startEvent = bmode_start;');
%         assignin('base','Control',Control);
% 
%         disp(['Jump to event: ' num2str(bmode_start)])
%         disp('Return to b-mode imaging.')
%     end
%     count = 1;
%     disp('New sweep');
% %     toggleSACallback(evalin('base','hObject'),evalin('base','eventdata'))
% end
% % 
% assignin('base','count',count)