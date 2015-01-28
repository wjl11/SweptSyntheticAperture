function sweep_probe()

sweep_range = evalin('base','SERIAL.sweep_range');
sweep_limits = evalin('base','SERIAL.sweep_limits');
scan_velocity = evalin('base','SERIAL.scan_velocity');
norm_velocity = evalin('base','SERIAL.norm_velocity');
step = evalin('base','SERIAL.step');
acc_fnc = evalin('base','SERIAL.acc_fnc');
rs232Toggle = evalin('base','rs232Toggle');

if rs232Toggle == 1
    s_port = evalin('base','s_port');
    % move transducer to home and start turn table scanning movement

    % set positioning velocity
    fprintf(s_port,['Set Velocity ' num2str(norm_velocity)]);
    if strcmpi(fscanf(s_port),'Ok')
    else error('Failed to configure velocity. Aborting operation.')
    end
    pause(0.5)
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
        disp('Current position invalid. Moving TDR to home.');
        if tt_pos > sweep_range(1)
            tt_dir = 'CW';
        elseif tt_pos < sweep_range(1)
            tt_dir = 'CCW';
        else
            tt_dir = 'CCW';
        end

        fprintf(s_port,['Set Velocity ' num2str(norm_velocity)]);
        if strcmpi(fscanf(s_port),'Ok')
        else error('Failed to configure velocity. Aborting operation.')
        end

        cmd = ['GoTo ' tt_dir ' ' num2str(sweep_range(1))];
    %     disp(cmd)
    %     disp('*** Press any key to send command ***')
    %     pause()
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

    cmd = ['GoTo ' tt_dir ' ' num2str(dest)];
    % disp(cmd)
    % disp('*** Press any key to send command ***')
    % pause()
    fprintf(s_port,cmd);
    if strcmpi(fscanf(s_port),'Ok')
    else error('Failed to send command. Aborting operation.')
    end

    assignin('base','tt_pos',tt_pos);
    assignin('base','tt_dir',tt_dir);
else 
    warning('RS232 debug mode ON.')
end