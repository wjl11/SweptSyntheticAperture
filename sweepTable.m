function sweepTable()

sweep_range = evalin('base','TABLE.sweep_range');
scan_velocity = evalin('base','TABLE.scan_velocity');
norm_velocity = evalin('base','TABLE.norm_velocity');
rs232Toggle = evalin('base','SETUP.rs232Toggle');

if rs232Toggle == 1
    s_port = evalin('base','s_port');
    % move transducer to home and start turn table scanning movement

    % set positioning velocity
    fprintf(s_port,['Set Velocity ' num2str(norm_velocity)]);
    if strcmpi(fscanf(s_port),'Ok')
    else error('Failed to configure velocity. Aborting operation.')
    end
    
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
        warning('Current position invalid. Moving TDR to home.');
        % move tdr to home position
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