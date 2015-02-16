function table_mvCenter(hObject, eventdata)

sweep_range = evalin('base','TABLE.sweep_range');
norm_velocity = evalin('base','TABLE.norm_velocity');
rs232Toggle = evalin('base','SETUP.rs232Toggle');

mid_range = round((sweep_range(2)-sweep_range(1))/2)+sweep_range(1);
if rs232Toggle == 1
    s_port = evalin('base','s_port');
        
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
    elseif tt_pos == mid_range
        disp('Position in center. Move to sweep position.')
        tt_dir = 'CW';
        dest = sweep_range(1);
        cmd = ['GoTo ' tt_dir ' ' num2str(dest)];
        
        disp(cmd)
        disp('*** Press any key to send command ***')
        pause()
        fprintf(s_port,cmd);
        if strcmpi(fscanf(s_port),'Ok')
        else error('Failed to send command. Aborting operation.')
        end
    end

    assignin('base','tt_pos',tt_pos);
    assignin('base','tt_dir',tt_dir);
else
    warning('RS232 debug mode ON.')
end

