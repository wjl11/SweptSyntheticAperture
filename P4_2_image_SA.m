function P4_2_image_SA(RData)
persistent frameNum
sweep_range = evalin('base','SERIAL.sweep_range');
scan_range = evalin('base','SERIAL.scan_range');
sweep_limits = evalin('base','SERIAL.sweep_limits');
scan_velocity = evalin('base','SERIAL.scan_velocity');
norm_velocity = evalin('base','SERIAL.norm_velocity');
step = evalin('base','SERIAL.step');
acc_fnc = evalin('base','SERIAL.acc_fnc');

SAframes = evalin('base','SAframes');
SArows = evalin('base','rowsPerFrameSA');
SAcols = 64;
SAdata = evalin('base','SAdata');
tic
if isempty(frameNum)
    frameNum = 1;
end

if isempty(SAdata.rf)
    SAdata.rf = zeros(SArows,SAcols,SAframes);
    SAdata.t = zeros(length(clock),SAframes);
    SAdata.rf(:,:,frameNum) = RData(:,[1:32 97:128]);
    SAdata.t(:,frameNum) = clock;
else
    SAdata.rf(:,:,frameNum) = RData(:,[1:32 97:128]);
    SAdata.t(:,frameNum) = clock;
end

% disp(['SA frame ' num2str(frameNum) '/' num2str(SAframes)])
if frameNum == SAframes
    % save all frames
%     disp(['Saving ' num2str(SAframes) ' SA frames.'])
%     save_start = evalin('base','save_start'); 
%     Control(1).Parameters = {'Parameters',1,'startEvent',save_start};
%     evalin('base','Resource.Parameters.startEvent = save_start;');
%     assignin('base','Control',Control);
    frameNum = 1;
else 
    % continue collecting frames
    SA_start = evalin('base','SA_start');
    Control = evalin('base','Control');
    Control(1).Command = 'set&Run';
    Control(1).Parameters = {'Parameters',1,'startEvent',SA_start};
    evalin('base','Resource.Parameters.startEvent = SA_start;');
    assignin('base','Control',Control);
    frameNum = frameNum+1;
end
assignin('base','SAdata',SAdata);
toc