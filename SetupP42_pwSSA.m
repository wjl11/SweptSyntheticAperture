clear all; close all;
addpath(genpath('/home/wjl11/Matlab Simulator/'));
% simulation mode y/n
simToggle = 0; % 1 [sim], 0 [probe]
rs232Toggle = 1; % 1 [table com on], 0 [table com off]

% general imaging parameters
maxVoltage = 50;
numEl = 64; % number of physical elements (64 elements in the P4-2, only half channels of vsx used)

% general SSA parameters
ssaPRT = 1e3; % time between SSA planewave acquisitions
disp(['Plane wave SA frame rate: ' num2str(1/(ssaPRT*1e-3)) ' kHz']);
rowsPerFrameSA = 4096; 
rs232wait = 0; % wait time for rs232 sweep command [ms] (max at 104.8 ms)
ssaEndDepthMM = 150;
ssaStartDepthMM =[];

% planewave SSA parameters

% angled planewave SSA paramters
pwAngles = [-20:20:20];
na = length(pwAngles);

% bmode parameters
nr = 128; % no of ray lines in b-mode scan
rowsPerFrameBmode = nr*4096; % rows of RF data stored
PRF = 23; % determines time between b-mode frames
bmodeFocusMM = 100;
bmodeEndDepthMM = 150;
txfnum = 2;
tPulse = 250; % time between phased array acqs [us]
c = 1540;

% pw imaging paramters
pwPRT = 1e4; % time between planewave imaging frames

disp(['Fundamental frame rate: ' num2str(1/(ssaEndDepthMM*2/c)) ' kHz']);

% Frame buffers
imgBuffer = 10;
bmodeSaveFrames = 2;

% Turn table and serial parameters
SERIAL.sweep_range = [0.0 45.0];
SERIAL.sweep_limits = [0.0 45.0]; % DO NOT ALTER
SERIAL.scan_velocity = 5.0;
SERIAL.norm_velocity = 3.0;
SERIAL.step = 1.0;
SERIAL.acc_fnc = 3;
% 0 - impulse
% 1 - flat
% 2 - ramp
% 3 - sin 
if SERIAL.sweep_range(1) < SERIAL.sweep_limits(1) || ...
        SERIAL.sweep_range(2) > SERIAL.sweep_limits(2) || ...
        SERIAL.sweep_range(1) > SERIAL.sweep_range(2) ||...
        SERIAL.sweep_limits(1) > SERIAL.sweep_limits(2) 
    error('Error in TDR sweep bounds.')
end

%calculation of frames per sweep
ssaFrames = round(((SERIAL.sweep_range(2)-SERIAL.sweep_range(1))/360.0)/(SERIAL.scan_velocity/60)/(ssaPRT*1e-6))+500;
disp(['SA frames per acquisition: ' num2str(ssaFrames)])

% Save properties
saveLabel = input('Save file label: ','s');

IM_STATE = 'pa';
SSA_TYPE = [];
% pa = phased array
% pw = plane wave

% specify media points
% Media.MP(1,:) = [0,0,100,1.0]; % specify point in media [x,y,z,reflectivity]
% Media.MP(2,:) = [0,-5,100,1.0];
% Media.MP(3,:) = [0,5,100,1.0];
% Media.numPoints = 3;

% Specify Moving Media.  Use point targets in middle of PData.
% Set up Media points
Media.MP(1,:) = [-45,0,30,1.0];
Media.MP(2,:) = [-15,0,30,1.0];
Media.MP(3,:) = [15,0,30,1.0];
Media.MP(4,:) = [45,0,30,1.0];
Media.MP(5,:) = [-15,0,60,1.0];
Media.MP(6,:) = [-15,0,90,1.0];
Media.MP(7,:) = [-15,0,120,1.0];
Media.MP(8,:) = [-15,0,150,1.0];
Media.MP(9,:) = [-45,0,120,1.0];
Media.MP(10,:) = [15,0,120,1.0];
Media.MP(11,:) = [45,0,120,1.0];
Media.MP(12,:) = [-10,0,69,1.0];
Media.MP(13,:) = [-5,0,75,1.0];
Media.MP(14,:) = [0,0,78,1.0];
Media.MP(15,:) = [5,0,80,1.0];
Media.MP(16,:) = [10,0,81,1.0];
Media.MP(17,:) = [-75,0,120,1.0];
Media.MP(18,:) = [75,0,120,1.0];
Media.MP(19,:) = [-15,0,180,1.0];
Media.numPoints = 19;
Media.function = 'movePoints';

% specify system parameters
Resource.Parameters.numTransmit = 128;
Resource.Parameters.numRcvChannels = 256; % no of receive channels for 4 board VDAS unit (changed from 128 to 256 for vantage system)
% Resource.Parameters.connector = 1; % tdr connector to use for 4 board VDAS unit
Resource.Parameters.speedOfSound = c;
Resource.Parameters.simulateMode = simToggle;
if Resource.Parameters.simulateMode ~= 0
    disp('Simulation mode: ON')
else
    warning('Simulation mode: OFF')
end
if rs232Toggle ~= 1
    disp('Serial communication: OFF')
else
    warning('Serial communication: ON')
end

% computeTrans function used to define the transducer struct
freqMHz = 2.50;
Trans.name = 'P4-2';
Trans.frequency = freqMHz; % not needed if using default center freq -
% assume center frequency
Trans.units = 'mm';
Trans = computeTrans(Trans);
Trans.maxHighVoltage = maxVoltage;

% Set up SFormat structure array for bmode
aperture = numEl*Trans.spacing;
SFormat(1).transducer = 'P4-2';
SFormat(1).scanFormat = 'VAPX'; 
SFormat(1).theta = -pi/4;
SFormat(1).radius = (aperture/2)/tan(-SFormat(1).theta); % dist. to virt. apex
SFormat(1).numRays = nr;      % no. of Rays
SFormat(1).FirstRayLoc = Trans.ElementPos(1,1:3);   % x,y,z
SFormat(1).rayDelta = 2*(-SFormat(1).theta)/(nr-1);  % spacing in radians(sector) or dist. between rays
SFormat(1).startDepth = 0;
SFormat(1).endDepth = round(bmodeEndDepthMM/1000/(c/(Trans.frequency*1e6))); % Acquisition depth in wavelengths

% Set up SFormat structure array for pw imaging
SFormat(2).transducer = 'P4-2';
SFormat(2).scanFormat = 'RLIN'; 
SFormat(2).theta = 0;
SFormat(2).radius = 0; 
SFormat(2).numRays = 1;
SFormat(2).FirstRayLoc = [0,0,0];   % x,y,z
SFormat(2).rayDelta = aperture;  % spacing in radians(sector) or dist. between rays
SFormat(2).startDepth = 0; %[edit] revise start and end depth to reduce stored data for ssa
SFormat(2).endDepth = round(bmodeEndDepthMM/1000/(c/(Trans.frequency*1e6))); % [edit]

% Set up PData structure for bmode. Pixels representing phased array image
PData(1).sFormat = 1;
PData(1).pdelta = 0.875;
PData(1).Size(1) = 10 + ceil((SFormat(1).endDepth-SFormat(1).startDepth)/PData(1).pdelta);
PData(1).Size(2) = 10 + ceil(2*(SFormat(1).endDepth + SFormat(1).radius)*sin(-SFormat(1).theta)/PData(1).pdelta);
PData(1).Size(3) = 1;
PData(1).Origin = [-(PData(1).Size(2)/2)*PData(1).pdelta,0,0];

% Set up PData structure for pw imaging. Pixels representing flash image
PData(2).sFormat = 2;
PData(2).pdeltaX = Trans.spacing;
PData(2).pdeltaZ = 0.5; % smaller values = higher resolution (longer reconstruction time)
PData(2).Size(1) = ceil((SFormat(2).endDepth-SFormat(2).startDepth)/PData(2).pdeltaZ);
PData(2).Size(2) = ceil(aperture/PData(2).pdeltaX);
PData(2).Size(3) = 1;
PData(2).Origin = [-Trans.spacing*(numEl-1)/2,0,SFormat(2).startDepth];

% resource buffer for image display (bmode and PW)
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = rowsPerFrameBmode; 
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = imgBuffer;

% b-mode image buffer and image display window
Resource.ImageBuffer(1).datatype = 'double';
Resource.ImageBuffer(1).rowsPerFrame = PData(1).Size(1); % display all pixels
Resource.ImageBuffer(1).colsPerFrame = PData(1).Size(2);
Resource.ImageBuffer(1).numFrames = imgBuffer;

Resource.DisplayWindow(1).Title = 'Bmode Display';
Resource.DisplayWindow(1).pdelta = 0.3;
Resource.DisplayWindow(1).Position = [200,200, ...
    ceil(PData(1).Size(2)*(PData(1).pdelta/Resource.DisplayWindow(1).pdelta)), ...
    ceil(PData(1).Size(1)*PData(1).pdelta/Resource.DisplayWindow(1).pdelta)];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),PData(1).Origin(3)];
Resource.DisplayWindow(1).Colormap = gray(256);

% pw image buffer and image display window
Resource.ImageBuffer(2).datatype = 'double';
Resource.ImageBuffer(2).rowsPerFrame = ceil((rowsPerFrameSA/8)/PData(2).pdeltaZ); % display all pixels
Resource.ImageBuffer(2).colsPerFrame = PData(2).Size(2);
Resource.ImageBuffer(2).numFrames = imgBuffer;

% Resource.DisplayWindow(2).Title = 'PW Display';
% Resource.DisplayWindow(2).pdelta = 0.3;
%     ScrnSize = get(0,'ScreenSize');
%     DwWidth = ceil(PData(2).Size(2)*PData(2).pdeltaX/Resource.DisplayWindow(2).pdelta);
%     DwHeight = ceil(PData(2).Size(1)*PData(2).pdeltaZ/Resource.DisplayWindow(2).pdelta);
% Resource.DisplayWindow(2).Position = [250,(ScrnSize(4)-(DwHeight+150))/2,DwWidth,DwHeight];
% Resource.DisplayWindow(2).ReferencePt = [PData(2).Origin(1),PData(2).Origin(3)];
% Resource.DisplayWindow(2).Colormap = gray(256);

% planewave SA resource buffer
Resource.RcvBuffer(2).datatype = 'int16';
Resource.RcvBuffer(2).rowsPerFrame = rowsPerFrameSA; % number of samples of RF data per channel, default sample rate for A/D 
% is 4x the center frequency of the probe
Resource.RcvBuffer(2).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(2).numFrames = ssaFrames;

% buffer for saving phased array bmode data
Resource.RcvBuffer(3).datatype = 'int16';
Resource.RcvBuffer(3).rowsPerFrame = rowsPerFrameBmode; 
Resource.RcvBuffer(3).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(3).numFrames = bmodeSaveFrames;

% transmit waveform for both SA and bmode
TW.type = 'parametric';
TW.Parameters = [Trans.frequency,.67,2,1];

% TW(1).Parameters = [36,17,2,1];   % A, B, C, D [edit]
% A,B,C,D in terms of master clock (180 MHz T = 5.556 ns)
% A = clocks in half cycle based on f0
% B = clocks during which transmit drivers are active
% C = number of half cycles in transmit waveform (2 = single cycle)
% D = polarity of first half cycle (1 = positive)

% tx definitions for b-mode imaging
tx_foc = round(bmodeFocusMM/1000/(c/(Trans.frequency*1e6))); %focus in wl
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0,0,0], ...
                   'focus', 0, ... % in wavelengths; 
                   'Steer', [0,0], ...
                   'Apod', ones(1,numEl), ...  
                   'Delay', zeros(1,numEl)), 1, nr+1+na); 

% specify phased array bmode tx attributes 
bmodeAngles = SFormat(1).theta:SFormat(1).rayDelta:(SFormat(1).theta + (nr-1)*SFormat(1).rayDelta);
origin = SFormat(1).radius*tan(bmodeAngles);
for n = 1:nr
    TX(n).Origin = [origin(n),0.0,0.0];
    TX(n).focus = tx_foc;
    TX(n).Steer = [bmodeAngles(n),0.0];
    TX(n).Delay = computeTXDelays(TX(n));
end

% specify planewave tx attributes
TX(nr+1).Delay = computeTXDelays(TX(nr+1));

% specify angled planewave tx attributes
for n = 1:na
    TX(n+nr+1).Steer = [degtorad(pwAngles(n)),0.0]; %[edit] currently in radians
    TX(n+nr+1).Delay = computeTXDelays(TX(n+nr+1));
end

% specify transmit power controllers (optional)
TPC(1).name = 'B-mode';
TPC(1).maxHighVoltage = maxVoltage;
TPC(2).name = 'Swept Synthetic Aperture';
TPC(2).maxHighVoltage = maxVoltage;

% specify TGC
TGC.CntrlPts = [450,550,650,710,770,830,890,950];
TGC.rangeMax = SFormat(1).endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% rcv defintions for b-mode 
maxAcqLengthBmode = sqrt(aperture^2 + SFormat(1).endDepth^2 - 2*aperture*SFormat(1).endDepth*cos(SFormat(1).theta-pi/2)) - SFormat(1).startDepth;
maxAcqLengthPW = sqrt(SFormat(2).endDepth^2+aperture^2)-SFormat(2).startDepth;
wlsPer128 = 128/(2*4); % wavelengths in 128 samples for 4 samplesPerWave
Receive = repmat(struct('Apod', ones(1,numEl), ...
                        'startDepth', 0, ...
                        'endDepth', 0, ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'samplesPerWave', 4, ...
                        'mode', 0, ...
                        'InputFilter', repmat([0.0036,0.0127,0.0066,-0.0881,-0.2595,0.6484],Resource.Parameters.numRcvChannels,1), ...
                        'callMediaFunc', 1),1,(nr+1)*Resource.RcvBuffer(1).numFrames+Resource.RcvBuffer(2).numFrames+nr*Resource.RcvBuffer(3).numFrames); % call media function [edit]

% receive attributes for guidance b-mode (multiple acqs required for each frame)
for i = 1:Resource.RcvBuffer(1).numFrames
    k = nr*(i-1);
    for j = 1:nr 
        Receive(k+j).startDepth = SFormat(1).startDepth;
        Receive(k+j).endDepth = SFormat(1).startDepth + wlsPer128*ceil(maxAcqLengthBmode/wlsPer128);
        Receive(k+j).framenum = i;
        Receive(k+j).acqNum = j; 
        Receive(k+j).bufnum = 1; % added to specify b-mode buf [edit]
    end
end

% receive attributes for guidance pw imaging
for i = 1:Resource.RcvBuffer(1).numFrames
    k = i+nr*Resource.RcvBuffer(1).numFrames;
    Receive(k).startDepth = SFormat(2).startDepth;
    Receive(k).endDepth = SFormat(2).startDepth + wlsPer128*ceil(maxAcqLengthPW/wlsPer128);
    Receive(k).framenum = i;
    Receive(k).bufnum = 1; 
end

% receive attributes for SA (1 acq for multiple frames to be post-processed)
for i = 1:Resource.RcvBuffer(2).numFrames
    k = i+(nr+1)*Resource.RcvBuffer(1).numFrames;
    Receive(k).startDepth = SFormat(2).startDepth;
    Receive(k).endDepth = SFormat(2).startDepth + wlsPer128*ceil(maxAcqLengthPW/wlsPer128);
    Receive(k).framenum = i;      
    Receive(k).bufnum = 2; 
end

% receive attributes for b-mode save (multiple acqs required for each frame)
for i = 1:Resource.RcvBuffer(3).numFrames
    k = nr*(i-1)+(nr+1)*Resource.RcvBuffer(1).numFrames+Resource.RcvBuffer(2).numFrames;
    for j = 1:nr 
        Receive(k+j).startDepth = SFormat(1).startDepth;
        Receive(k+j).endDepth = SFormat(1).startDepth + wlsPer128*ceil(maxAcqLengthBmode/wlsPer128);
        Receive(k+j).framenum = i;
        Receive(k+j).acqNum = j; 
        Receive(k+j).bufnum = 3; % added to specify save b-mode buffer
    end
end

% image reconstruction for guidance b-mode 
Recon = repmat(struct('senscutoff', 0.6, ...
               'pdatanum', [], ...
               'newFrameTimeout',1000,...
               'IntBufDest', [0,0], ...
               'ImgBufDest', [], ...
               'RINums', []),1,2);

Recon(1).pdatanum = 1;
Recon(1).rcvBufFrame = -1;
Recon(1).ImgBufDest = [1,-1];
Recon(1).RINums = 1:nr;
           
% specific image reconstruction params for pw imaging
Recon(2).pdatanum = 2;
Recon(2).rcvBufFrame = -1; % changed from [2,-1]
Recon(2).ImgBufDest = [2,-1];
Recon(2).RINums = nr+1;

% define ReconInfo structures
ReconInfo = repmat(struct('mode', 0, ...  % replace data.
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'regionnum', 0), 1, nr+1);

% set specific ReconInfo attributes for phased array
for i = 1:nr
    ReconInfo(i).txnum = i;
    ReconInfo(i).rcvnum = i;
    ReconInfo(i).regionnum = i;
end

% set specific ReconInfo attributes for planewave
ReconInfo(nr+1).txnum = nr+1;
ReconInfo(nr+1).rxnum = 1;
ReconInfo(nr+1).regionnum = 0;

pers = 10; % persistence level 
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'norm',1,...        % normalization method(1 means fixed)
                         'pgain',0.5,...            % pgain is image processing gain
                         'persistMethod','none',...
                         'persistLevel',pers,...
                         'interp',1,...      % method of interpolation (1=4pt interp)
                         'compression',0.5,...      % X^0.5 normalized to output word size
                         'mappingMode','full',...
                         'display',1,...     % display image after processing
                         'displayWindow',1};
 
% save all b mode and SA data
Process(2).classname = 'External';
Process(2).method = 'P4_2_save_bmode';
Process(2).Parameters = {'srcbuffer','receive',...
    'srcbufnum',3,...
    'srcframenum',0,... % 0 = all frames, -1 = last frame
    'dstbuffer','none'};

Process(3).classname = 'External';
Process(3).method = 'P4_2_save_SA';
Process(3).Parameters = {'srcbuffer','receive',...
    'srcbufnum',2,...
    'srcframenum',0,... 
    'dstbuffer','none'};

Process(4).classname = 'External';
Process(4).method = 'P4_2_image_SA';
Process(4).Parameters = {'srcbuffer','receive',...
    'srcbufnum',2,...
    'srcframenum',-1,... 
    'dstbuffer','none'};

Process(5).classname = 'Image';
Process(5).method = 'imageDisplay';
Process(5).Parameters = {'imgbufnum',2,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',2,...    % number of PData structure to use
                         'norm',1,...        % normalization method(1 means fixed)
                         'pgain',1.0,...            % pgain is image processing gain
                         'persistMethod','simple',...
                         'persistLevel',30,...
                         'interp',1,...      % method of interpolation (1=4pt interp)
                         'compression',0.5,...      % X^0.5 normalized to output word size
                         'reject',2,...
                         'mappingMode','full',...
                         'display',1,...     % display image after processing
                         'displayWindow',1};

Process(6).classname = 'External';
Process(6).method = 'sweep_probe';
Process(6).Parameters = {'srcbuffer','none',...
    'srcframenum','none',... 
    'dstbuffer','none'};

Process(7).classname = 'External';
Process(7).method = 'P4_2_save_bmode_complete';
Process(7).Parameters = {'srcbuffer','receive',...
    'srcbufnum',3,...
    'srcframenum',0,... % 0 = all frames, -1 = last frame
    'dstbuffer','none'};

% bmode sequence control
SeqControl(1).command = 'timeToNextAcq';
SeqControl(1).argument = tPulse;  % 200 usec between ray lines
SeqControl(2).command = 'timeToNextAcq';
SeqControl(2).argument = round(1/PRF*1e6 - 128*tPulse);

% change to Profile 2 (SSA)
SeqControl(3).command = 'setTPCProfile';
SeqControl(3).condition = 'next';
SeqControl(3).argument = 2;

% change to Profile 1 (B-mode)
SeqControl(4).command = 'setTPCProfile';
SeqControl(4).condition = 'next';
SeqControl(4).argument = 1;

% time between SSA acquisitions
SeqControl(5).command = 'timeToNextAcq';
SeqControl(5).argument = ssaPRT;

% wait time for rs232 command [edit] disabled
SeqControl(6).command = 'noop';
SeqControl(6).argument = rs232wait;

% trigger for SSA acq
SeqControl(7).command = 'triggerOut'; 

SeqControl(8).command = 'returnToMatlab';

% time between pw imaging frames
SeqControl(9).command = 'timeToNextAcq';
SeqControl(9).argument = pwPRT;

nsc = 10;

% ****************** guidance b mode imaging routine **********************
n = 1;
bmode_image = n;
Event(n).info = 'Switch to b-mode acquire';
Event(n).tx = 0;         
Event(n).rcv = 0;       
Event(n).recon = 0;      
Event(n).process = 0;
Event(n).seqControl = 4; % switch to b-mode profile
n = n+1;

bmode_loop = n;
for i = 1:Resource.RcvBuffer(1).numFrames
    for j = 1:nr
        Event(n).info = 'Acquire b-mode ray line';
        Event(n).tx = j; 
        Event(n).rcv = nr*(i-1)+j;   
        Event(n).recon = 0;      
        Event(n).process = 0;    
        Event(n).seqControl = 1; % time to next acq (tPulse)
        n = n+1;
    end
    
    % Replace last events SeqControl for inter-frame timeToNextAcq.
    Event(n-1).seqControl = [2,nsc]; 
    SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
    nsc = nsc+1;

%     Event(n).info = 'Wait for transfer complete';
%     Event(n).tx = 0;         
%     Event(n).rcv = 0;        
%     Event(n).recon = 0;      
%     Event(n).process = 0;    
%     Event(n).seqControl = [nsc]; % wait for data to be transferred
%     SeqControl(nsc).command = 'waitForTransferComplete';
%     SeqControl(nsc).argument = nsc-1;
%     nsc = nsc+1;
%     n = n+1;    

    Event(n).info = 'Recon and process'; 
    Event(n).tx = 0;         
    Event(n).rcv = 0;        
    Event(n).recon = 1;      % image reconstruction
    Event(n).process = 1;    % image processing using process(1)  
    if floor(i/2) == i/2     % exit to Matlab every 2nd frame to check GUI
        Event(n).seqControl = 8;
    else
        Event(n).seqControl = 0;
    end
    n = n+1;
end

Event(n).info = 'Jump back to guidance bmode start';
Event(n).tx = 0;        
Event(n).rcv = 0;       
Event(n).recon = 0;     
Event(n).process = 0; 
Event(n).seqControl = nsc; % jump back to event 2 to continue b-mode imaging 
    SeqControl(nsc).command = 'jump';
    SeqControl(nsc).argument = bmode_loop;
    nsc = nsc+1;    
n = n+1;

% ****************** guidance planewave imaging routine *******************
pw_image = n;
Event(n).info = 'Switch to SA acquire';
Event(n).tx = 0;         
Event(n).rcv = 0;       
Event(n).recon = 0;      
Event(n).process = 0;
Event(n).seqControl = 3; % switch to SA profile
n = n+1;

pw_loop = n;
for i = 1:Resource.RcvBuffer(1).numFrames
    Event(n).info = 'Acquire PW image';
    Event(n).tx = nr+1; 
    Event(n).rcv = i+nr*Resource.RcvBuffer(1).numFrames;   
    Event(n).recon = 0;      
    Event(n).process = 0;    
    Event(n).seqControl = [9,nsc];
        SeqControl(nsc).command = 'transferToHost';
        nsc = nsc+1;
    n = n+1;

    Event(n).info = 'Recon and process'; 
    Event(n).tx = 0;         
    Event(n).rcv = 0;        
    Event(n).recon = 2;      % image reconstruction
    Event(n).process = 5;    % image processing using process(1)  
    if floor(i/5) == i/5
        Event(n).seqControl = 8;
    else
        Event(n).seqControl = 0;
    end
    n = n+1;
end

Event(n).info = 'Jump back to PW imagine start';
Event(n).tx = 0;        
Event(n).rcv = 0;       
Event(n).recon = 0;     
Event(n).process = 0; 
Event(n).seqControl = nsc;
    SeqControl(nsc).command = 'jump';
    SeqControl(nsc).argument = pw_loop;
    nsc = nsc+1;
n = n+1;

%*************************** b-mode acquisition ***************************
bmode_acq = n;
Event(n).info = 'Switch to b-mode acquire';
Event(n).tx = 0;         
Event(n).rcv = 0;       
Event(n).recon = 0;      
Event(n).process = 0;
Event(n).seqControl = 4; % switch to b-mode profile
n = n+1;

for i = 1:Resource.RcvBuffer(3).numFrames
    for j = 1:nr                 
        Event(n).info = 'Acquire ray line';
        Event(n).tx = j;         
        Event(n).rcv = nr*(i-1)+j+(nr+1)*Resource.RcvBuffer(1).numFrames+Resource.RcvBuffer(2).numFrames;        % use rcv for only single frame   
        Event(n).recon = 0;      
        Event(n).process = 0;    
        Event(n).seqControl = 1; % seqCntrl
        n = n+1;
    end

% Replace last events SeqControl (switch to SA) with transfer to host.
Event(n-1).seqControl = [2,nsc];
    SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
    SeqControl(nsc).condition = 'waitForProcessing'; % wait for the channel to be processed
    SeqControl(nsc).argument = nsc;
    nsc = nsc+1;

Event(n).info = 'Wait for transfer complete';
Event(n).tx = 0;         
Event(n).rcv = 0;        
Event(n).recon = 0;      
Event(n).process = 0;    
Event(n).seqControl = [nsc,nsc+1]; % wait for data to be transferred, mark processed
    SeqControl(nsc).command = 'waitForTransferComplete';
    SeqControl(nsc).argument = nsc-1;
    SeqControl(nsc+1).command = 'markTransferProcessed';
    SeqControl(nsc+1).argument = nsc-1;
    nsc = nsc+2;
n = n+1;

end

Event(n).info = 'Save b-mode image';
Event(n).tx = 0; 
Event(n).rcv = 0; 
Event(n).recon = 0; 
Event(n).process = 7; 
Event(n).seqControl = nsc;
    SeqControl(nsc).command = 'returnToMatlab';
    nsc = nsc+1;
n = n+1;

Event(n).info = 'Jump back to guidance b-mode'; % [edit] remove if B-mode needed within sequence
Event(n).tx = 0;        % no TX
Event(n).rcv = 0;       % no Rcv
Event(n).recon = 0;     % no Recon
Event(n).process = 0; 
Event(n).seqControl = nsc;
    SeqControl(nsc).command = 'jump';
    SeqControl(nsc).argument = bmode_image;
    nsc = nsc+1;
n = n+1;

%*********************** planewave SSA acquisition ************************
% 
SSA_acq = n;

% 1st bmode frame before sweep
Event(n).info = 'Switch to b-mode acquire';
Event(n).tx = 0;         
Event(n).rcv = 0;       
Event(n).recon = 0;      
Event(n).process = 0;
Event(n).seqControl = 4; % switch to b-mode profile
n = n+1;

for j = 1:nr                 
    Event(n).info = 'Acquire ray line';
    Event(n).tx = j;         
    Event(n).rcv = j+(nr+1)*Resource.RcvBuffer(1).numFrames+Resource.RcvBuffer(2).numFrames;        % use rcv for only single frame   
    Event(n).recon = 0;      
    Event(n).process = 0;    
    Event(n).seqControl = 1; % seqCntrl
    n = n+1;
end
Event(n-1).seqControl = [2,nsc];
    SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
    SeqControl(nsc).condition = 'waitForProcessing'; % wait for the channel to be processed
    SeqControl(nsc).argument = nsc;
    nsc = nsc+1;

Event(n).info = 'Wait for transfer complete';
Event(n).tx = 0;         
Event(n).rcv = 0;        
Event(n).recon = 0;      
Event(n).process = 0;    
Event(n).seqControl = [nsc,nsc+1]; % wait for data to be transferred, mark processed
    SeqControl(nsc).command = 'waitForTransferComplete';
    SeqControl(nsc).argument = nsc-1;
    SeqControl(nsc+1).command = 'markTransferProcessed';
    SeqControl(nsc+1).argument = nsc-1;
    nsc = nsc+2;
n = n+1;

% start sweep and take SSA frames

Event(n).info = 'Trigger sweep';
Event(n).tx = 0;         
Event(n).rcv = 0;       
Event(n).recon = 0;      
Event(n).process = 6;
Event(n).seqControl = nsc;
    SeqControl(nsc).command = 'returnToMatlab';
    nsc = nsc+1;
n = n+1;

Event(n).info = 'Switch to SA acquire';
Event(n).tx = 0;         
Event(n).rcv = 0;       
Event(n).recon = 0;      
Event(n).process = 0;
Event(n).seqControl = 3; % switch to SA profile
n = n+1;

% Event(n).info = 'Wait for turn table command'; % uncomment for rs232 wait
% Event(n).tx = 0;         
% Event(n).rcv = 0;       
% Event(n).recon = 0;      
% Event(n).process = 0;
% Event(n).seqControl = 6; 
% n = n+1;

SSA_loop = n;
for i = 1:Resource.RcvBuffer(2).numFrames
    Event(n).info = 'Send trigger out.';
    Event(n).tx = 0;
    Event(n).rcv = 0; 
    Event(n).recon = 0; 
    Event(n).process = 0; 
    Event(n).seqControl = 7; 
    n = n+1;
    
    Event(n).info = 'Acquire SA RF Data.';
    Event(n).tx = nr+1;
    Event(n).rcv = i+(nr+1)*Resource.RcvBuffer(1).numFrames; 
    Event(n).recon = 0; % no reconstruction
    Event(n).process = 0; % no processing
    Event(n).seqControl = [5,nsc]; % wait acq time transfer data to host
        SeqControl(nsc).command = 'transferToHost';
        nsc = nsc+1;
    n = n+1;
    
%     Event(n).info = 'debug';
%     Event(n).tx = 0;
%     Event(n).rcv = 0; 
%     Event(n).recon = 0; 
%     Event(n).process = 4; 
%     Event(n).seqControl = [8]; 
%     n = n+1;
end
    
Event(n).info = 'Wait and sync.';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = [nsc nsc+1 8];
    SeqControl(nsc).command = 'waitForTransferComplete';
    SeqControl(nsc).argument = nsc-1;
    SeqControl(nsc+1).command = 'sync';
    nsc = nsc+2;
n = n+1;

% 2nd bmode frame after sweep
Event(n).info = 'Switch to b-mode acquire';
Event(n).tx = 0;         
Event(n).rcv = 0;       
Event(n).recon = 0;      
Event(n).process = 0;
Event(n).seqControl = 4; % switch to b-mode profile
n = n+1;

for j = 1:nr                 
    Event(n).info = 'Acquire ray line';
    Event(n).tx = j;         
    Event(n).rcv = j+nr+(nr+1)*Resource.RcvBuffer(1).numFrames+Resource.RcvBuffer(2).numFrames;        % use rcv for only single frame   
    Event(n).recon = 0;      
    Event(n).process = 0;    
    Event(n).seqControl = 1; % seqCntrl
    n = n+1;
end
Event(n-1).seqControl = [2,nsc];
    SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
    SeqControl(nsc).condition = 'waitForProcessing'; % wait for the channel to be processed
    SeqControl(nsc).argument = nsc;
    nsc = nsc+1;

Event(n).info = 'Wait for transfer complete';
Event(n).tx = 0;         
Event(n).rcv = 0;        
Event(n).recon = 0;      
Event(n).process = 0;    
Event(n).seqControl = [nsc,nsc+1]; % wait for data to be transferred, mark processed
    SeqControl(nsc).command = 'waitForTransferComplete';
    SeqControl(nsc).argument = nsc-1;
    SeqControl(nsc+1).command = 'markTransferProcessed';
    SeqControl(nsc+1).argument = nsc-1;
    nsc = nsc+2;
n = n+1;

Event(n).info = 'Save b-mode image';
Event(n).tx = 0; 
Event(n).rcv = 0; 
Event(n).recon = 0; 
Event(n).process = 2; 
Event(n).seqControl = nsc;
    SeqControl(nsc).command = 'returnToMatlab';
    nsc = nsc+1;
n = n+1;
% end of bmode saving

Event(n).info = 'Save SA frames and return to imaging'; 
Event(n).tx = 0;         
Event(n).rcv = 0;        
Event(n).recon = 0;      
Event(n).process = 3;
Event(n).seqControl = nsc;
    SeqControl(nsc).command = 'returnToMatlab';
    nsc = nsc+1;
n = n+1;

Event(n).info = 'Jump back to guidance b-mode'; % [edit] remove if B-mode needed within sequence
Event(n).tx = 0;        % no TX
Event(n).rcv = 0;       % no Rcv
Event(n).recon = 0;     % no Recon
Event(n).process = 0; 
Event(n).seqControl = nsc;
    SeqControl(nsc).command = 'jump';
    SeqControl(nsc).argument = bmode_image;
    nsc = nsc+1;
n = n+1;

%*********************** steered planewave SSA acquisition ****************
steerSSA_acq = n;

% 1st bmode frame before sweep
Event(n).info = 'Switch to b-mode acquire';
Event(n).tx = 0;         
Event(n).rcv = 0;       
Event(n).recon = 0;      
Event(n).process = 0;
Event(n).seqControl = 4; % switch to b-mode profile
n = n+1;

for j = 1:nr                 
    Event(n).info = 'Acquire ray line';
    Event(n).tx = j;         
    Event(n).rcv = j+(nr+1)*Resource.RcvBuffer(1).numFrames+Resource.RcvBuffer(2).numFrames;        % use rcv for only single frame   
    Event(n).recon = 0;      
    Event(n).process = 0;    
    Event(n).seqControl = 1; % seqCntrl
    n = n+1;
end
Event(n-1).seqControl = [2,nsc];
    SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
    SeqControl(nsc).condition = 'waitForProcessing'; % wait for the channel to be processed
    SeqControl(nsc).argument = nsc;
    nsc = nsc+1;

Event(n).info = 'Wait for transfer complete';
Event(n).tx = 0;         
Event(n).rcv = 0;        
Event(n).recon = 0;      
Event(n).process = 0;    
Event(n).seqControl = [nsc,nsc+1]; % wait for data to be transferred, mark processed
    SeqControl(nsc).command = 'waitForTransferComplete';
    SeqControl(nsc).argument = nsc-1;
    SeqControl(nsc+1).command = 'markTransferProcessed';
    SeqControl(nsc+1).argument = nsc-1;
    nsc = nsc+2;
n = n+1;

% start sweep and take SSA frames
Event(n).info = 'Trigger sweep';
Event(n).tx = 0;         
Event(n).rcv = 0;       
Event(n).recon = 0;      
Event(n).process = 6;
Event(n).seqControl = nsc;
    SeqControl(nsc).command = 'returnToMatlab';
    nsc = nsc+1;
n = n+1;

Event(n).info = 'Switch to SA acquire';
Event(n).tx = 0;         
Event(n).rcv = 0;       
Event(n).recon = 0;      
Event(n).process = 0;
Event(n).seqControl = 3; % switch to SA profile
n = n+1;

% Event(n).info = 'Wait for turn table command'; % uncomment for rs232 wait
% Event(n).tx = 0;         
% Event(n).rcv = 0;       
% Event(n).recon = 0;      
% Event(n).process = 0;
% Event(n).seqControl = 6; 
% n = n+1;

steerSSA_loop = n;
j=1;
steerAngles = zeros(1,Resource.RcvBuffer(2).numFrames);
for i = 1:Resource.RcvBuffer(2).numFrames
    Event(n).info = 'Send trigger out.';
    Event(n).tx = 0;
    Event(n).rcv = 0; 
    Event(n).recon = 0; 
    Event(n).process = 0; 
    Event(n).seqControl = 7;
    n = n+1;
    
    Event(n).info = 'Acquire SA RF Data.';
    Event(n).tx = nr+1+j;
    Event(n).rcv = i+(nr+1)*Resource.RcvBuffer(1).numFrames; 
    Event(n).recon = 0; 
    Event(n).process = 0; 
    Event(n).seqControl = [5,nsc]; % wait acq time transfer data to host
        SeqControl(nsc).command = 'transferToHost';
        nsc = nsc+1;
    n = n+1;
    
    steerAngles(i) = pwAngles(j); % save angles of scan
    j = j+1;
    if j > na; j = 1; end
end

Event(n).info = 'Wait and sync.';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = [nsc nsc+1 8];
    SeqControl(nsc).command = 'waitForTransferComplete';
    SeqControl(nsc).argument = nsc-1;
    SeqControl(nsc+1).command = 'sync';
    nsc = nsc+2;
n = n+1;

% 2nd bmode frame after sweep
Event(n).info = 'Switch to b-mode acquire';
Event(n).tx = 0;         
Event(n).rcv = 0;       
Event(n).recon = 0;      
Event(n).process = 0;
Event(n).seqControl = 4; % switch to b-mode profile
n = n+1;

for j = 1:nr                 
    Event(n).info = 'Acquire ray line';
    Event(n).tx = j;         
    Event(n).rcv = j+nr+(nr+1)*Resource.RcvBuffer(1).numFrames+Resource.RcvBuffer(2).numFrames;        % use rcv for only single frame   
    Event(n).recon = 0;      
    Event(n).process = 0;    
    Event(n).seqControl = 1; % seqCntrl
    n = n+1;
end
Event(n-1).seqControl = [2,nsc];
    SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
    SeqControl(nsc).condition = 'waitForProcessing'; % wait for the channel to be processed
    SeqControl(nsc).argument = nsc;
    nsc = nsc+1;

Event(n).info = 'Wait for transfer complete';
Event(n).tx = 0;         
Event(n).rcv = 0;        
Event(n).recon = 0;      
Event(n).process = 0;    
Event(n).seqControl = [nsc,nsc+1]; % wait for data to be transferred, mark processed
    SeqControl(nsc).command = 'waitForTransferComplete';
    SeqControl(nsc).argument = nsc-1;
    SeqControl(nsc+1).command = 'markTransferProcessed';
    SeqControl(nsc+1).argument = nsc-1;
    nsc = nsc+2;
n = n+1;

Event(n).info = 'Save b-mode image';
Event(n).tx = 0; 
Event(n).rcv = 0; 
Event(n).recon = 0; 
Event(n).process = 2; 
Event(n).seqControl = nsc;
    SeqControl(nsc).command = 'returnToMatlab';
    nsc = nsc+1;
n = n+1;

Event(n).info = 'Save SA frames and return to imaging'; 
Event(n).tx = 0;         
Event(n).rcv = 0;        
Event(n).recon = 0;      
Event(n).process = 3;
Event(n).seqControl = nsc;
    SeqControl(nsc).command = 'returnToMatlab';
    nsc = nsc+1;
n = n+1;

Event(n).info = 'Jump back to guidance b-mode'; % [edit] remove if B-mode needed within sequence
Event(n).tx = 0;        % no TX
Event(n).rcv = 0;       % no Rcv
Event(n).recon = 0;     % no Recon
Event(n).process = 0; 
Event(n).seqControl = nsc;
    SeqControl(nsc).command = 'jump';
    SeqControl(nsc).argument = bmode_image;
    nsc = nsc+1;
n = n+1;

%*********************** diverging SSA acquisition ************************


%************************ GUI elements ************************************

UI(1).Control = {'UserC3','Style','VsPushButton','Tag','scanSSA','Label','SSA Scan'};
UI(1).Callback = text2cell('%CB_SSA');

UI(2).Control = {'UserB5','Style','VsPushButton','Tag','testSweep','Label','Test Sweep'};
UI(2).Callback = text2cell('%CB#2');

UI(3).Control = {'UserB2','Style','VsPushButton','Tag','closePort','Label','Reset RS232'};
UI(3).Callback = text2cell('%CB#3');

UI(4).Control = {'UserB1','Style','VsPushButton','Tag','abortSweep','Label','Abort Sweep'};
UI(4).Callback = text2cell('%CB#4');

% UI(5).Control = {'UserB5','Style','VsPushButton','Tag','stepCW','Label','Step CW'};
% UI(5).Callback = text2cell('%CB#5');
% 
% UI(6).Control = {'UserB6','Style','VsPushButton','Tag','stepCCW','Label','Step CCW'};
% UI(6).Callback = text2cell('%CB#6');
% 
% UI(7).Control = {'UserB7','Style','VsPushButton','Tag','setOrigin','Label','Set Origin'};
% UI(7).Callback = text2cell('%CB#7');

UI(5).Control = {'UserC1','Style','VsPushButton','Tag','saveBmode','Label','Save B-mode'};
UI(5).Callback = text2cell('%CB#8');
    
UI(6).Control = {'UserB4','Style','VsPushButton','Tag','mvCenter','Label','Center Probe'};
UI(6).Callback = text2cell('%CB#9');

UI(7).Control = {'UserB3','Style','VsPushButton','Tag','pwImage','Label','Phased/Plane'};
UI(7).Callback = text2cell('%CB#10');

UI(8).Control = {'UserC2','Style','VsPushButton','Tag','scanaSSA','Label','aSSA Scan'};
UI(8).Callback = text2cell('%CB_aSSA');

%************************ output VSX **************************************
scriptName = 'P4-2_planewaveSSA';
disp(['filename =''' scriptName ''';VSX'])
save(scriptName);
return 

%CB_SSA
startSACallback(hObject,eventdata)
return
%CB_SSA

%CB_aSSA
startAngledSACallback(hObject,eventdata)
return
%CB_aSSA

%CB#2
testSweepCallback(hObject,eventdata)
return
%CB#2

%CB#3
closeSerialCallback(hObject,eventdata)
return
%CB#3

%CB#4
abortSweepCallback(hObject,eventdata)
return
%CB#4

%CB#5
stepCWCallback(hObject,eventdata)
return
%CB#5

%CB#6
stepCCWCallback(hObject,eventdata)
return
%CB#6

%CB#7
setOriginCallback(hObject,eventdata)
return
%CB#7

%CB#8
saveBmodeCallback(hObject,eventdata)
return
%CB#8

%CB#9
mvCenterCallback(hObject,eventdata)
return
%CB#9

%CB#10
imageToggleCallback(hObject,eventdata)
return
%CB#10