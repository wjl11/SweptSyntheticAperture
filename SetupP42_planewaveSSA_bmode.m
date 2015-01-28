clear all; close all;
% General imaging parameters
maxVoltage = 50;
numEl = 64; % number of physical elements (64 elements in the P4-2, only half channels of vsx used)

% SA parameters
SAPRT = 50e3; % time between SA planewave acquisitions
nFramesSweep = 50;
rowsPerFrameSA = 4096; 

% B-mode parameters
nr = 128; % no of ray lines in b-mode scan
rowsPerFrameBmode = nr*4096; % rows of RF data stored **** changed from nr*4096 [edit]
PRF = 6.67; % time between b-mode frames
foc_mm = 60;
endDepth_mm = 80;
txfnum = 2;
tPulse = 200; % time between phased array acqs [us]
c = 1540; 

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
% Media.function = 'movePoints';

% specify system parameters
Resource.Parameters.numTransmit = 128;
Resource.Parameters.numRcvChannels = 128; % no of receive channels for 4 board VDAS unit
% Resource.Parameters.connector = 1; % tdr connector to use for 4 board VDAS unit
Resource.Parameters.speedOfSound = c;
Resource.Parameters.simulateMode = 1;
if Resource.Parameters.simulateMode ~= 0
    disp('Run in simulation mode.')
else
    disp('Attach P4-2 transducer.')
end

% computeTrans function used to define the transducer struct
freqMHz = 2.50;
Trans.name = 'P4-2';
Trans.frequency = freqMHz; % not needed if using default center freq -
% assume center frequency
Trans = computeTrans(Trans);
Trans.maxHighVoltage = maxVoltage;

% Set up SFormat structure array.
aperture = numEl*Trans.spacing;
SFormat.transducer = 'P4-2';
SFormat.scanFormat = 'VAPX'; 
SFormat.theta = -pi/4;
SFormat.radius = (aperture/2)/tan(-SFormat.theta); % dist. to virt. apex
SFormat.numRays = nr;      % no. of Rays
SFormat.FirstRayLoc = Trans.ElementPos(1,1:3);   % x,y,z
SFormat.rayDelta = 2*(-SFormat.theta)/(nr-1);  % spacing in radians(sector) or dist. between rays
SFormat.startDepth = 0;
SFormat.endDepth = round(endDepth_mm/1000/(c/(Trans.frequency*1e6))); % Acquisition depth in wavelengths

% Set up PData structure. Pixels representing phased array image
PData.sFormat = 1;
PData.pdelta = 0.875;
PData.Size(1) = 10 + ceil((SFormat.endDepth-SFormat.startDepth)/PData.pdelta);
PData.Size(2) = 10 + ceil(2*(SFormat.endDepth + SFormat.radius)*sin(-SFormat.theta)/PData.pdelta);
PData.Size(3) = 1;
PData.Origin = [-(PData.Size(2)/2)*PData.pdelta,0,0];

% b-mode resource and image buffer and image display window
bmode_nframes = 2; %**** may need to change [edit]
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = rowsPerFrameBmode; 
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = bmode_nframes; % [edit]

Resource.ImageBuffer(1).datatype = 'double';
Resource.ImageBuffer(1).rowsPerFrame = PData.Size(1); % display all pixels
Resource.ImageBuffer(1).colsPerFrame = PData.Size(2);
Resource.ImageBuffer(1).numFrames = bmode_nframes;

Resource.DisplayWindow(1).Title = 'Image Display';
Resource.DisplayWindow(1).pdelta = 0.3;
Resource.DisplayWindow(1).Position = [200,200, ...
    ceil(PData.Size(2)*(PData.pdelta/Resource.DisplayWindow(1).pdelta)), ...
    ceil(PData.Size(1)*PData.pdelta/Resource.DisplayWindow(1).pdelta)];
Resource.DisplayWindow(1).ReferencePt = [PData.Origin(1),PData.Origin(3)];
Resource.DisplayWindow(1).Colormap = gray(256);

% planewave SA resource buffer
Resource.RcvBuffer(2).datatype = 'int16';
Resource.RcvBuffer(2).rowsPerFrame = rowsPerFrameSA; % number of samples of RF data per channel, default sample rate for A/D 
% is 4x the center frequency of the probe
Resource.RcvBuffer(2).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(2).numFrames = nFramesSweep; % may need to buffer [edit]

% transmit waveform for both SA and bmode
TW(1).type = 'parametric';
TW(1).Parameters = [36,17,2,1];   % A, B, C, D
% A,B,C,D in terms of master clock (180 MHz T = 5.556 ns)
% A = clocks in half cycle based on f0
% B = clocks during which transmit drivers are active
% C = number of half cycles in transmit waveform (2 = single cycle)
% D = polarity of first half cycle (1 = positive)

% tx definitions for b-mode imaging
tx_foc = round(foc_mm/1000/(c/(Trans.frequency*1e6))); %focus in wl
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'focus', 0, ... % in wavelengths; 
                   'Steer', [0.0,0.0], ...
                   'Apod', ones(1,numEl), ...  
                   'Delay', zeros(1,numEl)), 1, nr+1); 

% % tx definitions for planewave SA imaging               
% TX(nr+1) = struct('waveform', 1, ...
%                    'focus', 0, ... % in wavelengths; 
%                    'Apod', ones(1,numEl), ...  
%                    'Delay', computeTXDelays(TX(nr+1)));

% specify unique phased array tx attributes 
angles = SFormat.theta:SFormat.rayDelta:(SFormat.theta + (nr-1)*SFormat.rayDelta);
origin = SFormat.radius*tan(angles);
for n = 1:nr
    TX(n).Origin = [origin(n),0.0,0.0];
    TX(n).focus = tx_foc;
    TX(n).Steer = [angles(n),0.0];
    TX(n).Delay = computeTXDelays(TX(n));
end
clear origin angles

TX(nr+1).Delay = computeTXDelays(TX(nr+1));

% specify transmit power controllers (optional)
TPC(1).name = 'B-mode';
TPC(1).maxHighVoltage = maxVoltage;
TPC(2).name = 'Synthetic Aperture';
TPC(2).maxHighVoltage = maxVoltage;

% specify TGC
TGC.CntrlPts = [450,550,650,710,770,830,890,950]; % [First one for L12_5 is 400]
TGC.rangeMax = SFormat.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% rcv defintions for b-mode 
maxAcqLength = sqrt(aperture^2 + SFormat.endDepth^2 - 2*aperture*SFormat.endDepth*cos(SFormat.theta-pi/2)) - SFormat.startDepth;
wlsPer128 = 128/(2*4); % wavelengths in 128 samples for 4 samplesPerWave
Receive = repmat(struct('Apod', ones(1,numEl), ...
                        'startDepth', SFormat.startDepth, ...
                        'endDepth', SFormat.startDepth + wlsPer128*ceil(maxAcqLength/wlsPer128), ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'samplesPerWave', 4, ...
                        'mode', 0, ...
                        'InputFilter', repmat([0.0036,0.0127,0.0066,-0.0881,-0.2595,0.6484],Resource.Parameters.numRcvChannels,1), ...
                        'callMediaFunc', 1),1,nr*Resource.RcvBuffer(1).numFrames+Resource.RcvBuffer(2).numFrames); % call media function [edit]

% receive attributes for guidance b-mode (multiple acqs required for each frame)
for i = 1:Resource.RcvBuffer(1).numFrames
    k = nr*(i-1);
    for j = 1:nr % changed from SFormat.numRays [edit]
        Receive(k+j).framenum = i;
        Receive(k+j).acqNum = j; 
        Receive(k+j).bufnum = 1; % added to specify b-mode buf [edit]
    end
end

% receive attributes for SA (1 acq for multiple frames to be post-processed)
for i = 1:Resource.RcvBuffer(2).numFrames
    k = i+nr*Resource.RcvBuffer(1).numFrames;
    Receive(k).framenum = i;
    Receive(k).acqNum = 1;       
    Receive(k).bufnum = 2; 
end

% image reconstruction for guidance b-mode
Recon = struct('senscutoff', 0.6, ...
               'pdatanum', 1, ...
               'rcvBufFrame',[1,-1], ...
               'IntBufDest', [0,0], ...
               'ImgBufDest', [1,-1], ...
               'RINums', 1:nr);

% define ReconInfo structures.
ReconInfo = repmat(struct('mode', 0, ...  % replace data.
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'regionnum', 0), 1, nr);

% set specific ReconInfo attributes.
for i = 1:nr
    ReconInfo(i).txnum = i;
    ReconInfo(i).rcvnum = i;
    ReconInfo(i).regionnum = i;
end

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
    'srcbufnum',1,...
    'srcframenum',0,...
    'dstbuffer','none'};

Process(3).classname = 'External';
Process(3).method = 'P4_2_save_SA';
Process(3).Parameters = {'srcbuffer','receive',...
    'srcbufnum',2,...
    'srcframenum',0,...
    'dstbuffer','none'};

% potential method for triggering frame acq with rs232 pulse [edit]
% Process(4).classname = 'External';
% Process(4).method = 'wait_rs232';
% Process(4).Parameters = {'srcbuffer','none',...
%     'srcbufnum',

n_process = 3;
% process for each SA frame
for i = 1:Resource.RcvBuffer(2).numFrames
    Process(n_process+i).classname = 'External';
    Process(n_process+i).method = 'P4_2_process_SA';
    Process(n_process+i).Parameters = {'srcbuffer','receive',... % name of buffer to process
        'srcbufnum',2,... 
        'srcframenum',i,...   
        'dstbuffer','none'}; 
end

% Specify SeqControl structure arrays.
%  - Jump back to start.
SeqControl(1).command = 'timeToNextAcq';
SeqControl(1).argument = tPulse;  % 200 usec between ray lines
SeqControl(2).command = 'timeToNextAcq';
SeqControl(2).argument = round(1/PRF*1e6 - 128*tPulse);
SeqControl(3).command = 'jump';
SeqControl(3).argument = 2;

% -- Change to Profile 2 (SA)
SeqControl(4).command = 'setTPCProfile';
SeqControl(4).condition = 'next';
SeqControl(4).argument = 2;

% -- Change to Profile 1 (B-mode)
SeqControl(5).command = 'setTPCProfile';
SeqControl(5).condition = 'next';
SeqControl(5).argument = 1;

% -- Time between planewave acquisitions
SeqControl(6).command = 'timeToNextAcq';
SeqControl(6).argument = SAPRT;
% SeqControl(6).argument = round(1/SAPRF*1e6 - 128*tPulse); % ****  SAPRF not valid here [edit] 
nsc = 7;

% ****************** guidance b mode imaging routine **********************
n = 1;
bmode_start = n;
Event(n).info = 'Switch to b-mode acquire';
Event(n).tx = 0;         
Event(n).rcv = 0;       
Event(n).recon = 0;      
Event(n).process = 0;
Event(n).seqControl = 5; % switch to b-mode profile
n = n+1;

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

    Event(n).info = 'Wait for transfer complete';
    Event(n).tx = 0;         
    Event(n).rcv = 0;        
    Event(n).recon = 0;      
    Event(n).process = 0;    
    Event(n).seqControl = [nsc]; % wait for data to be transferred
    SeqControl(nsc).command = 'waitForTransferComplete';
    SeqControl(nsc).argument = nsc-1;
    nsc = nsc+1;
    n = n+1;    

    Event(n).info = 'Recon and process'; 
    Event(n).tx = 0;         
    Event(n).rcv = 0;        
    Event(n).recon = 1;      % image reconstruction
    Event(n).process = 1;    % image processing using process(1)
    Event(n).seqControl = 0;    
    if floor(i/2) == i/2     % exit to Matlab every 2nd frame to check GUI
        Event(n).seqControl = nsc;
        SeqControl(nsc).command = 'returnToMatlab';
        nsc = nsc+1;
    end
    n = n+1;
end

Event(n).info = 'Jump back';
Event(n).tx = 0;        
Event(n).rcv = 0;       
Event(n).recon = 0;     
Event(n).process = 0; 
Event(n).seqControl = 3; % jump back to event 2 to continue b-mode imaging 
n = n+1;

bmode_end = n;

%*************************** save b-mode frame ****************************
for j = 1:nr                 
    Event(n).info = 'Acquire ray line';
    Event(n).tx = j;         
    Event(n).rcv = j;        % use rcv for only single frame   
    Event(n).recon = 0;      
    Event(n).process = 0;    
    Event(n).seqControl = 1; % seqCntrl
    n = n+1;
end

% Replace last events SeqControl (switch to SA) with transfer to host.
Event(n-1).seqControl = [2,4,nsc];
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
Event(n).process = 2; % process
Event(n).seqControl = nsc;
SeqControl(nsc).command = 'returnToMatlab';
nsc = nsc+1;
n = n+1;

%*********************** plane wave SA acquisition ************************
SA_start = n;
for i = 1:Resource.RcvBuffer(2).numFrames
    Event(n).info = 'Acquire SA RF Data.';
    Event(n).tx = nr+1;
    Event(n).rcv = i+nr*Resource.RcvBuffer(1).numFrames; 
    Event(n).recon = 0; % no reconstruction
    Event(n).process = 0; % no processing
    Event(n).seqControl = [6,nsc]; % wait acq time transfer data to host
        SeqControl(nsc).command = 'transferToHost';
        nsc = nsc+1;
    n = n+1;
    
    Event(n).info = 'Call external Processing function';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = n_process+i;
    Event(n).seqControl = [nsc nsc+1];
        SeqControl(nsc).command = 'waitForTransferComplete';
        SeqControl(nsc).argument = nsc-1;
        SeqControl(nsc+1).command = 'sync';
        nsc = nsc+2;
    n = n+1;
end

Event(n).info = 'Save SA frames'; 
Event(n).tx = 0;         
Event(n).rcv = 0;        
Event(n).recon = 0;      
Event(n).process = 3;
Event(n).seqControl = nsc;
SeqControl(nsc).command = 'returnToMatlab';
nsc = nsc+1;
n = n+1;

Event(n).info = 'Jump back to SA';
Event(n).tx = 0;        % no TX
Event(n).rcv = 0;       % no Rcv
Event(n).recon = 0;     % no Recon
Event(n).process = 0; 
Event(n).seqControl = nsc;
    SeqControl(nsc).command = 'jump';
    SeqControl(nsc).argument = SA_start;
    nsc = nsc+1;
n = n+1;

%************************ GUI elements ************************************

UI(1).Control = {'UserB1','Style','VsToggleButton','Tag','toggleSA','Label','Toggle SA'};
UI(1).Callback = text2cell('%CB#1');

%************************ output VSX **************************************
scriptName = 'P4-2_planewaveSSA';
display(['filename =''' scriptName ''';VSX'])
save(scriptName);

return 
%CB#1
toggleSACallback(hObject,eventdata)
return
%CB#1
