% ***** P4-2 Swept Synthetic Aperture Image Acquisition Routine *****
% Version 12.4 (debug in progress)
%
% Overhaul from original version 10.1 framework 
%     Git SHA v10.1: c26a26bf91c412d7f06841dca3f03c2e8a702fe9 (working)
%     Git SHA v11: d3349b3c815baa076c7b8e26f63e52a84b85e509 
%
% SA overhaul from version 12.4
%     Git SHA v12.4: 45b02aea59d05754a60e41e0eeed760918c62fe3 (working)
%
%
% Latest revision 03/01/15 (debugging in process) 
% -Debug SA acquisition routine (time zero)
% -Optimize save routine using fwritef
%
%
% Revision 12 02/17/15 (debugging in progress) - by Will Long
% -Modifications to improve code readibility 
% -Revisions to reduce redundancy of turn table rs232 communication
% -Implementation of manual routine
% -Focused Tx SSA removed from script
% -Eliminated B-mode acquisitions before and after SSA acquisition
% -Revised Tx and Rcv structure indexing
% -Implemented full SA acquisition
% -v7.3 saving implemented

vsx_path = '/home/wjl11/Matlab Simulator/';
addpath(genpath(vsx_path));

clear all
%%%%%%%%%%%%%%%%%%%%%%
% Imaging Case Setup %
%%%%%%%%%%%%%%%%%%%%%%
SETUP.simToggle = 0;            % 1 [simulation], 0 [probe connected]
SETUP.rs232Toggle = 1;          % 1 [table com on], 0 [table com off]
SETUP.scanType = input('Scan type [manual/turntable]: ','s');   
                                % 'manual' [manual scanning]
                                % 'turntable' [turn table scanning]
                            
%%%%%%%%%%%%%%%%%%%%%%
% General Parameters %
%%%%%%%%%%%%%%%%%%%%%%
maxVoltage = 50;            % max voltage used for transmit
numEl = 64;                 % number of physical elements (1/2 tot channel)
c = 1540;                   % speed of sound

%%%%%%%%%%%%%%%%%%
% SSA Parameters %
%%%%%%%%%%%%%%%%%%
SSA.PRT = 1e3;              % time between SSA acquisitions [us]
SSA.rowsPerFrame = 4096;
SSA.endDepthMM = 150;
SSA.startDepthMM = [];
SSA.nFrames = 1000;         % WARNING: value overridden in turntable mode
SSA.frameBuffer = 550;      % extra frames to pad turntable acquisition

% STEERED:
SSA.pwAngles = -10:10:10;
SSA.nAngles = length(SSA.pwAngles);

% FOCUSED:
SSA.focusMM = 90;           % (not implemented in version 11)

%%%%%%%%%%%%%%%%%
% SA Parameters %
%%%%%%%%%%%%%%%%%
SA.txFnum = 0;
SA.numEl = 1;
SA.PRF = 10;
SA.nRay = numEl-SA.numEl+1;
SA.rowsPerFrame = SA.nRay*4096;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Phased B-mode Parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PHASED_B.nRay = 128;
PHASED_B.rowsPerFrame = PHASED_B.nRay*4096;
PHASED_B.PRF = 23;          % pulse repetition frequency [hz]
PHASED_B.focusMM = 100;
PHASED_B.endDepthMM = 150;  
PHASED_B.tPulse = 250;      % time between phased tx [us]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Planewave B-mode parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PW_B.endDepthMM = SSA.endDepthMM;
PW_B.PRT = 1e4;             % time between pw tx [us]

IM_STATE = 'pa';            % default guidance imaging mode
                            % pa [phased], pw [planewave]

%********************* END OF GENERAL PARAMETER SETUP *********************

% DYNAMIC FLAGS
SSA_TYPE = [];

% BUFFER SIZES
imgBuffer = 10;
bmodeSaveFrames = 2;
saSaveFrames = 1;

% TURN TABLE SWEEP LIMITS AND SETTINGS
TABLE.sweep_range = [0.0 45.0];
TABLE.sweep_limits = [0.0 45.0]; % DO NOT ALTER
TABLE.scan_velocity = 5.0;
TABLE.norm_velocity = 3.0;
TABLE.step = 1.0;
TABLE.acc_fnc = 3;         % 0 [impulse], 1 [flat], 2 [ramp], 3 [sinusoid] 
if TABLE.sweep_range(1) < TABLE.sweep_limits(1) || ...
        TABLE.sweep_range(2) > TABLE.sweep_limits(2) || ...
        TABLE.sweep_range(1) > TABLE.sweep_range(2) ||...
        TABLE.sweep_limits(1) > TABLE.sweep_limits(2) 
    error('Error in TDR sweep bounds.')
end

switch SETUP.scanType
    case 'manual'

    case 'turntable'
        SSA.nFrames = ...
            round(((TABLE.sweep_range(2)-TABLE.sweep_range(1))/360.0)...
            /(TABLE.scan_velocity/60)/(SSA.PRT*1e-6))+SSA.frameBuffer;
end

% CONSOLE OUTPUTS
disp(['Plane wave SA frame rate: ' num2str(1/(SSA.PRT*1e-3)) ' kHz']);
disp(['SSA frames per acquisition: ' num2str(SSA.nFrames)])
% disp(['Fundamental frame rate: ' num2str(1/(SSA.endDepthMM*2/c)) ' kHz']);

% DATA LABEL
saveLabel = input('Data file unique identifier: ','s');

%********************* BEGIN VSX PARAMETER DEFINITION *********************

% SIMULATED MEDIUM DEFINITION
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

% SYSTEM PARAMETERS
Resource.Parameters.numTransmit = 128;
Resource.Parameters.numRcvChannels = 256;       % 256 for Vantage system
Resource.Parameters.speedOfSound = c;
Resource.Parameters.simulateMode = SETUP.simToggle;

if Resource.Parameters.simulateMode ~= 0, disp('Simulation mode: ON')
else warning('Simulation mode: OFF'), end

if SETUP.rs232Toggle ~= 1, disp('Serial communication: OFF')
else warning('Serial communication: ON'), end

% DEFINE TRANSDUCER
freqMHz = 2.50;
Trans.name = 'P4-2';
Trans.frequency = freqMHz;          % not needed if using default f0
% Trans.units = 'mm'; 
Trans.units = 'wavelengths'; 
Trans = computeTrans(Trans);
Trans.maxHighVoltage = maxVoltage;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Guidance Imaging Display Structures %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SFORMAT FOR DISPLAY (PHASED ARRAY)
aperture = numEl*Trans.spacing;
SFormat(1).transducer = 'P4-2';
SFormat(1).scanFormat = 'VAPX'; 
SFormat(1).theta = -pi/8;
SFormat(1).radius = (aperture/2)/tan(-SFormat(1).theta);                            % distance to virtual apex
SFormat(1).numRays = PHASED_B.nRay;
SFormat(1).FirstRayLoc = Trans.ElementPos(1,1:3);                                   % x,y,z
SFormat(1).rayDelta = 2*(-SFormat(1).theta)/(PHASED_B.nRay-1);                      % spacing in radians(sector) between rays
SFormat(1).startDepth = 0;
SFormat(1).endDepth = round(PHASED_B.endDepthMM/1000/(c/(Trans.frequency*1e6)));    % acquisition depth [wavelengths]

% SFORMAT FOR DISPLAY (PLANEWAVE)
SFormat(2).transducer = 'P4-2';
SFormat(2).scanFormat = 'RLIN'; 
SFormat(2).theta = 0;
SFormat(2).radius = 0; 
SFormat(2).numRays = 1;
SFormat(2).FirstRayLoc = [0,0,0];                                                   % x,y,z
SFormat(2).rayDelta = aperture;                                                     % distance between rays given by aperture dim
SFormat(2).startDepth = 0;
SFormat(2).endDepth = round(PW_B.endDepthMM/1000/(c/(Trans.frequency*1e6)));

% PIXEL DATA (PHASED ARRAY)
PData(1).sFormat = 1;
PData(1).pdelta = 0.875;
PData(1).Size(1) = 10 + ceil((SFormat(1).endDepth-SFormat(1).startDepth)/PData(1).pdelta);
PData(1).Size(2) = 10 + ceil(2*(SFormat(1).endDepth + SFormat(1).radius)*sin(-SFormat(1).theta)/PData(1).pdelta);
PData(1).Size(3) = 1;
PData(1).Origin = [-(PData(1).Size(2)/2)*PData(1).pdelta,0,0];

% PIXEL DATA (PLANEWAVE)
PData(2).sFormat = 2;
PData(2).pdeltaX = Trans.spacing;
PData(2).pdeltaZ = 0.5;                                                             
PData(2).Size(1) = ceil((SFormat(2).endDepth-SFormat(2).startDepth)/PData(2).pdeltaZ);
PData(2).Size(2) = ceil(aperture/PData(2).pdeltaX);
PData(2).Size(3) = 1;
PData(2).Origin = [-Trans.spacing*(numEl-1)/2,0,SFormat(2).startDepth];

% NOTE: smaller values for pdetla give higher resolution however result in 
% longer image reconstruction time (decrease FR)

% RESOURCE BUFFER 1 (DISPLAY AND RECON USE ONLY)
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = PHASED_B.rowsPerFrame; 
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = imgBuffer;

% IMAGE BUFFER (PHASED ARRAY)
Resource.ImageBuffer(1).datatype = 'double';
Resource.ImageBuffer(1).rowsPerFrame = PData(1).Size(1);
Resource.ImageBuffer(1).colsPerFrame = PData(1).Size(2);
Resource.ImageBuffer(1).numFrames = imgBuffer;

% IMAGE BUFFER (PLANEWAVE)
Resource.ImageBuffer(2).datatype = 'double';
Resource.ImageBuffer(2).rowsPerFrame = ceil((SSA.rowsPerFrame/8)/PData(2).pdeltaZ);
Resource.ImageBuffer(2).colsPerFrame = PData(2).Size(2);
Resource.ImageBuffer(2).numFrames = imgBuffer;

% GUIDANCE IMAGING DISPLAY WINDOW (PHASED ARRAY AND PLANEWAVE)
Resource.DisplayWindow(1).Title = 'Guidance Imaging Display';
Resource.DisplayWindow(1).pdelta = 0.3;
Resource.DisplayWindow(1).Position = [200,200, ...
    ceil(PData(1).Size(2)*(PData(1).pdelta/Resource.DisplayWindow(1).pdelta)), ...
    ceil(PData(1).Size(1)*PData(1).pdelta/Resource.DisplayWindow(1).pdelta)];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),PData(1).Origin(3)];
Resource.DisplayWindow(1).Colormap = gray(256);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data Acquisition Buffers %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% RESOURCE BUFFER 2 (ALL SSA FRAMES)
Resource.RcvBuffer(2).datatype = 'int16';
Resource.RcvBuffer(2).rowsPerFrame = SSA.rowsPerFrame; 
Resource.RcvBuffer(2).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(2).numFrames = SSA.nFrames;

% NOTE: rows per frame indicates the number of samples of RF data collected 
% per channel. Default sample rate for A/D is 4x the center frequency of 
% the probe

% RESOURCE BUFFER 3 (ALL BMODE FRAMES)
Resource.RcvBuffer(3).datatype = 'int16';
Resource.RcvBuffer(3).rowsPerFrame = PHASED_B.rowsPerFrame; 
Resource.RcvBuffer(3).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(3).numFrames = bmodeSaveFrames;

% RESOURCE BUFFER 4 (ALL SA FRAMES)
Resource.RcvBuffer(4).datatype = 'int16';
Resource.RcvBuffer(4).rowsPerFrame = SA.rowsPerFrame; 
Resource.RcvBuffer(4).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(4).numFrames = saSaveFrames;

%%%%%%%%%%%%%%%%%%%%%%%
% Transmit Structures %
%%%%%%%%%%%%%%%%%%%%%%%

% TRANSMIT WAVEFORM DEFINITION (USED FOR ALL IMAGING MODES)
TW.type = 'parametric';
TW.Parameters = [Trans.frequency,.67,2,1]; % *** different from NB code ***

% NOTE:
% A,B,C,D in terms of master clock (180 MHz T = 5.556 ns)
% A = clocks in half cycle based on f0
% B = clocks during which transmit drivers are active
% C = number of half cycles in transmit waveform (2 = single cycle)
% D = polarity of first half cycle (1 = positive)
% **TW(1).Parameters = [36,17,2,1];   % A, B, C, D [old definition of TW]

% TRANSMIT STRUCTURE DEFINITION
numTxDef = PHASED_B.nRay+1+SSA.nAngles+SA.nRay;
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0,0,0], ...
                   'focus', 0, ... % in wavelengths; 
                   'Steer', [0,0], ...
                   'Apod', ones(1,numEl), ...  
                   'Delay', zeros(1,numEl)), 1, numTxDef); 

tx_i = 0;

% PHASED ARRAY IMAGING TX 
% - individual tx structs defined for 128 rays
tx_foc = round(PHASED_B.focusMM/1000/(c/(Trans.frequency*1e6)));    % [wavelengths]
bmodeAngles = SFormat(1).theta:SFormat(1).rayDelta:(SFormat(1).theta...
    + (PHASED_B.nRay-1)*SFormat(1).rayDelta);
origin = SFormat(1).radius*tan(bmodeAngles);

paTxStart = tx_i;     % start indices to be used in event structs
for n = 1:PHASED_B.nRay
    tx_i = tx_i+1;
    TX(tx_i).Origin = [origin(n),0.0,0.0];
    TX(tx_i).focus = tx_foc;
    TX(tx_i).Steer = [bmodeAngles(n),0.0];
    TX(tx_i).Delay = computeTXDelays(TX(tx_i));
end

% PLANEWAVE TX (PW SSA & PW GUIDANCE) 
% - single tx stuct defined
pwTxStart = tx_i;
tx_i = tx_i+1;
TX(tx_i).Delay = computeTXDelays(TX(tx_i));

% STEERED PW SSA TX 
% - individual tx structs defined for each steering angle
spwTxStart = tx_i;
for n = 1:SSA.nAngles
    tx_i = tx_i+1;
    TX(tx_i).Steer = [degtorad(SSA.pwAngles(n)),0.0]; % [edit] currently in radians
    TX(tx_i).Delay = computeTXDelays(TX(tx_i));
end

% FULL SA TX 
saTxStart = tx_i;
SA.txFocus = SA.txFnum*SA.numEl*Trans.spacing;
for n = 1:SA.nRay
    tx_i = tx_i+1;
    TX(tx_i).Origin = [SFormat(1).FirstRayLoc(1)+(n-1+floor(SA.numEl/2))*Trans.spacing, 0.0, 0.0];
    TX(tx_i).focus = SA.txFocus;
    TX(tx_i).Apod(:) = 0;
    TX(tx_i).Apod(n:n+SA.numEl-1) = 1.0;
    TX(tx_i).Delay = computeTXDelays(TX(tx_i));
end

% TRANSMIT POWER CONTROLLER 
TPC(1).name = 'B-mode';
TPC(1).maxHighVoltage = maxVoltage;
TPC(2).name = 'Swept Synthetic Aperture';
TPC(2).maxHighVoltage = maxVoltage;

% TIME GAIN CONTROL
TGC.CntrlPts = [450,550,650,710,770,830,890,950];
TGC.rangeMax = SFormat(1).endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

%%%%%%%%%%%%%%%%%%%%%%
% Receive Structures %
%%%%%%%%%%%%%%%%%%%%%%

numRcvDef = ...
    PHASED_B.nRay*Resource.RcvBuffer(1).numFrames ...   % frames for phased array imaging 
    + Resource.RcvBuffer(1).numFrames ...               % frames for pw imaging
    + Resource.RcvBuffer(2).numFrames ...               % frames for SSA acquisition
    + PHASED_B.nRay*Resource.RcvBuffer(3).numFrames...  % frames for phased array acquisition
    + SA.nRay*Resource.RcvBuffer(4).numFrames;          % frames for SA acquisition

% RECEIVE STRUCTURE DEFINITION
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
                        'callMediaFunc', 1),1,numRcvDef); 
rcv_i = 0;

% PHASED ARRAY GUIDANCE RCV 
maxAcqLengthBmode = sqrt(aperture^2 + SFormat(1).endDepth^2 ...
    - 2*aperture*SFormat(1).endDepth*cos(SFormat(1).theta-pi/2)) ...
    - SFormat(1).startDepth;

maxAcqLengthPW = sqrt(SFormat(2).endDepth^2+aperture^2)-SFormat(2).startDepth;
wlsPer128 = 128/(2*4);                  % wavelengths in 128 samples for 4 samplesPerWave

paGuideRcvStart = rcv_i;
for i = 1:Resource.RcvBuffer(1).numFrames
    k = PHASED_B.nRay*(i-1)+paGuideRcvStart;
    for j = 1:PHASED_B.nRay
        rcv_i = rcv_i+1;
        Receive(k+j).startDepth = SFormat(1).startDepth;
        Receive(k+j).endDepth = SFormat(1).startDepth + wlsPer128*ceil(maxAcqLengthBmode/wlsPer128);
        Receive(k+j).framenum = i;
        Receive(k+j).acqNum = j; 
        Receive(k+j).bufnum = 1;        % data stored resource buffer 1 (guidance imaging)
    end
end

% PLANEWAVE GUIDANCE RCV
pwGuideRcvStart = rcv_i;
for i = 1:Resource.RcvBuffer(1).numFrames
    k = i+pwGuideRcvStart;
    rcv_i = rcv_i+1;
    Receive(k).startDepth = SFormat(2).startDepth;
    Receive(k).endDepth = SFormat(2).startDepth + wlsPer128*ceil(maxAcqLengthPW/wlsPer128);
    Receive(k).framenum = i;
    Receive(k).bufnum = 1; 
end

% SSA ACQUISITION RCV 
% - single acquisition for multiple frames to be post-processed
ssaRcvStart = rcv_i;
for i = 1:Resource.RcvBuffer(2).numFrames
    k = i+ssaRcvStart;
    rcv_i = rcv_i+1;
    Receive(k).startDepth = SFormat(2).startDepth;
    Receive(k).endDepth = SFormat(2).startDepth + wlsPer128*ceil(maxAcqLengthPW/wlsPer128);
    Receive(k).framenum = i;      
    Receive(k).bufnum = 2; 
end

% PHASED ARRAY ACQUISITION RCV
paAcqRcvStart = rcv_i;
for i = 1:Resource.RcvBuffer(3).numFrames
    k = PHASED_B.nRay*(i-1)+paAcqRcvStart;
    for j = 1:PHASED_B.nRay
        rcv_i = rcv_i+1;
        Receive(k+j).startDepth = SFormat(1).startDepth;
        Receive(k+j).endDepth = SFormat(1).startDepth + wlsPer128*ceil(maxAcqLengthBmode/wlsPer128);
        Receive(k+j).framenum = i;
        Receive(k+j).acqNum = j; 
        Receive(k+j).bufnum = 3; % added to specify save b-mode buffer
    end
end

% SA ACQUISITION RCV
saRcvStart = rcv_i;
for i = 1:Resource.RcvBuffer(4).numFrames
    k = SA.nRay*(i-1)+saRcvStart;
    for j = 1:SA.nRay
        rcv_i = rcv_i+1;
        Receive(k+j).startDepth = SFormat(1).startDepth;
        Receive(k+j).endDepth = SFormat(1).startDepth + wlsPer128*ceil(maxAcqLengthBmode/wlsPer128);
        Receive(k+j).framenum = i;
        Receive(k+j).acqNum = j;
        Receive(k+j).bufnum = 4;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Image Reconstruction Structures %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PHASED ARRAY RECONSTRUCTION PARAMS
Recon = repmat(struct('senscutoff', 0.6, ...
               'pdatanum', [], ...
               'newFrameTimeout',1000,...
               'IntBufDest', [0,0], ...
               'ImgBufDest', [], ...
               'RINums', []),1,2);

Recon(1).pdatanum = 1;
Recon(1).rcvBufFrame = -1;
Recon(1).ImgBufDest = [1,-1];
Recon(1).RINums = 1:PHASED_B.nRay;
           
% PLANEWAVE RECONSTRUCTION PARAMS
Recon(2).pdatanum = 2;
Recon(2).rcvBufFrame = -1; 
Recon(2).ImgBufDest = [2,-1];
Recon(2).RINums = PHASED_B.nRay+1;

% RECON INFO STRUCTURE DEFINITION
ReconInfo = repmat(struct('mode', 0, ...    % replace data
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'regionnum', 0), 1, PHASED_B.nRay+1);

% PHASED ARRAY RECONIFO
for i = 1:PHASED_B.nRay
    ReconInfo(i).txnum = i;
    ReconInfo(i).rcvnum = i;
    ReconInfo(i).regionnum = i;
end

% PLANEWAVE RECONINFO
ReconInfo(PHASED_B.nRay+1).txnum = PHASED_B.nRay+1;
ReconInfo(PHASED_B.nRay+1).rxnum = 1;
ReconInfo(PHASED_B.nRay+1).regionnum = 0;

%%%%%%%%%%%%%%%%%%%%%%%
% PROCESS DEFINITIONS %
%%%%%%%%%%%%%%%%%%%%%%%
p_i = 1;

% PHASED ARRAY IMAGE DISPLAY 
pIdx.displayPhased = p_i;
pers = 10;                                      % persistence level 
Process(p_i).classname = 'Image';
Process(p_i).method = 'imageDisplay';
Process(p_i).Parameters = {'imgbufnum',1,...      % number of buffer to process.
                         'framenum',-1,...      % (-1 => lastFrame)
                         'pdatanum',1,...       % number of PData structure to use
                         'norm',1,...           % normalization method(1 means fixed)
                         'pgain',0.5,...        % pgain is image processing gain
                         'persistMethod','none',...
                         'persistLevel',pers,...
                         'interp',1,...         % method of interpolation (1=4pt interp)
                         'compression',0.5,...  % X^0.5 normalized to output word size
                         'mappingMode','full',...
                         'display',1,...        % display image after processing
                         'displayWindow',1};
p_i = p_i+1;

% PLANEWAVE IMAGE DISPLAY 
pIdx.displayPW = p_i;
Process(p_i).classname = 'Image';
Process(p_i).method = 'imageDisplay';
Process(p_i).Parameters = {'imgbufnum',2,...   % number of buffer to process.
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
p_i = p_i+1;
 
% SAVE B-MODE FRAMES (EXTERNAL)
pIdx.saveBmode = p_i;
Process(p_i).classname = 'External';
Process(p_i).method = 'saveBmode';
Process(p_i).Parameters = {'srcbuffer','receive',...
    'srcbufnum',3,...
    'srcframenum',0,... 
    'dstbuffer','none'};
p_i = p_i+1;
% NOTE: RData sent to external function determined by srcframenum
% 0 [all frames], -1 [last frame]

% SAVE B-MODE FRAMES AND IMAGE (EXTERNAL)
pIdx.saveBmodeFull = p_i;
Process(p_i).classname = 'External';
Process(p_i).method = 'saveBmodeFull';
Process(p_i).Parameters = {'srcbuffer','receive',...
    'srcbufnum',3,...
    'srcframenum',0,... 
    'dstbuffer','none'};
p_i = p_i+1;

% SAVE SSA FRAMES (EXTERNAL)
pIdx.saveSSA = p_i;
Process(p_i).classname = 'External';
Process(p_i).method = 'saveSSA';
Process(p_i).Parameters = {'srcbuffer','receive',...
    'srcbufnum',2,...
    'srcframenum',0,... 
    'dstbuffer','none'};
p_i = p_i+1;

% DISPLAY SSA FRAME (EXTERNAL) - for debugging purposed only
pIdx.displaySSA = p_i;
Process(p_i).classname = 'External';
Process(p_i).method = 'displaySSA';
Process(p_i).Parameters = {'srcbuffer','receive',...
    'srcbufnum',2,...
    'srcframenum',-1,... 
    'dstbuffer','none'};
p_i = p_i+1;

% SWEEP TURNTABLE ARM (EXTERNAL)
pIdx.sweepTable = p_i;
Process(p_i).classname = 'External';
Process(p_i).method = 'sweepTable';
Process(p_i).Parameters = {'srcbuffer','none',...
    'srcframenum','none',... 
    'dstbuffer','none'};
p_i = p_i+1;

% SAVE SSA FRAMES (EXTERNAL)
pIdx.saveSA = p_i;
Process(p_i).classname = 'External';
Process(p_i).method = 'saveSA';
Process(p_i).Parameters = {'srcbuffer','receive',...
    'srcbufnum',4,...
    'srcframenum',0,... 
    'dstbuffer','none'};
p_i = p_i+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SEQ CONTROL DEFINITIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PHASED ARRAY TIMING
nsc = 1;
scIdx.delayLinePA = nsc;                     
SeqControl(nsc).command = 'timeToNextAcq';
SeqControl(nsc).argument = PHASED_B.tPulse;  
nsc = nsc+1;

scIdx.delayFramePA = nsc; 
SeqControl(nsc).command = 'timeToNextAcq';
SeqControl(nsc).argument = round((1/PHASED_B.PRF)*1e6 - 128*PHASED_B.tPulse);
nsc = nsc+1;

% SWITCH TPC CONTROLLER 2 (SSA)
scIdx.profileSSA = nsc; 
SeqControl(nsc).command = 'setTPCProfile';
SeqControl(nsc).condition = 'next';
SeqControl(nsc).argument = 2;
nsc = nsc+1;

% SWITCH TPC CONTROLLER (B-mode)
scIdx.profileBmode = nsc; 
SeqControl(nsc).command = 'setTPCProfile';
SeqControl(nsc).condition = 'next';
SeqControl(nsc).argument = 1;
nsc = nsc+1;

% SSA TIMING
scIdx.delayFrameSSA = nsc; 
SeqControl(nsc).command = 'timeToNextAcq';
SeqControl(nsc).argument = SSA.PRT;
nsc = nsc+1;

% WAIT TIME FOR DELAYING SEQUENCE
scIdx.waitNOOP = nsc; 
SeqControl(nsc).command = 'noop';
SeqControl(nsc).argument = 524287;
nsc = nsc+1;

% OUTPUT TRIGGER 
scIdx.triggerOut = nsc; 
SeqControl(nsc).command = 'triggerOut';
nsc = nsc+1;

% PLANEWAVE TIMING
scIdx.delayFramePW = nsc; 
SeqControl(nsc).command = 'timeToNextAcq'; 
SeqControl(nsc).argument = PW_B.PRT;
nsc = nsc+1;

% PLANEWAVE TIMING
scIdx.delayFrameSA = nsc; 
SeqControl(nsc).command = 'timeToNextAcq'; 
SeqControl(nsc).argument = round(1/SA.PRF*1e6 - 64*PHASED_B.tPulse);
nsc = nsc+1;

% ****************** GUIDANCE PHASED ARRAY BMODE EVENTS *******************
n = 1;

bmode_image = n;

Event(n).info = 'Switch to b-mode acquire';
Event(n).tx = 0;         
Event(n).rcv = 0;       
Event(n).recon = 0;      
Event(n).process = 0;
Event(n).seqControl = scIdx.profileBmode; 
n = n+1;

bmode_loop = n;
for i = 1:Resource.RcvBuffer(1).numFrames
    for j = 1:PHASED_B.nRay
        Event(n).info = 'Acquire b-mode ray line';
        Event(n).tx = paTxStart + j; 
        Event(n).rcv = paGuideRcvStart+PHASED_B.nRay*(i-1)+j;   
        Event(n).recon = 0;      
        Event(n).process = 0;    
        Event(n).seqControl = scIdx.delayLinePA; 
        n = n+1;
    end
    
    % Replace last event SeqControl to invoke delay between frames.
    Event(n-1).seqControl = [scIdx.delayFramePA,nsc]; 
    SeqControl(nsc).command = 'transferToHost'; 
    nsc = nsc+1;

    Event(n).info = 'Recon and process'; 
    Event(n).tx = 0;         
    Event(n).rcv = 0;        
    Event(n).recon = 1;      
    Event(n).process = pIdx.displayPhased; 
    Event(n).seqControl = 0;
    if floor(i/2) == i/2 
        Event(n).seqControl = nsc;
        SeqControl(nsc).command = 'returnToMatlab';
        nsc = nsc+1;
    end
    n = n+1;
end

Event(n).info = 'Jump back to guidance bmode start';
Event(n).tx = 0;        
Event(n).rcv = 0;       
Event(n).recon = 0;     
Event(n).process = 0; 
Event(n).seqControl = nsc;  
    SeqControl(nsc).command = 'jump';
    SeqControl(nsc).argument = bmode_loop;
    nsc = nsc+1;    
n = n+1;

% ****************** GUIDANCE PLANEWAVE BMODE EVENTS **********************
pw_image = n;

Event(n).info = 'Switch to SA acquire';
Event(n).tx = 0;         
Event(n).rcv = 0;       
Event(n).recon = 0;      
Event(n).process = 0;
Event(n).seqControl = scIdx.profileSSA;
n = n+1;

pw_loop = n;
for i = 1:Resource.RcvBuffer(1).numFrames
    Event(n).info = 'Acquire PW image';
    Event(n).tx = pwTxStart+1; 
    Event(n).rcv = pwGuideRcvStart+i;   
    Event(n).recon = 0;      
    Event(n).process = 0;    
    Event(n).seqControl = [scIdx.delayFramePW,nsc];
        SeqControl(nsc).command = 'transferToHost';
        nsc = nsc+1;
    n = n+1;

    Event(n).info = 'Recon and process'; 
    Event(n).tx = 0;         
    Event(n).rcv = 0;        
    Event(n).recon = 2;      
    Event(n).process = pIdx.displayPW;    
    Event(n).seqControl = 0;
    if floor(i/5) == i/5
        Event(n).seqControl = nsc;
        SeqControl(nsc).command = 'returnToMatlab';
        nsc = nsc+1;
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

%******************** PHASED ARRAY ACQUISITION EVENTS *********************
bmode_acq = n;

Event(n).info = 'Switch to b-mode acquire';
Event(n).tx = 0;         
Event(n).rcv = 0;       
Event(n).recon = 0;      
Event(n).process = 0;
Event(n).seqControl = scIdx.profileBmode;
n = n+1;

for i = 1:Resource.RcvBuffer(3).numFrames
    for j = 1:PHASED_B.nRay                 
        Event(n).info = 'Acquire ray line';
        Event(n).tx = paTxStart+j;         
        Event(n).rcv = paAcqRcvStart+PHASED_B.nRay*(i-1)+j;   
        Event(n).recon = 0;      
        Event(n).process = 0;    
        Event(n).seqControl = scIdx.delayLinePA;
        n = n+1;
    end

Event(n-1).seqControl = [scIdx.delayFramePA,nsc];
    SeqControl(nsc).command = 'transferToHost';
    SeqControl(nsc).condition = 'waitForProcessing';
    SeqControl(nsc).argument = nsc;
    nsc = nsc+1;

Event(n).info = 'Wait for transfer complete';
Event(n).tx = 0;         
Event(n).rcv = 0;        
Event(n).recon = 0;      
Event(n).process = 0;    
Event(n).seqControl = [nsc,nsc+1]; 
    SeqControl(nsc).command = 'waitForTransferComplete';
    SeqControl(nsc).argument = nsc-1;
    SeqControl(nsc+1).command = 'markTransferProcessed';
    SeqControl(nsc+1).argument = nsc-1;
    nsc = nsc+2;
n = n+1;

end

Event(n).info = 'Save full B-mode data (RF and image)';
Event(n).tx = 0; 
Event(n).rcv = 0; 
Event(n).recon = 0; 
Event(n).process = pIdx.saveBmodeFull; 
Event(n).seqControl = nsc;
    SeqControl(nsc).command = 'returnToMatlab';
    nsc = nsc+1;
n = n+1;

Event(n).info = 'Jump back to guidance b-mode'; 
Event(n).tx = 0;       
Event(n).rcv = 0;
Event(n).recon = 0;     
Event(n).process = 0; 
Event(n).seqControl = nsc;
    SeqControl(nsc).command = 'jump';
    SeqControl(nsc).argument = bmode_image;
    nsc = nsc+1;
n = n+1;

%************************ SSA ACQUISITION EVENTS **************************
% IMAGING CASE 1: PLANEWAVE SSA
% IMAGING CASE 2: STEERED PLANEWAVE SSA

for ic = 1:2
    
%   CASE SPECIFIC PARAMETERS [EDIT]
    if ic == 1,             
        SSA_acq = n;
    elseif ic == 2,         
        steerSSA_acq = n;
    end

    Event(n).info = 'Switch to SA acquire';
    Event(n).tx = 0;         
    Event(n).rcv = 0;       
    Event(n).recon = 0;      
    Event(n).process = 0;
    Event(n).seqControl = scIdx.profileSSA;
    n = n+1;

    switch SETUP.scanType
        case 'manual'

        case 'turntable'
            Event(n).info = 'Trigger probe sweep';
            Event(n).tx = 0;         
            Event(n).rcv = 0;       
            Event(n).recon = 0;      
            Event(n).process = pIdx.sweepTable;
            Event(n).seqControl = nsc;
                SeqControl(nsc).command = 'returnToMatlab';
                nsc = nsc+1;
            n = n+1;
    end

    [Event,n] = genWaitEvents(Event,n,scIdx.waitNOOP,3);

%   CASE SPECIFIC PARAMETERS [EDIT]
    if ic == 1,         
        SSA_loop = n;
    elseif ic ==2,
        steerSSA_loop = n;
        j=1;
        steerAngles = zeros(1,Resource.RcvBuffer(2).numFrames);
    end
    
    for i = 1:Resource.RcvBuffer(2).numFrames
        if i > 1
            Event(n).info = 'Send trigger out.';
            Event(n).tx = 0;
            Event(n).rcv = 0; 
            Event(n).recon = 0; 
            Event(n).process = 0; 
            Event(n).seqControl = scIdx.triggerOut; 
            n = n+1;
            
%           CASE SPECIFIC PARAMETERS [EDIT]
            if ic == 1
                Event(n).info = 'Acquire planewave SSA RF';
                Event(n).tx = pwTxStart+1;
                Event(n).rcv = i+ssaRcvStart; 
                Event(n).recon = 0; 
                Event(n).process = 0; 
                Event(n).seqControl = [scIdx.delayFrameSSA,nsc]; 
                    SeqControl(nsc).command = 'transferToHost';
                    nsc = nsc+1;
                n = n+1;
            elseif ic == 2
                Event(n).info = 'Acquire steered planewave SSA RF';
                Event(n).tx = spwTxStart+j;
                Event(n).rcv = i+ssaRcvStart; 
                Event(n).recon = 0; 
                Event(n).process = 0; 
                Event(n).seqControl = [scIdx.delayFrameSSA,nsc]; 
                    SeqControl(nsc).command = 'transferToHost';
                    nsc = nsc+1;
                n = n+1;

                steerAngles(i-1) = SSA.pwAngles(j); 
                j = j+1;
                if j > SSA.nAngles; j = 1; end
            end
        end
    end

    Event(n).info = 'Wait and sync.';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = [nsc nsc+1 nsc+2];
        SeqControl(nsc).command = 'waitForTransferComplete';
        SeqControl(nsc).argument = nsc-1;
        SeqControl(nsc+1).command = 'sync';
        SeqControl(nsc+2).command = 'returnToMatlab';
        nsc = nsc+3;
    n = n+1;

    Event(n).info = 'Save SSA frames'; 
    Event(n).tx = 0;         
    Event(n).rcv = 0;        
    Event(n).recon = 0;      
    Event(n).process = pIdx.saveSSA;
    Event(n).seqControl = nsc;
        SeqControl(nsc).command = 'returnToMatlab';
        nsc = nsc+1;
    n = n+1;

    Event(n).info = 'Jump back to guidance b-mode';
    Event(n).tx = 0;        
    Event(n).rcv = 0;       
    Event(n).recon = 0;     
    Event(n).process = 0; 
    Event(n).seqControl = nsc;
        SeqControl(nsc).command = 'jump';
        SeqControl(nsc).argument = bmode_image;
        nsc = nsc+1;
    n = n+1;
end

%************************* SA ACQUISITION EVENTS **************************
SA_acq = n;

% Use same TCP profile as SSA

Event(n).info = 'Switch to SA acquire';
Event(n).tx = 0;         
Event(n).rcv = 0;       
Event(n).recon = 0;      
Event(n).process = 0;
Event(n).seqControl = scIdx.profileSSA;
n = n+1;
    
for i = 1:Resource.RcvBuffer(4).numFrames
    for j = 1:SA.nRay
        Event(n).info = 'Acquire ray line';
        Event(n).tx = saTxStart+j; 
        Event(n).rcv = saRcvStart+SA.nRay*(i-1)+j;
        Event(n).recon = 0; 
        Event(n).process = 0; 
        Event(n).seqControl = scIdx.delayLinePA; 
        n = n+1;
    end
    Event(n-1).seqControl = [scIdx.delayFrameSA,nsc];
        SeqControl(nsc).command = 'transferToHost'; 
        nsc = nsc+1;
end
Event(n-1).seqControl = [scIdx.delayFrameSA,nsc-1];
    SeqControl(nsc-1).condition = 'waitForProcessing'; 
    SeqControl(nsc-1).argument = nsc-1;

Event(n).info = 'Wait for transfer complete';
Event(n).tx = 0; 
Event(n).rcv = 0; 
Event(n).recon = 0; 
Event(n).process = 0; 
Event(n).seqControl = [nsc,nsc+1]; 
    SeqControl(nsc).command = 'waitForTransferComplete';
    SeqControl(nsc).argument = nsc-1;
    SeqControl(nsc+1).command = 'markTransferProcessed';
    SeqControl(nsc+1).argument = nsc-1;
    nsc = nsc+2;
n = n+1;

Event(n).info = 'Save SA frame';
Event(n).tx = 0; 
Event(n).rcv = 0; 
Event(n).recon = 0; 
Event(n).process = pIdx.saveSA; 
Event(n).seqControl = nsc;
    SeqControl(nsc).command = 'returnToMatlab';
    nsc = nsc+1;
n = n+1;

Event(n).info = 'Jump back to guidance b-mode';
Event(n).tx = 0;        
Event(n).rcv = 0;       
Event(n).recon = 0;     
Event(n).process = 0; 
Event(n).seqControl = nsc;
    SeqControl(nsc).command = 'jump';
    SeqControl(nsc).argument = bmode_image;
    nsc = nsc+1;
n = n+1;

%************************ GUI elements ************************************
ui = 1;

UI(ui).Control = {'UserB3','Style','VsPushButton','Tag','pwImage','Label','Toggle B-mode'};
UI(ui).Callback = text2cell('%CB_imagingToggle');
ui = ui+1;

UI(ui).Control = {'UserB2','Style','VsPushButton','Tag','saveBmodeFull','Label','B-mode'};
UI(ui).Callback = text2cell('%CB_saveBmodeFull');
ui = ui+1;

UI(ui).Control = {'UserB1','Style','VsPushButton','Tag','saveSA','Label','SA'}; % change UI location
UI(ui).Callback = text2cell('%CB_saveSA');
ui = ui+1;

switch SETUP.scanType
    case 'manual'
        UI(ui).Control = {'UserC2','Style','VsPushButton','Tag','pwSSA','Label','pw SSA'};
        UI(ui).Callback = text2cell('%CB_pwSSA');
        ui = ui+1;

        UI(ui).Control = {'UserC1','Style','VsPushButton','Tag','steerSSA','Label','steer SSA'};
        UI(ui).Callback = text2cell('%CB_steerSSA');
        ui = ui+1;
        
    case 'turntable'
        UI(ui).Control = {'UserC2','Style','VsPushButton','Tag','pwSSA','Label','pw SSA'};
        UI(ui).Callback = text2cell('%CBtable_pwSSA');
        ui = ui+1;

        UI(ui).Control = {'UserC1','Style','VsPushButton','Tag','steerSSA','Label','steer SSA'};
        UI(ui).Callback = text2cell('%CBtable_steerSSA');
        ui = ui+1;

        UI(ui).Control = {'UserB6','Style','VsPushButton','Tag','testSweep','Label','Test Sweep'};
        UI(ui).Callback = text2cell('%CB_testSweep');
        ui = ui+1;

        UI(ui).Control = {'UserB4','Style','VsPushButton','Tag','resetRS232','Label','Reset RS232'};
        UI(ui).Callback = text2cell('%CB_resetRS232');
        ui = ui+1;
        
        UI(ui).Control = {'UserB5','Style','VsPushButton','Tag','mvCenter','Label','Center Probe'};
        UI(ui).Callback = text2cell('%CB_mvCenter');
        ui = ui+1;
end

%************************ output VSX **************************************
scriptName = sprintf('P4-2_SSA_%s',SETUP.scanType); 
% 'P4-2_planewaveSSA';
disp(['filename =''' scriptName ''';VSX'])
save(scriptName);
return 

%CB_pwSSA
init_pwSSA(hObject,eventdata)
return
%CB_pwSSA

%CB_steerSSA
init_steerSSA(hObject,eventdata)
return
%CB_steerSSA

%CB_saveBmodeFull
init_saveBmodeFull(hObject,eventdata)
return
%CB_saveBmodeFull

%CB_saveSA
init_saveSA(hObject,eventdata)
return
%CB_saveSA

%CB_imagingToggle
imagingToggle(hObject,eventdata)
return
%CB_imagingToggle

%CBtable_steerSSA
table_initRS232
init_steerSSA(hObject,eventdata)
return
%CBtable_steerSSA

%CBtable_pwSSA
table_initRS232
init_pwSSA(hObject,eventdata)
return
%CBtable_pwSSA

%CB_testSweep
table_initRS232
table_testSweep(hObject,eventdata)
return
%CB_testSweep

%CB_resetRS232
table_resetRS232(hObject,eventdata)
return
%CB_resetRS232

%CB_mvCenter
table_initRS232
table_mvCenter(hObject,eventdata)
return
%CB_mvCenter