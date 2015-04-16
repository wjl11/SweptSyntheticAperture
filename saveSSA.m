function saveSSA(RData)
label = evalin('base','saveLabel');
c = evalin('base','c');
Trans = evalin('base', 'Trans');
Receive = evalin('base', 'Receive');
Resource = evalin('base', 'Resource');
SFormat = evalin('base', 'SFormat');
TW = evalin('base', 'TW');
TX = evalin('base', 'TX');
rcv_i = evalin('base','ssaRcvStart')+1;

ssaPRT = evalin('base','SSA.PRT');
IM_STATE = evalin('base','IM_STATE');
SSA_TYPE = evalin('base','SSA_TYPE');

switch SSA_TYPE
    case 'pw'
        tx_i = evalin('base','pwTxStart')+1;
        txFocus = 0;
        rfdata.numElementsPerXmt = 2; % CHANGE IF PW USED IN FUTURE
    case 'steer_pw'
        tx_i = evalin('base','spwTxStart')+1;
        steerAngles = evalin('base','steerAngles');
        rfdata.steerAngles = steerAngles;
        txFocus = 0;
        rfdata.numElementsPerXmt = 64;
    case 'div'
        tx_i = evalin('base','sDivTxStart')+1;
        txFocus = evalin('base', 'divSSA.txFocus');
        rfdata.numElementsPerXmt = evalin('base', 'divSSA.numEl');
    case 'focus'
        tx_i = evalin('base','sFocTxStart')+1;
        txFocus = evalin('base', 'focSSA.txFocus');
        rfdata.numElementsPerXmt = evalin('base', 'focSSA.numEl');
end

persistent nframeSSA
if isempty(nframeSSA); nframeSSA = 0;end
dir = './data/';
if exist(dir,'file')~=7; mkdir(dir); end
name = ['ssa_' label '_' datestr(now,'yyyymmdd_HHMMSS') '_' num2str(nframeSSA)];
path = [dir name];

rf = RData(:,[1:32 97:128],2:end);
rfdata.c = c;
rfdata.numRcvChannels = 64;
% rfdata.numElementsPerXmt= 64;
rfdata.numFrames = Resource.RcvBuffer(2).numFrames;
rfdata.numXmtRxEvents = rfdata.numFrames;
rfdata.elementSpacingMM = Trans.spacingMm;
rfdata.XMTspacingMM = rfdata.elementSpacingMM;
rfdata.samplingRateMHz = Trans.frequency*Receive(rcv_i).samplesPerWave;
rfdata.frequencyMHz = Trans.frequency;
rfdata.frameRatekHz = 1/(ssaPRT*1e-3);
rfdata.timeZero = -(SFormat(2).startDepth+...
                Trans.lensCorrection*2+...
                TW(1).peak+...
                min(TX(tx_i).Delay))*Receive(rcv_i).samplesPerWave;
rfdata.focus = txFocus*c/(Trans.frequency*1e6);
rfdata.type = SSA_TYPE;


if strcmpi(label,'db') || strcmpi(label,'')
    disp('[DEBUG MODE] No file saved.')
else
    saveOpt = input('Save SSA data [y/n]? ','s')
    if strcmpi(saveOpt,'y')
        disp(['Saving SSA data to ' path '.mat']);
        save([path '.mat'],'rf','rfdata','-v7.3');
        disp(['SSA data saved to ' path '.mat']);
    else
        disp('SSA data not saved.');
    end
end

nframeSSA = nframeSSA+1;

IM_STATE = 'pa';
assignin('base','IM_STATE',IM_STATE);



