function saveSSA(RData)
label = evalin('base','saveLabel');
c = evalin('base','c');
Trans = evalin('base', 'Trans');
Receive = evalin('base', 'Receive');
Resource = evalin('base', 'Resource');
SFormat = evalin('base', 'SFormat');
TW = evalin('base', 'TW');
TX = evalin('base', 'TX');
txFocus = evalin('base', 'SSA.txFocus');
rcv_i = evalin('base','ssaRcvStart')+1;
tx_i = evalin('base','pwTxStart')+1;

ssaPRT = evalin('base','SSA.PRT');
steerAngles = evalin('base','steerAngles');
IM_STATE = evalin('base','IM_STATE');
SSA_TYPE = evalin('base','SSA_TYPE');


if strcmpi(SSA_TYPE,'steer_pw'); tx_i = evalin('base','spwTxStart')+1;
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
rfdata.numElementsPerXmt= 64;
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

if strcmpi(SSA_TYPE,'steer_pw'); rfdata.steerAngles = steerAngles;
end

if strcmpi(label,'db') || strcmpi(label,'')
    disp('[DEBUG MODE] No file saved.')
else
    disp(['Saving SA data to ' path '.mat']);
    save([path '.mat'],'rf','rfdata','-v7.3');
    disp(['SA data saved to ' path '.mat']);
end
nframeSSA = nframeSSA+1;

IM_STATE = 'pa';
assignin('base','IM_STATE',IM_STATE);



