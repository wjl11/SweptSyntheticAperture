function saveSSA(RData)
label = evalin('base','saveLabel');
c = evalin('base','c');
Trans = evalin('base', 'Trans');
Receive = evalin('base', 'Receive');
ssaPRT = evalin('base','SSA.PRT');
steerAngles = evalin('base','steerAngles');
IM_STATE = evalin('base','IM_STATE');
SSA_TYPE = evalin('base','SSA_TYPE');

persistent nframeSA
if isempty(nframeSA); nframeSA = 1;end
dir = './data/';
if exist(dir,'file')~=7; mkdir(dir); end
name = ['ssa_' label '_' datestr(now,'yyyymmdd_HHMMSS') '_' num2str(nframeSA)];
path = [dir name];

rf = RData(:,[1:32 97:128],:);
rfdata.c = c;
rfdata.numRcvChannels = 64;
rfdata.elementSpacingMM = Trans.spacingMm;
rfdata.XMTspacingMM = rfdata.elementSpacingMM;
rfdata.samplingRateMHz = Trans.frequency*Receive(1).samplesPerWave;
rfdata.frequencyMHz = Trans.frequency;
rfdata.frameRatekHz = 1/(ssaPRT*1e-3);
rfdata.type = SSA_TYPE;

if strcmpi(SSA_TYPE,'steer_pw'); rfdata.steerAngles = steerAngles;
end

if strcmpi(label,'db') || strcmpi(label,'')
    disp('[DEBUG MODE] No file saved.')
else
    disp(['Saving SA data to ' path '.mat']);
    save([path '.mat'],'rf','rfdata');
    disp(['SA data saved to ' path '.mat']);
end
nframeSA = nframeSA+1;

IM_STATE = 'pa';
assignin('base','IM_STATE',IM_STATE);



