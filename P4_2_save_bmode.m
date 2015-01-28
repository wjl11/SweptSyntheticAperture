function P4_2_save_bmode(RData)
label = evalin('base','saveLabel');
c = evalin('base','c');
nr = evalin('base','nr');
angles = evalin('base','bmodeAngles');
origin = evalin('base','origin');
Trans = evalin('base', 'Trans');
Receive = evalin('base', 'Receive');
txFocus = evalin('base', 'bmodeFocusMM');

persistent nframeBmode
if isempty(nframeBmode); nframeBmode = 1; end
dir = './data/';
if exist(dir,'file')~=7; mkdir(dir); end
name = ['bmode_' label '_' datestr(now,'yyyymmdd_HHMMSS') '_' num2str(nframeBmode)];
path = [dir name];

numRcvSamples = Receive(1).endSample-Receive(1).startSample+1;
tmp = RData(1:(numRcvSamples*nr),[1:32 97:128],1);
rf = reshape(tmp,[numRcvSamples,nr,64]);
rf = permute(rf,[1 3 2]);

rfdata.c = c;
rfdata.numRcvChannels = 64;
rfdata.numXmtRxEvents = nr;
rfdata.elementSpacingMM = Trans.spacingMm;
rfdata.XMTspacingMM = rfdata.elementSpacingMM;
rfdata.samplingRateMHz = Trans.frequency*Receive(1).samplesPerWave;
rfdata.frequencyMHz = Trans.frequency;
rfdata.focusMM = txFocus;
rfdata.theta = angles;
% rfdata.origin = origin;

if strcmpi(label,'db') || strcmpi(label,'')
    disp('[DEBUG MODE] No file saved.')
else
    disp(['Saving B-mode frame to ' path]);
    save([path '.mat'],'rf','rfdata');
    figure(2)
    print('-dpng',[path '_verasonics.png'])
    disp(['B-mode data saved to ' path]);
end

nframeBmode = nframeBmode+1;
