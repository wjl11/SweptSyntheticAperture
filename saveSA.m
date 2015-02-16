function saveSA(RData)
label = evalin('base','saveLabel');
c = evalin('base','c');
nr = evalin('base','SA.nRay');
Trans = evalin('base', 'Trans');
Receive = evalin('base', 'Receive');
Resource = evalin('base', 'Resource');
SFormat = evalin('base', 'SFormat');
TW = evalin('base', 'TW');
TX = evalin('base', 'TX');
txFocus = evalin('base', 'SA.focusMM');
rcv_i = evalin('base','saRcvStart');

persistent nframeSA
if isempty(nframeSA); nframeSA = 1; end
dir = './data/';
if exist(dir,'file')~=7; mkdir(dir); end
name = ['sa_' label '_' datestr(now,'yyyymmdd_HHMMSS') '_' num2str(nframeSA)];
path = [dir name];

numRcvSamples = Receive(rcv_i+1).endSample-Receive(rcv_i+1).startSample+1;
tmp = RData(1:(numRcvSamples*nr),[1:32 97:128],1);
tmp_rf = reshape(tmp,[numRcvSamples,nr,64]);
tmp_rf = permute(tmp_rf,[1 3 2]);
rf = tmp_rf;
clear tmp tmp_rf

rfdata.c = c;
rfdata.numElementsPerXmt= 64;
rfdata.numXmtRxEvents = nr;
rfdata.numFrames = Resource.RcvBuffer(4).numFrames;
rfdata.elementSpacingMM = Trans.spacingMm;
rfdata.XMTspacingMM = rfdata.elementSpacingMM;
rfdata.samplingRateMHz = Trans.frequency*Receive(rcv_i+1).samplesPerWave;
rfdata.frequencyMHz = Trans.frequency;
rfdata.focusMM = txFocus*c/(Trans.frequency*1e6);
rfdata.timeZero = -(SFormat(1).startDepth+...
                TW(1).peak+...
                min(TX(end).Delay))*Receive(rcv_i+1).samplesPerWave;

if strcmpi(label,'db') || strcmpi(label,'')
    disp('[DEBUG MODE] No file saved.')
else
    disp(['Saving SA frame to ' path]);
    save([path '.mat'],'rf','rfdata');
    disp(['SA data saved to ' path]);
end

nframeSA = nframeSA+1;
