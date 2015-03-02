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
txFocus = evalin('base', 'SA.txFocus');
rcv_i = evalin('base','saRcvStart')+1;
tx_i = evalin('base','saTxStart')+1;

persistent nframeSA
if isempty(nframeSA); nframeSA = 0; end
dir = './data/';
if exist(dir,'file')~=7; mkdir(dir); end
name = ['sa_' label '_' datestr(now,'yyyymmdd_HHMMSS') '_' num2str(nframeSA)];
path = [dir name];

numRcvSamples = Receive(rcv_i).endSample-Receive(rcv_i).startSample+1;
tmp = RData(1:(numRcvSamples*nr),[1:32 97:128],1);
tmp_rf = reshape(tmp,[numRcvSamples,nr,64]);
tmp_rf = permute(tmp_rf,[1 3 2]);
rf = tmp_rf;
clear tmp tmp_rf

rfdata.c = c;
rfdata.numRcvChannels = 64;
rfdata.numElementsPerXmt= 64;
rfdata.numXmtRxEvents = nr;
rfdata.numFrames = Resource.RcvBuffer(4).numFrames;
rfdata.elementSpacingMM = Trans.spacingMm;
rfdata.XMTspacingMM = rfdata.elementSpacingMM;
rfdata.samplingRateMHz = Trans.frequency*Receive(rcv_i).samplesPerWave;
rfdata.frequencyMHz = Trans.frequency;
rfdata.timeZero = -(SFormat(1).startDepth+...
                Trans.lensCorrection*2+...
                TW(1).peak+...
                min(TX(tx_i+nr-1).Delay))*Receive(rcv_i).samplesPerWave;
            
rfdata.focus = txFocus*c/(Trans.frequency*1e6);

if strcmpi(label,'db') || strcmpi(label,'')
    disp('[DEBUG MODE] No file saved.')
else
    disp(['Saving SA frame to ' path]);
    save([path '.mat'],'rf','rfdata','-v7.3');
    disp(['SA data saved to ' path]);
end

nframeSA = nframeSA+1;
