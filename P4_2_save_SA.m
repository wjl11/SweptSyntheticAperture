function P4_2_save_SA(RData)
label = evalin('base','saveLabel');
c = evalin('base','c');
Trans = evalin('base', 'Trans');
Receive = evalin('base', 'Receive');
ssaPRT = evalin('base','ssaPRT');
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
% return to imaging
% Control = evalin('base','Control');
% pw_image = evalin('base','pw_image');
% bmode_image = evalin('base','bmode_image');
% 
% if strcmpi(IM_STATE,'pw')
%     Control(1).Command = 'set&Run';
%     Control(1).Parameters = {'Parameters',1,'startEvent',pw_image};
%     evalin('base','Resource.Parameters.startEvent = pw_image;');
%     assignin('base','Control',Control);
% 
%     disp(['Jump to event: ' num2str(pw_image)])
%     disp('Switch to PW imaging mode.')
%     
% elseif strcmpi(IM_STATE,'pa')
%     Control(1).Command = 'set&Run';
%     Control(1).Parameters = {'Parameters',1,'startEvent',bmode_image};
%     evalin('base','Resource.Parameters.startEvent = bmode_image;');
%     assignin('base','Control',Control);
% 
%     disp(['Jump to event: ' num2str(bmode_image)])
%     disp('Switch to PA B-mode imaging mode.')
% else
%     error('Imaging state invalid.')
% end



