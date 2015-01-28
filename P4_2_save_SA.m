function P4_2_save_SA(RData)
SAVE_STATE = evalin('base','SAVE_STATE');
SAdata = evalin('base','SAdata');
label = evalin('base','save_label');
persistent nframeSA

if isempty(nframeSA)
    nframeSA = 1;
end

filename = ['SA_frame_' label '_' num2str(nframeSA) '.mat'];
disp(['Saving SA frame to ' filename]);
try
    save(filename,'SAdata');
    SAVE_STATE(2) = 1;
catch
    error('Failed to save SA RF.');
end

clear SAdata
SAdata.rf = [];
SAdata.th = [];
SAdata.t = [];

assignin('base','SAdata',SAdata);
assignin('base','SAVE_STATE',SAVE_STATE);

nframeSA = nframeSA+1;





