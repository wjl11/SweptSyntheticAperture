function P4_2_save_bmode(RData)
label = evalin('base','save_label');
SAVE_STATE = evalin('base','SAVE_STATE');
persistent nframeBmode

if isempty(nframeBmode)
    nframeBmode = 1;
end

filename = ['Bmode_frame_' label '_' num2str(nframeBmode) '.mat'];
disp(['Saving B-mode frame to ' filename]);
try
    save(filename,'RData');
    SAVE_STATE(1) = 1;
catch
    error('Failed to save b-mode RF.');
end

assignin('base','SAVE_STATE',SAVE_STATE);
nframeBmode = nframeBmode+1;
