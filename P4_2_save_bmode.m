function P4_2_save_bmode(RData)
label = evalin('base','save_label');
persistent nframeBmode

if isempty(nframeBmode)
    nframeBmode = 1;
end

filename = ['Bmode_frame_' label '_' num2str(nframeBmode) '.mat'];
disp(['Saving B-mode frame to ' filename]);
rf = RData(:,[1:32 97:128],:);
try
    save(filename,'rf');
catch
    error('Failed to save b-mode RF.');
end

nframeBmode = nframeBmode+1;
