function P4_2_save_SA(RData)
label = evalin('base','save_label');
persistent nframeSA

if isempty(nframeSA)
    nframeSA = 1;
end

filename = ['SA_frame_' label '_' num2str(nframeSA) '.mat'];
disp(['Saving SA frame to ' filename]);
rf = RData(:,[1:32 97:128],:);
try
    save(filename,'rf');
catch
    error('Failed to save SA RF.');
end

nframeSA = nframeSA+1;





