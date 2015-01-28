function P4_2_save_SA(RData)
persistent nframe

if isempty(nframe)
    nframe = 1;
end

filename = ['SA_frame_' num2str(nframe) '.mat'];
disp(['Saving SA frame to ' filename]);
save(filename,'RData');

nframe = nframe+1;





