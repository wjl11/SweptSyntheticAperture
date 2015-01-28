function P4_2_image_SA(RData)
persistent fig
if isempty(fig)
    fig = figure;
end

set(0,'CurrentFigure',fig);
imagesc(RData); colormap gray;
disp('Process SA frames')
