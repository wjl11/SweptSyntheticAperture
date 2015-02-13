function showSSA(RData)
persistent fig
if isempty(fig)
    fig = figure;
end

set(0,'CurrentFigure',fig);
imagesc(RData); colormap gray;
disp('Process SSA frames')
