function frame=figure_frame(tl,handles,names,titleText,numColumns)
%FIGURE_FRAME  Single horizontal legend for a profile figure.
%   The panels themselves are drawn by profile_panel.m, so colours, line
%   styles, markers and ticks are already identical across groups.  This
%   routine only places the shared legend used by every figure of the paper:
%   one horizontal row, Times New Roman 9 pt, no box, in the north tile.
%
%   No title is drawn inside the figure: every figure is named by its LaTeX
%   caption, so an in-figure title would duplicate it.  The titleText input
%   is kept for call compatibility and ignored.
%
%   tl          tiledlayout handle of the figure
%   handles     line handles returned by profile_panel (one per curve)
%   names       legend labels, in the same order
%   titleText   unused, retained for call compatibility
%   numColumns  legend columns (four to six entries fit in one row)
if nargin<5 || isempty(numColumns), numColumns=numel(names); end
lg=legend(handles,names,'Orientation','horizontal','NumColumns',numColumns, ...
    'FontName','Times New Roman','FontSize',9,'Box','off');
lg.Layout.Tile='north';
frame=struct('Legend',lg,'Title','');
end
