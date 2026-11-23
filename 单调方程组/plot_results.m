function figs = plot_results(outdir)
%PLOT_RESULTS  Unified Dolan-More profiles for the eight monotone problems.
%
% The panels are drawn by profile_panel.m, i.e. with exactly the style block
% used by the Group A and Group B figures (same canvas, fonts, colours, line
% styles, ticks and legend conventions).  The proposed curve is drawn last and
% the markers of the five algorithms are staggered in alpha, so TSITTCGA and
% HTTCGP remain distinguishable where their profiles coincide.
% Each profile case is one problem--dimension pair. Its cost is the mean
% over the saved random starts (five in the formal experiment). A solver
% receives an infinite cost if any of those attempts failed; successful
% starts are never averaged selectively for a partially failed pair.

here=fileparts(mfilename('fullpath'));
if nargin<1 || isempty(outdir), outdir=fullfile(here,'results'); end
root=fileparts(here);
addpath(root,'-end');                     % profile_panel.m
assert(isfile(fullfile(root,'profile_panel.m')),'plot_results:MissingStyle', ...
    'profile_panel.m was not found in %s.',root);

d=load(fullfile(outdir,'compare.mat'));
folder=fullfile(fileparts(outdir),'figures');
if ~isfolder(folder), mkdir(folder); end
master={'TSITTCGA','AT2DFP','ADCGPM','TSIM-DFPM-PRP','HTTCGP'};
names=d.names; nm=numel(names);
assert(nm==numel(master) && all(ismember(master,names)), ...
    'plot_results:UnknownAlgorithms','Unexpected Group C algorithm list.');
order=cellfun(@(q)find(strcmp(names,q),1),master);   % legend order
names=names(order);
display={'TSITTCGA (this paper)','AT2DFP [15]','ADCGPM [16]','TSIM-DFPM-PRP [29]','HTTCGP [10]'};
close(findall(groot,'Type','figure','Tag','TSITTCGA_large_profiles'));
f=figure('Visible','off','Color','w','Units','centimeters','Position',[2 2 16 5.9], ...
    'Name','Group C: complete algorithms','NumberTitle','off','Tag','TSITTCGA_large_profiles');
tl=tiledlayout(f,1,3,'TileSpacing','compact','Padding','compact');
fields={'Itr','NF','CPU'};
titles={'(a) Itr','(b) NF','(c) Tcpu'};
np=numel(d.inst); nrep=size(d.Success,2);
assert(size(d.Success,1)==np && size(d.Success,3)==nm, ...
    'plot_results:ShapeMismatch','Unexpected success-array dimensions.');
success=reshape(all(d.Success,2),np,nm); success=success(:,order);
fprintf('Group C profiles: %d problem-dimension pairs, mean over %d starts per pair.\n',np,nrep);
fprintf('Group C complete-pair success fractions %s\n',mat2str(mean(success,1),4));
for m=1:3
    ax=nexttile; hold(ax,'on');
    runs=d.(fields{m});
    assert(isequal(size(runs),size(d.Success)), ...
        'plot_results:ShapeMismatch','Unexpected %s dimensions.',fields{m});
    assert(all(isfinite(runs(d.Success))) && all(runs(d.Success)>0), ...
        'plot_results:InvalidSuccessfulCost','A successful run has an invalid %s cost.',fields{m});
    cost=reshape(mean(runs,2),np,nm); cost=cost(:,order);
    assert(isequal(size(cost),size(success)),'Inconsistent Group C arrays.');
    cost(~success | ~isfinite(cost) | cost<=0)=Inf;
    best=min(cost,[],2); ratios=cost./best; ratios(~isfinite(ratios))=Inf;
    h=profile_panel(ax,ratios,display,titles{m},m==1);
    legend(ax,'off');
    if m==1, firstHandle=h; end
    fprintf('Group C %s best fractions %s\n',fields{m}, ...
        mat2str(mean(abs(ratios-1)<1e-12,1),4));
end
figure_frame(tl,firstHandle,display,sprintf('Group C: five algorithms, %d problem-dimension pairs',np),3);
exportOptions={'Units','centimeters','Width',f.Position(3),'Height',f.Position(4), ...
    'Padding','figure','PreserveAspectRatio','on'};
exportgraphics(f,fullfile(folder,'large_profiles.pdf'),'ContentType','vector', ...
    'BackgroundColor','white',exportOptions{:});
exportgraphics(f,fullfile(folder,'large_profiles.png'),'Resolution',300, ...
    'BackgroundColor','white',exportOptions{:});
figs=f;
if usejava('desktop'), f.Visible='on'; else, close(f); end
end


