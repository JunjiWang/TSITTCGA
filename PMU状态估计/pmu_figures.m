function pmu_figures(which)
%PMU_FIGURES  Dolan-More performance profiles for the PMU comparison.
%   pmu_figures            IEEE14 + IEEE118 (main comparison, six methods)
%   pmu_figures('300')     IEEE300 alone, once its run has finished
here=fileparts(mfilename('fullpath'));oldpath=path;
cleanup=onCleanup(@()path(oldpath));addpath(here,'-begin'); %#ok<NASGU>
if nargin<1,which='main';end
c=pmu_parameters();
switch which
    case 'main'
        systems={'ieee14','ieee118'};
    case '300'
        systems={'ieee300'};
    otherwise
        error('Use main or 300.');
end
sources=pmu_efficiency_sources(systems);
figdir=fullfile(here,'figures');if ~isfolder(figdir),mkdir(figdir);end
names=c.names;
titles={'(a) Itr + 1','(b) Time'};
for si=1:numel(systems)
    sys=systems{si};
    d=load(sources(si).Path,'rows');T=struct2table(d.rows);
    cases=unique(T.Case(strcmp(T.System,sys)));
    cost=inf(numel(cases),numel(names),2);
    for j=1:numel(cases)
        for k=1:numel(names)
            z=T(strcmp(T.System,sys)&T.Case==cases(j)&strcmp(T.Method,names{k}),:);
            if isempty(z)||~all(z.Success),continue;end
            % A method can meet the tolerance at its initial point (Itr=0).
            % Shift every method's iteration cost by one for finite ratios;
            % the raw iteration counts and tabulated statistics are unchanged.
            cost(j,k,1)=mean(z.Itr)+1;
            cost(j,k,2)=mean(z.TotalTime);
        end
    end
    fig=figure('Visible','off','Color','w','Units','centimeters','Position',[2 2 18 6.8]);
    tl=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
    handles=cell(1,2);
    for m=1:2
        ratios=cost(:,:,m)./min(cost(:,:,m),[],2);
        ratios(~isfinite(ratios))=Inf;
        ax=nexttile;hold(ax,'on');
        handles{m}=profile_panel(ax,ratios,names,titles{m},m==1);
        legend(ax,'off');
    end
    lg=legend(handles{1},names,'Orientation','horizontal','NumColumns',6, ...
        'FontName','Times New Roman','FontSize',8,'Box','off');
    lg.Layout.Tile='south';
    title(tl,sprintf('%s: six methods, 45 cases each',sys),'FontWeight','normal','FontSize',9);
    out=fullfile(figdir,['pmu_profiles_' sys '.png']);
    exportgraphics(fig,out,'Resolution',300,'BackgroundColor','white');
    exportgraphics(fig,strrep(out,'.png','.pdf'),'ContentType','vector','BackgroundColor','white');
    close(fig);
    fprintf('saved %s\n',out);
end
end

function handles=profile_panel(ax,ratio,names,panelTitle,showYLabel)
%PROFILE_PANEL  One Dolan-More performance-profile panel in the paper's style.
%
%   handles = profile_panel(ax, ratio, names, panelTitle, showYLabel)
%
%   ratio       n-by-nm matrix of cost ratios (Inf for failed runs)
%   names       cellstr of algorithm labels, in legend order
%   panelTitle  e.g. '(a) Iterations'
%
%   Shared style block (every profile figure of the paper): same canvas,
%   fonts, colours, line width (1.4 pt), marker size (3.8 pt), ticks and
%   legend conventions.
%
%   Each curve is drawn on the exact breakpoints alpha=log2(ratio), so the
%   staircase is the true profile rather than a fixed-step approximation,
%   and the markers are evaluated at the same function.  A curve stops at the
%   alpha where it first reaches rho_s=1 (a constant saturated tail adds no
%   information); algorithms that never reach 1 are drawn to the right edge.

nm=numel(names);
[colors,styles,markers]=profile_style(names);
alphaUpper=profile_upper(ratio);
rho=@(j,a) mean(ratio(:,j)<=2^a*(1+1e-12));
proposed=find(ismember(names,{'TSITTCGA','Proposed','Warm start + increase'}),1);
drawOrder=[setdiff(1:nm,proposed,'stable'),proposed];
markerStep=max(1,alphaUpper/5);
handles=gobjects(1,nm);
for j=drawOrder
    % Exact breakpoints of this curve, clipped to the plotted range.
    rj=ratio(:,j); rj=rj(isfinite(rj) & rj>0);
    brk=unique([0; log2(rj(:))]);
    brk=brk(brk>=0 & brk<=alphaUpper);
    yb=arrayfun(@(a)rho(j,a),brk);
    sat=brk(find(yb>=1-1e-12,1));
    if isempty(sat), sat=alphaUpper; end
    aa=[brk(brk<sat); sat];
    ya=arrayfun(@(a)rho(j,a),aa);
    stairs(ax,aa,ya,'Color',colors(j,:),'LineStyle',styles{j}, ...
        'LineWidth',1.4,'HandleVisibility','off');
    offset=(j-1)*markerStep/(nm+1);
    alphaMark=offset:markerStep:sat;
    if isempty(alphaMark), alphaMark=0; end
    if alphaMark(end)<sat, alphaMark=[alphaMark,sat]; end
    ym=arrayfun(@(a)rho(j,a),alphaMark);
    plot(ax,alphaMark,ym,markers{j},'Color',colors(j,:),'MarkerFaceColor','w', ...
        'MarkerSize',3.8,'LineWidth',.8,'LineStyle','none','HandleVisibility','off');
    handles(j)=plot(ax,NaN,NaN,'Color',colors(j,:),'LineStyle',styles{j}, ...
        'LineWidth',1.4,'Marker',markers{j},'MarkerSize',3.8,'MarkerFaceColor','w');
end
tickStep=max(1,ceil(alphaUpper/6));
set(ax,'XScale','linear','XTick',0:tickStep:alphaUpper,'FontName','Times New Roman', ...
    'FontSize',8.5,'Box','on','TickDir','out','LineWidth',.8,'YTick',0:.25:1, ...
    'GridAlpha',.10,'XMinorGrid','off','YMinorGrid','off');
xlim(ax,[0 alphaUpper]); ylim(ax,[0 1.025]); grid(ax,'on');
pbaspect(ax,[1.2 1 1]);
title(ax,panelTitle,'FontWeight','normal','FontSize',9);
xlabel(ax,'\alpha = log_2(\tau)','Interpreter','tex');
if showYLabel, ylabel(ax,'\rho_s(\alpha)','Interpreter','tex'); end
legend(ax,handles,names,'Location','southeast','Box','off','Color','white', ...
    'FontName','Times New Roman','FontSize',7.2,'Interpreter','none');
end

function upper=profile_upper(ratio)
finite=ratio(isfinite(ratio) & ratio>0);
if isempty(finite)
    upper=1;
else
    upper=ceil(log2(max(2,max(finite)*1.01)));
end
end

function [colors,styles,markers]=profile_style(names)
% Fixed colour / line style / marker for every label used in the paper.
blue=[0 114 178]/255; pink=[204 121 167]/255; orange=[213 94 0]/255;
gray=[70 70 70]/255; green=[0 158 115]/255;
table={ ...
    'TSITTCGA', blue, '-', 'o'; ...
    'AT2DFP', pink, ':', 'd'; ...
    'ADCGPM', orange, '--', 's'; ...
    'TSIM-DFPM-PRP', gray, '--', '>'; ...
    'HTTCGP', green, '-.', '^'; ...
    'Huber-IRLS', [0.49 0.18 0.56], ':', 'p'};
nm=numel(names); colors=zeros(nm,3); styles=cell(1,nm); markers=cell(1,nm);
for j=1:nm
    k=find(strcmp(table(:,1),names{j}),1);
    assert(~isempty(k),'profile_panel:UnknownLabel','Unknown algorithm label: %s.',names{j});
    colors(j,:)=table{k,2}; styles{j}=table{k,3}; markers{j}=table{k,4};
end
end




