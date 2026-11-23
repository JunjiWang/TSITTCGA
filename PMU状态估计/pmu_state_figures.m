function metrics=pmu_state_figures(sys)
%PMU_STATE_FIGURES Compare saved TSITTCGA, Huber-IRLS and WLS estimates.
%   pmu_state_figures            IEEE14 case: full currents, 4 bad components
%   pmu_state_figures('ieee118') same idea on IEEE118
%Produces fig_bad_data (injected versus estimated bad data) and
%fig_state_estimate (true and estimated bus voltage states).
here=fileparts(mfilename('fullpath'));
if nargin<1||isempty(sys),sys='ieee14';end
c=pmu_parameters();
d=load(fullfile(here,'data',['input_' sys '.mat']));
if isfield(d,'instances'),cs=d.instances;else,cs=d.instance_records;end
% Pick the representative case used before: full current set, middle-level
% contamination, third seed.
target=[];
maxcur=max(cellfun(@(z)z.Currents,cs));
fullset=cs(cellfun(@(z)z.Currents==maxcur,cs));
lv=unique(cellfun(@(z)z.BadData,fullset));
lvl=lv(ceil(numel(lv)/2));
for k=1:numel(cs)
    a=cs{k};
    if a.Currents==maxcur && a.BadData==lvl && a.SeedIndex==3
        target=k;break
    end
end
assert(isscalar(target),'PMU:FigureSelection','Expected one full-current, middle-contamination, third-seed case.');
k=target;inst=cs{k};
if isfield(inst,'BadData'),nb=inst.BadData;else,nb=inst.nbad;end
% Use the completed current-profile paired study that supplies the error table.
f=fullfile(here,'results',['wls_paired_' c.profile '.mat']);
assert(isfile(f),'Run the paired WLS accuracy study before drawing this figure: %s',f);
R=load(f,'rows','endpoints','completed');assert(R.completed,'Accuracy study is incomplete.');
rows=struct2table(R.rows);
selected=strcmp(rows.System,sys)&rows.Case==k ...
    &strcmp(rows.Scenario,'SparseBadData');
ix=find(selected&strcmp(rows.Method,'TSITTCGA'));
ii=find(selected&strcmp(rows.Method,'Huber-IRLS'));
iw=find(selected&strcmp(rows.Method,'WLS'));
assert(isscalar(ix)&&isscalar(ii)&&isscalar(iw),'Expected one saved endpoint per method.');
res=R.endpoints{ix};irls=R.endpoints{ii};wls=R.endpoints{iw};
for endpoint={res,irls}
    r=endpoint{1};m=inst.m;
    assert(r.success&&r.residual<=c.target&&numel(r.s)==2*m ...
        &&norm(min(r.s,0))<=c.feasibility_target,'Saved robust endpoint is not feasible and converged.');
    assert(all(isfinite([r.xi(:);r.e(:)])),'Nonfinite saved robust estimate.');
    assert(norm(r.e(:)-(r.s(1:m)-r.s(m+1:end)))<=1e-12,'Saved split and error estimates differ.');
    fit=inst.y-r.e-inst.G*r.xi;
    residual=norm([min(r.s(1:m),c.lambda-fit);min(r.s(m+1:end),c.lambda+fit)]);
    assert(residual<=c.target&&abs(residual-r.residual)<=1e-9,'Saved robust residual mismatch.');
end
assert(wls.success&&all(isfinite(wls.xi)),'Invalid saved WLS estimate.');
n=inst.n;xi=res.xi;xi_irls=irls.xi;xls=wls.xi;
v=inst.true_xi(1:n)+1i*inst.true_xi(n+1:end);
vh=xi(1:n)+1i*xi(n+1:end);
vi=xi_irls(1:n)+1i*xi_irls(n+1:end);
vl=xls(1:n)+1i*xls(n+1:end);
bus=(1:n)';blue=[0 .38 .68];green=[0 .48 .32];orange=[.83 .35 .12];
scope=sprintf('%s, %d current phasors, %d contaminated components',strrep(upper(sys),'IEEE','IEEE '),inst.Currents,nb);
figdir=fullfile(here,'figures');if ~isfolder(figdir),mkdir(figdir);end
format=@(ax,nn) set(ax,'FontName','Times New Roman','FontSize',9.5,'XGrid','on','YGrid','on','GridAlpha',.15,'Box','on','TickDir','in','XTickLabelRotation',0,'XLim',[.5 nn+.5],'LineWidth',1.0);
% State estimation figure
f=figure('Visible','off','Color','w','Units','centimeters','Position',[2 2 16 8.8]);
tl=tiledlayout(f,2,2,'TileSpacing','compact','Padding','compact');
% Title omitted: the LaTeX caption names the figure.
ax=nexttile(tl);hold(ax,'on');
h(1)=plot(ax,bus,abs(v),'k-','LineWidth',1.0);
h(3)=plot(ax,bus,abs(vi),'-.d','Color',green,'MarkerFaceColor','w','MarkerIndices',2:2:n,'MarkerSize',4.2,'LineWidth',1.0);
h(2)=plot(ax,bus,abs(vh),'--s','Color',blue,'MarkerFaceColor',blue,'MarkerIndices',1:2:n,'MarkerSize',3.3,'LineWidth',1.0);
h(4)=plot(ax,bus,abs(vl),':^','Color',orange,'MarkerFaceColor','w','MarkerIndices',1:3:n,'MarkerSize',3.8,'LineWidth',1.0);
xlabel(ax,'Bus');ylabel(ax,'Voltage magnitude (p.u.)');title(ax,'(a) Voltage magnitude','FontWeight','normal','FontSize',10);format(ax,n);
mag=abs([v,vh,vi,vl]);lo=min(mag,[],'all');hi=max(mag,[],'all');
pad=.08*max(hi-lo,eps);ylim(ax,[lo-pad,hi+pad]);
lg=legend(ax,h,{'True voltage','TSITTCGA','Huber–IRLS','WLS'}, ...
    'Orientation','horizontal','NumColumns',4,'Box','off','FontSize',9);lg.Layout.Tile='north';
ax=nexttile(tl);hold(ax,'on');
plot(ax,bus,angle(v)*180/pi,'k-','LineWidth',1.0);
plot(ax,bus,angle(vi)*180/pi,'-.d','Color',green,'MarkerFaceColor','w','MarkerIndices',2:2:n,'MarkerSize',4.2,'LineWidth',1.0);
plot(ax,bus,angle(vh)*180/pi,'--s','Color',blue,'MarkerFaceColor',blue,'MarkerIndices',1:2:n,'MarkerSize',3.3,'LineWidth',1.0);
plot(ax,bus,angle(vl)*180/pi,':^','Color',orange,'MarkerFaceColor','w','MarkerIndices',1:3:n,'MarkerSize',3.8,'LineWidth',1.0);
xlabel(ax,'Bus');ylabel(ax,'Voltage angle (deg)');title(ax,'(b) Voltage angle','FontWeight','normal','FontSize',10);format(ax,n);
phase=angle([v,vh,vi,vl])*180/pi;lo=min(phase,[],'all');hi=max(phase,[],'all');
pad=.08*max(hi-lo,eps);ylim(ax,[lo-pad,hi+pad]);
ax=nexttile(tl);b=bar(ax,bus,[abs(vh-v),abs(vi-v),abs(vl-v)],'grouped','BarWidth',.85,'EdgeColor','none');
b(1).FaceColor=blue;b(2).FaceColor=green;b(3).FaceColor=orange;
xlabel(ax,'Bus');ylabel(ax,'|V_{est}-V_{true}| (p.u.)');title(ax,'(c) Voltage phasor error','FontWeight','normal','FontSize',10);format(ax,n);
ax=nexttile(tl);b=bar(ax,bus,[abs(angle(vh.*conj(v))),abs(angle(vi.*conj(v))),abs(angle(vl.*conj(v)))]*180/pi,'grouped','BarWidth',.85,'EdgeColor','none');
b(1).FaceColor=blue;b(2).FaceColor=green;b(3).FaceColor=orange;
xlabel(ax,'Bus');ylabel(ax,'Absolute angle error (deg)');
title(ax,'(d) Angle error','FontWeight','normal','FontSize',10);format(ax,n);
save_fig(f,figdir,sprintf('fig_state_estimate_%s',sys));
% Bad-data figure
f=figure('Visible','off','Color','w','Units','centimeters','Position',[2 2 16 6.2]);
tl=tiledlayout(f,1,1,'TileSpacing','compact','Padding','compact');
ax=nexttile(tl);hold(ax,'on');idx=(1:inst.m)';
truth=inst.bad(:)/inst.sigma;estimate=res.e(:)/inst.sigma;estimate_irls=irls.e(:)/inst.sigma;
stem(ax,idx,truth,'k','Marker','none','LineWidth',.9);
stem(ax,idx,estimate_irls,'Color',green,'LineStyle','-.','Marker','none','LineWidth',.9);
stem(ax,idx,estimate,'Color',blue,'LineStyle','--','Marker','none','LineWidth',.9);
a=truth~=0;b=abs(estimate)>10;bi=abs(estimate_irls)>10;
plot(ax,idx(a),truth(a),'ko','MarkerFaceColor','w','MarkerSize',5.2,'LineStyle','none','HandleVisibility','off');
plot(ax,idx(bi),estimate_irls(bi),'d','Color',green,'MarkerFaceColor','w','MarkerSize',5.2,'LineStyle','none','HandleVisibility','off');
plot(ax,idx(b),estimate(b),'s','Color',blue,'MarkerFaceColor',blue,'MarkerSize',3.2,'LineStyle','none','HandleVisibility','off');
xlabel(ax,'Measurement component');ylabel(ax,'Gross error / \sigma');
% Title omitted: the LaTeX caption names the figure.
hs(1)=plot(ax,nan,nan,'k-o','MarkerFaceColor','w','MarkerSize',5.2,'LineWidth',.9);
hs(2)=plot(ax,nan,nan,'--s','Color',blue,'MarkerFaceColor',blue,'MarkerSize',3.2,'LineWidth',.9);
hs(3)=plot(ax,nan,nan,'-.d','Color',green,'MarkerFaceColor','w','MarkerSize',5.2,'LineWidth',.9);
lg=legend(ax,hs,{'Injected bad data','TSITTCGA','Huber–IRLS'}, ...
    'Orientation','horizontal','NumColumns',3,'Box','off','FontSize',9);
lg.Layout.Tile='north';
set(ax,'FontName','Times New Roman','FontSize',9.5,'XGrid','on','YGrid','on','GridAlpha',.15,'Box','on','TickDir','in','LineWidth',1.0);
xlim(ax,[.5 inst.m+.5]);
lo=min([truth;estimate;estimate_irls]);hi=max([truth;estimate;estimate_irls]);
pad=.08*max(hi-lo,eps);ylim(ax,[lo-pad,hi+pad]);
save_fig(f,figdir,sprintf('fig_bad_data_%s',sys));
phasor_ts=abs(vh-v);phasor_ls=abs(vl-v);
angle_ts=abs(angle(vh.*conj(v)))*180/pi;angle_ls=abs(angle(vl.*conj(v)))*180/pi;
metrics=table({sys},k,inst.Currents,nb,100*norm(xi-inst.true_xi)/norm(inst.true_xi), ...
    100*norm(xls-inst.true_xi)/norm(inst.true_xi),nnz(phasor_ts<phasor_ls), ...
    nnz(angle_ts<angle_ls),n,nnz(abs(estimate)>10), ...
    'VariableNames',{'System','Case','Currents','BadData','TSRelativeStateErrorPercent', ...
    'WLSRelativeStateErrorPercent','LowerPhasorErrorBuses','LowerAngleErrorBuses', ...
    'Buses','EstimatedComponentsAbove10Sigma'});
metrics.IRLSRelativeStateErrorPercent=100*norm(xi_irls-inst.true_xi)/norm(inst.true_xi);
metrics.IRLSLowerPhasorErrorBuses=nnz(abs(vi-v)<phasor_ls);
metrics.IRLSLowerAngleErrorBuses=nnz(abs(angle(vi.*conj(v)))*180/pi<angle_ls);
metrics.IRLSEstimatedComponentsAbove10Sigma=nnz(abs(estimate_irls)>10);
metrics.MaxAbsTSIRLSStateDifference=max(abs(xi-xi_irls));
metrics.MaxAbsTSIRLSGrossErrorDifference=max(abs(res.e-irls.e));
writetable(metrics,fullfile(figdir,sprintf('fig_state_metrics_%s.csv',sys)));
fprintf('case %d of %s, %d bad components; phasor/angle improvements %d/%d and %d/%d\n', ...
    k,sys,nb,metrics.LowerPhasorErrorBuses,n,metrics.LowerAngleErrorBuses,n);
end
function save_fig(f,out,name)
set(findall(f,'Type','axes'),'FontName','Times New Roman','FontSize',9.5,'LineWidth',.7);
exportgraphics(f,fullfile(out,[name '.png']),'Resolution',300,'BackgroundColor','white');
exportgraphics(f,fullfile(out,[name '.pdf']),'ContentType','vector','BackgroundColor','white');
close(f);
end

