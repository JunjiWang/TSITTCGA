function S=summarize_pmu(matfile)
% Arithmetic mean timing over repeats and cases; Itr/NF use common successes.
d=load(matfile,'rows','protocol','completed');T=struct2table(d.rows);rows=struct([]);
for sy=d.protocol.Config.systems
    for method=d.protocol.Config.names
        q=T(strcmp(T.System,sy{1})&strcmp(T.Method,method{1}),:);
        if isempty(q),continue;end
        ids=unique(q.Case);ok=false(numel(ids),1);tm=nan(size(ok));sol=tm;it=tm;nf=tm;er=tm;
        for j=1:numel(ids)
            z=q(q.Case==ids(j),:);ok(j)=all(z.Success);tm(j)=mean(z.TotalTime);
            sol(j)=mean(z.SolverWallTime);it(j)=z.Itr(1);nf(j)=z.NF(1);er(j)=z.RelativeStateError(1);
            assert(all(z.Success==z.Success(1))&&all(z.Itr==z.Itr(1)), 'Nondeterministic repeat.');
            assert(max(z.Residual)-min(z.Residual)<=1e-9,'Residual repeat mismatch.');
        end
        allmethods=d.protocol.Config.names;common=ids;
        for mm=allmethods
            z=T(strcmp(T.System,sy{1})&strcmp(T.Method,mm{1}),:);
            bad=unique(z.Case(~z.Success));common=setdiff(common,bad);
        end
        joint=ismember(ids,common);
        row=struct('System',sy{1},'Method',method{1},'Solved',nnz(ok),'Cases',numel(ids), ...
            'MeanTotalTimeAll',mean(tm), ...
            'MeanSolverTimeAll',mean(sol),'CommonSolved',nnz(joint), ...
            'MeanTotalTimeCommon',mean(tm(joint)),'MeanItrCommon',mean(it(joint)), ...
            'MeanNFCommon',mean(nf(joint)),'MeanStateErrorSolved',mean(er(ok)), ...
            'MaxResidualAll',max(q.Residual),'Completed',d.completed);
        if isempty(rows),rows=row;else,rows(end+1)=row;end %#ok<AGROW>
    end
end
S=struct2table(rows);[p,n]=fileparts(matfile);writetable(S,fullfile(p,[n '_summary.csv']));
end
