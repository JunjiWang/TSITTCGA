function [Tall,Tsummary,Tcomparison] = make_tables(outdir)
% Rebuild all published CSV summaries from the saved individual runs.
% Mean* requires all attempts in a problem--dimension row to succeed.
% Successful* columns describe each solver's own successful attempts only.
% CommonMean* uses the same all-solver successful cases for every solver.
% SuccessRate is a fraction in [0,1]; Attempts is the relevant denominator.
if nargin<1 || isempty(outdir)
    outdir=fullfile(fileparts(mfilename('fullpath')),'results');
end
datafile=fullfile(outdir,'compare.mat');
assert(isfile(datafile),'make_tables:MissingResults','Missing results: %s.',datafile);
required={'Itr','NF','CPU','G','Feasibility','Success','Status','names','inst'};
saved=load(datafile,required{:});
assert(all(isfield(saved,required)),'make_tables:IncompleteResults', ...
    'compare.mat does not contain all required run data.');
np=numel(saved.inst); nrep=size(saved.Itr,2); nm=numel(saved.names);
for field={'Itr','NF','CPU','G','Feasibility','Success','Status'}
    value=saved.(field{1});
    assert(size(value,1)==np && size(value,2)==nrep && size(value,3)==nm ...
        && numel(value)==np*nrep*nm,'make_tables:ShapeMismatch', ...
        'Unexpected dimensions for %s.',field{1});
end
assert(islogical(saved.Success),'make_tables:InvalidSuccess', ...
    'Success must be the logical, independently validated success array.');
assert(all(isfinite(saved.Itr(saved.Success))) ...
    && all(isfinite(saved.NF(saved.Success))) ...
    && all(isfinite(saved.CPU(saved.Success))),'make_tables:InvalidSuccess', ...
    'A successful run has nonfinite iteration, function-evaluation or time data.');

rows=cell(np*nrep*nm,10); q=0;
for i=1:np
    for r=1:nrep
        % Retain the same solver order as in the experiment driver.
        for j=mod(r-1+(0:nm-1),nm)+1
            q=q+1;
            rows(q,:)={saved.inst(i).name,saved.inst(i).n,r,saved.names{j}, ...
                saved.Itr(i,r,j),saved.NF(i,r,j),saved.CPU(i,r,j), ...
                saved.G(i,r,j),saved.Feasibility(i,r,j),saved.Status{i,r,j}};
        end
    end
end
Tall=cell2table(rows,'VariableNames', ...
    {'Problem','n','Start','Algorithm','Itr','NF','Time','NormF','Feasibility','Status'});

rows=cell(np*nm,19); q=0;
for i=1:np
    for j=1:nm
        q=q+1;
        ok=reshape(saved.Success(i,:,j),[],1);
        metrics=[reshape(saved.Itr(i,:,j),[],1), ...
            reshape(saved.NF(i,:,j),[],1),reshape(saved.CPU(i,:,j),[],1)];
        complete_mean=nan(1,3);
        if all(ok), complete_mean=mean(metrics,1); end
        successful_mean=nan(1,3); successful_std=nan(1,3); successful_median=nan(1,3);
        if any(ok)
            successful_mean=mean(metrics(ok,:),1);
            successful_median=median(metrics(ok,:),1);
            if nnz(ok)>1, successful_std=std(metrics(ok,:),0,1); end
        end
        rows(q,:)=[{saved.inst(i).name,saved.inst(i).n,saved.names{j}}, ...
            num2cell(complete_mean),{max(saved.G(i,:,j)),nnz(ok),nrep,nnz(ok)/nrep}, ...
            num2cell(successful_mean),num2cell(successful_std),num2cell(successful_median)];
    end
end
Tsummary=cell2table(rows,'VariableNames', ...
    {'Problem','n','Algorithm','MeanItr','MeanNF','MeanTime','MaxNormF', ...
    'Solved','Attempts','SuccessRate', ...
    'SuccessfulMeanItr','SuccessfulMeanNF','SuccessfulMeanTime', ...
    'SuccessfulStdItr','SuccessfulStdNF','SuccessfulStdTime', ...
    'SuccessfulMedianItr','SuccessfulMedianNF','SuccessfulMedianTime'});

success=reshape(saved.Success,np*nrep,nm);
common=all(success,2);
metrics={reshape(saved.Itr,np*nrep,nm),reshape(saved.NF,np*nrep,nm), ...
    reshape(saved.CPU,np*nrep,nm)};
common_mean=nan(nm,3);
if any(common)
    for j=1:3, common_mean(:,j)=mean(metrics{j}(common,:),1)'; end
end
attempts=np*nrep;
Tcomparison=table(string(saved.names(:)),sum(success,1)', ...
    common_mean(:,1),common_mean(:,2),common_mean(:,3), ...
    repmat(attempts,nm,1),sum(success,1)'/attempts,repmat(nnz(common),nm,1), ...
    'VariableNames',{'Algorithm','Solved','CommonMeanItr','CommonMeanNF', ...
    'CommonMeanTime','Attempts','SuccessRate','CommonCases'});

writetable(Tall,fullfile(outdir,'all_runs.csv'),'Encoding','UTF-8');
writetable(Tsummary,fullfile(outdir,'summary.csv'),'Encoding','UTF-8');
writetable(Tcomparison,fullfile(outdir,'comparison_summary.csv'),'Encoding','UTF-8');
end
