function out=pmu_map(s,mode,inst,lam)
%PMU_MAP  Reduced optimality system of the PMU state-estimation problem.
%
%   F(s) = min{ s , K s + c }  with  s = (p,q),  e = p - q,
%   K = [P -P; -P P],  c = [lam*1 - P*y ; lam*1 + P*y],
%   and P the orthogonal projector onto range(G)^perp.
m=inst.m;
if mode==2
    out=max(s,0);
    return
end
p=s(1:m); q=s(m+1:end); e=p-q;
Pe=inst.proj(e);
Py=inst.Py;
out=[min(p, Pe+lam-Py); min(q, -Pe+lam+Py)];
end
