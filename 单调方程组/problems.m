function out=problems(x,mode,id)
x=x(:);
assert(numel(x)>=9,'These benchmark definitions are verified for n>=9.');
assert(mode==1 || mode==2);
switch id
    case 1, out=p1(x,mode);
    case 2, out=p2(x,mode);
    case 3, out=p3(x,mode);
    case 4, out=p4(x,mode);
    case 5, out=p5(x,mode);
    case 6, out=p6(x,mode);
    case 7, out=p7(x,mode);
    case 8, out=p8(x,mode);
    otherwise, error('Unknown problem id.');
end
end

function out = p1(x, mode)
n = numel(x);
if mode == 1
    out = exp(x) - 1;
else
    out = max(x, 0);
end
end

function out = p2(x, mode)
n = numel(x);
if mode == 1
    % Vectorised form of
    %   F_1 = x_1 - exp(cos((x_1+x_2)/(n+1))),
    %   F_i = x_i - exp(cos((x_{i-1}+x_i+x_{i+1})/(n+1))),
    %   F_n = x_n - exp(cos((x_{n-1}+x_n)/(n+1))).
    % The elementwise sums use the same left-to-right order, so the values
    % are identical to the scalar loop.
    s = x;
    s(1) = x(1) + x(2);
    s(2:n-1) = x(1:n-2) + x(2:n-1) + x(3:n);
    s(n) = x(n-1) + x(n);
    out = x - exp(cos(s/(n+1)));
else
    out = max(x, 0);
end
end

function out = p3(x, mode)
n = numel(x);
if mode == 1
    out = x - sin(abs(x-1));
else
    v=x+1;
    z=max(v,0);
    if sum(z)>2*n
        u=sort(v,'descend');
        css=cumsum(u)-2*n;
        active=find(u-css./(1:n)'>0,1,'last');
        theta=css(active)/active;
        z=max(v-theta,0);
    end
    out=z-1;
end
end

function out=p4(x,mode)
if mode==2, out=max(x,0); return; end
out=2.5*x-1;
out(1:end-1)=out(1:end-1)+x(2:end);
out(2:end)=out(2:end)+x(1:end-1);
end

function out=p5(x,mode)
if mode==2, out=max(x,0); return; end
q=x(1:end-1).^2+x(2:end).^2;
out=-ones(size(x));
out(1:end-1)=out(1:end-1)+x(1:end-1).*q;
out(2:end)=out(2:end)+x(2:end).*q;
end

function out = p6(x, mode)
% Jian, Ren, Yin and Jiang, JOTA 209:26 (2026), Problem 5.7.
% DOI: 10.1007/s10957-026-02970-y. The original constraint is R_+^n.
% This globally monotone mapping has the unique solution x=0.
if mode == 1
    out = exp(x) - 1;
    out(2:end) = out(2:end) + x(2:end);
else
    out = max(x, 0);
end
end

function out = p7(x, mode)
if mode == 1
    out = 2*x - sin(abs(x));
else
    out = max(x, 0);
end
end

function out = p8(x, mode)
if mode == 1
    out=2*x+expm1(x);
    out(1:end-1)=out(1:end-1)-x(2:end);
    out(2:end)=out(2:end)-x(1:end-1);
else
    out=max(x,0);
end
end


