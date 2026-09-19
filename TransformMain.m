function [T, r] = TransformMain(ix, n)
r = length(ix);
constrained = ix(:);
order = [constrained; setdiff((1:n).',constrained,'stable')];
I = eye(n);
T = I(order, :);
end
