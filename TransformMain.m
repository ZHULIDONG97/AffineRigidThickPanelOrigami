function [T, r] = TransformMain(ix, n)

r = length(ix);

constrained = ix(:);
allDofs = (1:n).';

free = setdiff(allDofs, constrained, 'stable');

order = [constrained; free];

I = eye(n);
T = I(order, :);

end