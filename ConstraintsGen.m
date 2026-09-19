function [cineq, ceq, state] = ConstraintsGen(beta,x0,y0,z0,B, ...
    rigidityDirections,rigidityReferences,rigidityOffsets, ...
    Cchi,Cfold,r,T,chiMargin,foldAngleMargin)
% Constraint evaluator: cineq <= 0, and ceq is the rigidity residual.
if nargin < 13 || isempty(chiMargin)
    chiMargin = 0;
end
if nargin < 14 || isempty(foldAngleMargin)
    foldAngleMargin = chiMargin;
end

mappedB = T \ B;
x = x0 + mappedB * beta(1:r);
y = y0 + mappedB * beta(r+1:2*r);
z = z0 + mappedB * beta(2*r+1:3*r);

% Evaluate only the independent A_e'*A_e-I constraints selected at setup.
ceq = AffineMetricResidualJacobian(beta,rigidityDirections, ...
    rigidityReferences,rigidityOffsets);

% The normalized triple product preserves the prescribed panel orientation.
nchi = size(Cchi,1);
chiValue = zeros(nchi,1);
chiLengthScale = zeros(nchi,1);
orientationIneq = zeros(nchi,1);
for i=1:nchi
    ia = Cchi(i,1);
    ib = Cchi(i,2);
    ic = Cchi(i,3);
    id = Cchi(i,4);
    rab = [x(ia)-x(ib); y(ia)-y(ib); z(ia)-z(ib)];
    rbc = [x(ib)-x(ic); y(ib)-y(ic); z(ib)-z(ic)];
    rdb = [x(id)-x(ib); y(id)-y(ib); z(id)-z(ib)];
    rab0 = [x0(ia)-x0(ib); y0(ia)-y0(ib); z0(ia)-z0(ib)];
    rbc0 = [x0(ib)-x0(ic); y0(ib)-y0(ic); z0(ib)-z0(ic)];
    rdb0 = [x0(id)-x0(ib); y0(id)-y0(ib); z0(id)-z0(ib)];
    chiLengthScale(i) = max([norm(rab0), norm(rbc0), norm(rdb0)]);
    chiScale = chiLengthScale(i)^3;
    chiValue(i) = dot(cross(rab,rbc),rdb) / chiScale;
    orientationIneq(i) = chiMargin / chiScale -  chiValue(i);
end

% Cfold rows: [refNodeFaceI creaseStart creaseEnd refNodeFaceJ sign].
nfold = size(Cfold,1);
foldAngle = zeros(nfold,1);
foldIneq = zeros(nfold,1);
for i=1:nfold
    i1 = Cfold(i,1);
    i2 = Cfold(i,2);
    i3 = Cfold(i,3);
    i4 = Cfold(i,4);
    rij = [x(i1)-x(i2), y(i1)-y(i2), z(i1)-z(i2)];
    rkj = [x(i3)-x(i2), y(i3)-y(i2), z(i3)-z(i2)];
    rkl = [x(i3)-x(i4), y(i3)-y(i4), z(i3)-z(i4)];
    m = cross(rkj,rij);
    n = cross(rkj,rkl);
    lrkj = norm(rkj);
    if lrkj <= eps
        error('ConstraintsGen:DegenerateFoldCrease', ...
            'Fold constraint row %d has zero crease length.', i);
    end
    foldAngle(i) = atan2(dot(cross(n,m),rkj) / lrkj,dot(n,m));
    foldIneq(i) = foldAngleMargin - Cfold(i,5) * foldAngle(i);
end

cineq = [orientationIneq; foldIneq];

if nargout > 2
    state.x = x;
    state.y = y;
    state.z = z;
    state.metricResidual = ceq;
    state.orientationChi = chiValue;
    state.orientationLengthScale = chiLengthScale;
    state.orientationInequality = orientationIneq;
    state.foldAngle = foldAngle;
    state.foldInequality = foldIneq;
    state.cineq = cineq;
    state.ceq = ceq;
end
end
