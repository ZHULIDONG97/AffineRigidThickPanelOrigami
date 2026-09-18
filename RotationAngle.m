function Rotangle = RotationAngle(nei,iRotation,x,y,z)
% Transformation from nodal coordinates to signed rotation angle.
Rotangle = zeros(nei,1);
for i = 1:nei
    i1 = iRotation(i,1);
    i2 = iRotation(i,2);
    i3 = iRotation(i,3);
    i4 = iRotation(i,4);

    x1 = x(i1);  y1 = y(i1);  z1 = z(i1);
    x2 = x(i2);  y2 = y(i2);  z2 = z(i2);
    x3 = x(i3);  y3 = y(i3);  z3 = z(i3);
    x4 = x(i4);  y4 = y(i4);  z4 = z(i4);

    rij = [x1-x2,y1-y2,z1-z2];
    rkj = [x3-x2,y3-y2,z3-z2];
    rkl = [x3-x4,y3-y4,z3-z4];

    m = cross(rkj,rij);
    n = cross(rkj,rkl);

    lrkj = norm(rkj);
    if lrkj == 0
        error('RotationAngle:DegenerateCrease', ...
            'Rotation row %d has zero crease length.', i);
    end

    signedSine = dot(cross(n,m),rkj) / lrkj;
    Rotangle(i) = atan2(signedSine,dot(n,m));
end
end
