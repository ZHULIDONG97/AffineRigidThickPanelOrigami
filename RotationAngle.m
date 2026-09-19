function Rotangle = RotationAngle(nei,iRotation,x,y,z)
% Transformation from nodal coordinates to signed rotation angle.
Rotangle = zeros(nei,1);
for i = 1:nei
    i1 = iRotation(i,1);
    i2 = iRotation(i,2);
    i3 = iRotation(i,3);
    i4 = iRotation(i,4);

    rij = [x(i1)-x(i2),y(i1)-y(i2),z(i1)-z(i2)];
    rkj = [x(i3)-x(i2),y(i3)-y(i2),z(i3)-z(i2)];
    rkl = [x(i3)-x(i4),y(i3)-y(i4),z(i3)-z(i4)];

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
