function dist = CalcDist(a, b)

ax = a(1);
ay = a(2);
bx = b(1);
by = b(2);

if length(a) == 2
    dist = sqrt((ax-bx)^2+(ay-by)^2);
elseif length(a) == 3
    az = a(3);
    bz = b(3);
    dist = sqrt((ax-bx)^2+(ay-by)^2+(az-bz)^2);
end

end