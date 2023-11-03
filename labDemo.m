wavelength = 633e-6;
distance = 10000;
apertureSize = 10;

[N, M, Q, L, l] = getParameters(2, apertureSize, wavelength, distance);

% aperture = getSingleSideBlock(N);
% aperture = aperture';
% aperture = getRectangle(N, M);
aperture = getTriangle(L, l, N);
% aperture = getPentagon(N, M);
% aperture = getHexagon(N, M);

[x, y, uz, Iz] = diffraction(aperture, L, wavelength, distance);

cutAmountPerSide = floor((N-M) * 0.95 / 2);
apertureCut = aperture(cutAmountPerSide:N-cutAmountPerSide, cutAmountPerSide:N-cutAmountPerSide);
IzCut = Iz(cutAmountPerSide:N-cutAmountPerSide, cutAmountPerSide:N-cutAmountPerSide);


figure(1);
imshow(apertureCut);
figure(2);
imshow(IzCut);

function aperture = getStackedSquares(N, M)
    aperture = zeros(N);
    aperture(floor(N/2) - floor(M/2):floor(N/2), floor(N/2) - floor(M/2):floor(N/2) + floor(M/2)) = 1;
    
    aperture(floor(N/2) - M:floor(N/2) - floor(M/2), floor(N/2) - floor(M/4):floor(N/2) + floor(M/4)) = 1;
    aperture(floor(N/2) - M - floor(M/4):floor(N/2) - M, floor(N/2) - floor(M/8):floor(N/2) + floor(M/8)) = 1;
end

function aperture = getSingleSideBlock(N)
    left = zeros(N/2, N);
    right = ones(N/2, N);
    aperture = [left; right];
end

function aperture = getRectangle(N, M)
    aperture = zeros(N);
    aperture(floor(N/2 - M/2):floor(N/2 + M/2), floor(N/2 - M/2):floor(N/2 + M/2)) = 1;
end

function aperture = getPentagon(N, M)
    radius = M;
    numSides = 5;
    theta = [linspace(0, 2*pi, numSides + 1)];
    x = radius * cos(theta);
    x = x + N/2;
    y = radius * sin(theta);
    y = y + N/2;
    aperture = poly2mask(x, y, N, N);
end

function aperture = getHexagon(N, M)
    m = 2^13/2;
    r = 300;
    x = [m+r, m+r*cosd(60), m + r * cosd(120), m+r*cosd(180), m+r*cosd(240), m+r*cosd(300)];
    y = [m, m+r*sind(60), m+r*sind(120), m+r*sind(180), m+r*sind(240), m+r*sind(300)];
    aperture = poly2mask(x, y, 2^13, 2^13);
end