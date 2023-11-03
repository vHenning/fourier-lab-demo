wavelength = 633e-6;
distanceNear = 1000;
distanceFar = 10000;
apertureSize = 10;

[N1, M1, Q1, L1, l1] = getParameters(2, apertureSize, wavelength, distanceNear);
[N2, M2, Q2, L2, l2] = getParameters(2, apertureSize, wavelength, distanceFar);

% apertureNear = getSingleSideBlock(N1);
% apertureFar = getSingleSideBlock(N2);
% apertureNear = apertureNear';
% apertureFar = apertureFar';
% apertureNear = getRectangle(N1, M1);
% apertureFar = getRectangle(N2, M2);
% apertureNear = getTriangle(L1, l1, N1);
% apertureFar = getTriangle(L2, l2, N2);
% apertureNear = getPentagon(N1, M1);
% apertureFar = getPentagon(N2, M2);
apertureNear = getHexagon(N1, M1);
apertureFar = getHexagon(N2, M2);


[xNear, yNear, uzNear, IzNear] = diffraction(apertureNear, L1, wavelength, distanceNear);
clear uzNear;
[xFar, yFar, uzFar, IzFar] = diffraction(apertureFar, L2, wavelength, distanceFar);
clear uzFar;

cutAmountPerSide = floor((N1-M1) * 0.95 / 2);
apertureCut = apertureNear(cutAmountPerSide:N1-cutAmountPerSide, cutAmountPerSide:N1-cutAmountPerSide);
IzNearCut = IzNear(cutAmountPerSide:N1-cutAmountPerSide, cutAmountPerSide:N1-cutAmountPerSide);

cutAmountPerSide = floor((N2-M2) * 0.95 / 2);
IzFarCut = IzFar(cutAmountPerSide:N2-cutAmountPerSide, cutAmountPerSide:N2-cutAmountPerSide);

subplot(2,2,1);
imshow(apertureCut);
title('Aperture');
subplot(2,2,2);
imshow(IzNearCut);
title('z = 1000 mm');
subplot(2,2,4);
imshow(IzFarCut);
title('z = 10 000 mm')

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
    m = N/2;
    r = M;
    x = [m+r, m+r*cosd(60), m + r * cosd(120), m+r*cosd(180), m+r*cosd(240), m+r*cosd(300)];
    y = [m, m+r*sind(60), m+r*sind(120), m+r*sind(180), m+r*sind(240), m+r*sind(300)];
    aperture = poly2mask(x, y, N, N);
end