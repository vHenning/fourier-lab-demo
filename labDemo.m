numberSquares = 3;
squareSize = 0.5;
squareDistance = 5;

wavelength = 633e-6;
distance = 10000;
apertureSize = 10;

[N, M, Q, L, l] = getParameters(2, apertureSize, wavelength, distance);

Q = 16;
N = M * Q;
L = l * Q;

aperture = zeros(N);
aperture(floor(N/2) - floor(M/2):floor(N/2) + floor(M/2), floor(N/2) - floor(M/2):floor(N/2) + floor(M/2)) = 1;

aperture(floor(N/2) - M:floor(N/2) - floor(M/2), floor(N/2) - floor(M/4):floor(N/2) + floor(M/4)) = 1;
aperture(floor(N/2) - M - floor(M/4):floor(N/2) - M, floor(N/2) - floor(M/8):floor(N/2) + floor(M/8)) = 1;

imshow(aperture);

