wavelength = 633e-6;
distance = 10000;
apertureSize = 10;

[N, M, Q, L, l] = getParameters(2, apertureSize, wavelength, distance);

aperture = getTriangle(L, l, N);

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