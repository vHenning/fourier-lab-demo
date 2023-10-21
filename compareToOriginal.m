% Declare our dimensions
wavelength = 1;
apertureLength = 5;
distance1 = 100;

[N1, M1, Q1, L1, l1] = getParameters(apertureLength, wavelength, distance1);

N1 = 16 * M1;

aperture11 = zeros(1, N1);
aperture11(floor((N1/2)-(M1/2)):floor((N1/2)+(M1/2))) = 1;

[x11, uz11, Iz11] = oneDDiffraction(aperture11, N1, wavelength, distance1);

plot(Iz11);