% EELE 581 1D Diffraction
% We want to calculate the optical field behind a 1-dimensional aperture.
% We do so using the angular spectrum/transfer function approach.
clear;
% Define our parameters
apertureLength = 5;       % l [mm]
distance = 100;           % z [mm]
wavelength = 1;           % lambda [mm (nm x10^-6)]

% Find adequate sampling values
% We use figure 5.5 (b) from Goodman's Intro to Fourier Optics for this.
% We assume log(Q(log(NF)) and log(M(log(NF)) are linear.
fresnelNumber = (apertureLength / 2)^2 / (wavelength * distance);

slopeQ = -0.6;
offsetQ = 2.1;

slopeM = 0.4;
offsetM = 2.8;

% log(Q) = slopeQ * log(fresnelNumber) + offsetQ gives us:
Q = power(10, offsetQ) * power(fresnelNumber, slopeQ)
% And the same for M
M = power(10, offsetM) * power(fresnelNumber, slopeM)

% zeroPaddingRatio = 16;     % Q [-]
% samplesPerAperture = 1000; % M [-]
zeroPaddingRatio = floor(Q);
samplesPerAperture = floor(M);
Q = 16;
samples = 2^16;
M = floor(samples / Q);
% samples = zeroPaddingRatio * samplesPerAperture;           % N [-]
% totalWidth = apertureLength * samples / samplesPerAperture; % L [mm]

% Make the aperture function. This is just a rect.
aperture = zeros(1, samples);
aperture((samples / 2) - (samplesPerAperture / 2) : (samples / 2) + (samplesPerAperture / 2)) = 1;

% Calculate the discrete fourier transform of the aperture and shift it so
% the DC part is in the middle. We could also shift the transfer function,
% however, this approach was chosen to be consistent with the transfer
% function formula used in Goodman's Intro to Fourier Optics.
angSpec = fft(aperture);
angSpecShifted = fftshift(angSpec);

% Create our discrete frequency steps for fx (k). This is a 1D example so
% we only need fx.
k = 1:samples;

% Create our transfer function as described in Goodman's Intro to Fourier
% Optics
transferFunction = exp(1j * ((2 * pi * distance)/wavelength) * sqrt(1 - (wavelength / totalWidth)^2 * ((k - (samples / 2)).^2)));

% Multiply element wise with our shifted angular spectrum
convolutedShifted = angSpecShifted .* transferFunction;

% Shift our angular spectrum back so DC is at the beginning and do an
% inverse fourier transform
convoluted = ifftshift(convolutedShifted);
pattern = ifft(convoluted);

intensity = pattern .* conj(pattern);

figure(1);
plot(phase(transferFunction));
figure(2);
plot(abs(convolutedShifted));
figure(3);
plot(intensity);
figure(4);
plot(aperture);
