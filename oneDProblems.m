% Declare our dimensions
wavelength = 1;
apertureLength = 100;
distance1 = 100;
distance2 = 10000;

% Get our sampling values
[N1, M1, Q1, L1, l1] = getParameters(apertureLength, wavelength, distance1);
[N2, M2, Q2, L2, l2] = getParameters(apertureLength, wavelength, distance2);

% How long is one step
deltaX1 = L1/N1;
deltaX2 = L2/N2;

% Our x (only the aperture, no zero padding
input1 = 0:deltaX1:apertureLength;
input2 = 0:deltaX2:apertureLength;

% How many samples does one side of zero padding contain
halfZeroPadding1 = floor((N1 - M1) / 2);
halfZeroPadding2 = floor((N2 - M2) / 2);

% Aperture 1
aperture11 = zeros(1, N1);
aperture11(floor((N1/2)-(M1/2)):floor((N1/2)+(M1/2))) = 1;

aperture12 = zeros(1, N2);
aperture12(floor((N2/2)-(M2/2)):floor((N2/2)+(M2/2))) = 1;

% Aperture 2
aperture21NoPadding = 0.5 + 0.5 * cos((2 * pi * input1) / 10);
aperture22NoPadding = 0.5 + 0.5 * cos((2 * pi * input2) / 10);

aperture21 = [zeros(1, halfZeroPadding1), aperture21NoPadding, zeros(1, halfZeroPadding1)];
aperture22 = [zeros(1, halfZeroPadding2), aperture22NoPadding, zeros(1, halfZeroPadding2)];

% Aperture 3
aperture31NoPadding = exp((pi/2) * 1j * sign(cos((2 * pi * input1) / 10)));
aperture32NoPadding = exp((pi/2) * 1j * sign(cos((2 * pi * input2) / 10)));

aperture31 = [zeros(1, halfZeroPadding1), aperture21NoPadding, zeros(1, halfZeroPadding1)];
aperture32 = [zeros(1, halfZeroPadding2), aperture22NoPadding, zeros(1, halfZeroPadding2)];

% Calculate all diffractions
[x11, uz11, Iz11] = oneDDiffraction(aperture11, N1, wavelength, distance1);
[x12, uz12, Iz12] = oneDDiffraction(aperture12, N2, wavelength, distance2);
[x21, uz21, Iz21] = oneDDiffraction(aperture21, N1, wavelength, distance1);
[x22, uz22, Iz22] = oneDDiffraction(aperture22, N2, wavelength, distance2);
[x31, uz31, Iz31] = oneDDiffraction(aperture31, N1, wavelength, distance1);
[x32, uz32, Iz32] = oneDDiffraction(aperture32, N2, wavelength, distance2);

% Plot diffractions. First cut some of the zero padding
cutEdgeAmount = 0.45;
[x11, Iz11] = cutZeros(x11, Iz11, cutEdgeAmount);
[x12, Iz12] = cutZeros(x12, Iz12, cutEdgeAmount);
[x21, Iz21] = cutZeros(x21, Iz21, cutEdgeAmount);
[x22, Iz22] = cutZeros(x22, Iz22, cutEdgeAmount);
[x31, Iz31] = cutZeros(x31, Iz31, cutEdgeAmount);
[x32, Iz32] = cutZeros(x32, Iz32, cutEdgeAmount);

% Save plots as png
plotDiffraction(x11, Iz11, "Aperture 1, z = 100", "aperture11.png");
plotDiffraction(x12, Iz12, "Aperture 1, z = 10000", "aperture12.png");
plotDiffraction(x21, Iz21, "Aperture 2, z = 100", "aperture21.png");
plotDiffraction(x22, Iz22, "Aperture 2, z = 10000", "aperture22.png");
plotDiffraction(x31, Iz31, "Aperture 3, z = 100", "aperture31.png");
plotDiffraction(x32, Iz32, "Aperture 3, z = 10000", "aperture32.png");

function [x, y] = cutZeros(x, y, percent)
    samples = size(x, 2);
    x = x(floor(samples * percent):floor(samples * (1-percent)));
    y = y(floor(samples * percent):floor(samples * (1-percent)));
end

function plotDiffraction(x, y, plotTitle, plotName)
    fig = figure("Visible", "off");
    plot(x, y);
    title(plotTitle);
    print(plotName, '-dpng', ['-r' '2000']);
    close(fig);
end

