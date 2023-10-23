% Declare our dimensions
wavelength = 1;
apertureLength = 100;
distance1 = 100;
distance2 = 10000;

% Get our sampling values
[N1, M1, Q1, L1, l1] = getParameters(2, apertureLength, wavelength, distance1);
[N2, M2, Q2, L2, l2] = getParameters(2, apertureLength, wavelength, distance2);

% Aperture 1
disp("Aperture 1 Distance 1");
aperture = getAperture1(N1, M1);
plotDiffraction(aperture, L1, wavelength, distance1, "Aperture 1, z = 100", "aperture2D11.png");

disp("Aperture 1 Distance 2");
aperture = getAperture1(N2, M2);
plotDiffraction(aperture, L2, wavelength, distance2, "Aperture 1, z = 1000", "aperture2D12.png");

% Aperture 2
disp("Aperture 2 Distance 1");
aperture = getAperture2(N1, M1);
plotDiffraction(aperture, L1, wavelength, distance1, "Aperture 2, z = 100", "aperture2D21.png");

disp("Aperture 2 Distance 2");
aperture = getAperture2(N2, M2);
plotDiffraction(aperture, L2, wavelength, distance2, "Aperture 2, z = 1000", "aperture2D22.png");

% Aperture 3
disp("Aperture 3 Distance 1");
aperture = getAperture3(N1, M1, L1/N1);
plotDiffraction(aperture, L1, wavelength, distance1, "Aperture 3, z = 100", "aperture2D31.png");

disp("Aperture 3 Distance 2");
aperture = getAperture3(N2, M2, L2/N2);
plotDiffraction(aperture, L2, wavelength, distance2, "Aperture 3, z = 1000", "aperture2D32.png");

% Aperture 4
disp("Aperture 4 Distance 1");
aperture = getAperture4(N1, M1, L1/N1);
plotDiffraction(aperture, L1, wavelength, distance1, "Aperture 4, z = 100", "aperture2D41.png");

disp("Aperture 4 Distance 2");
aperture = getAperture4(N2, M2, L2/N2);
plotDiffraction(aperture, L2, wavelength, distance2, "Aperture 4, z = 1000", "aperture2D42.png");

% Helper functions for plotting the diffraction pattern
function plotDiffraction(aperture, length, wavelength, distance, title, name)
    [x, y, uz, Iz] = diffraction(aperture, length, wavelength, distance);
    Iz = cutZeros(Iz, 0.35);
    savePlot(Iz, title, name)
end

function intensity = cutZeros(intensity, percent)
    samples = size(intensity, 2);
    intensity = intensity(floor(samples * percent):floor(samples * (1-percent)), floor(samples * percent):floor(samples * (1-percent)));
end

function savePlot(intensity, plotTitle, plotName)
    fig = figure("Visible", "off");
    imshow(intensity);
    title(plotTitle);
    print(plotName, '-dpng', ['-r' '2000']);
    close(fig);
end

% Helper functions that create our apertures
function aperture = getAperture1(N, M)
    zeroPaddingAmount = floor((N - M) / 2);
    aperture = zeros(N);
    aperture(zeroPaddingAmount : N - zeroPaddingAmount, zeroPaddingAmount : N - zeroPaddingAmount) = 1;
end

function aperture = getAperture2(N, M)
    % Create completely opaque aperture
    aperture = zeros(N);

    % Make center transparent
    halfN = floor(N/2);
    for i = 1:N
        for j = 1:N
            iShift = i - halfN;
            jShift = j - halfN;

            % If radius is smaller than M we want it to be transparent
            if (sqrt(iShift^2 + jShift^2) < M)
                aperture(i, j) = 1;
            end
        end
    end
end

function aperture = getAperture3(N, M, deltaX)
    % Create completely opaque aperture
    aperture = zeros(N);

    % Make center transparent
    halfN = floor(N/2);
    for i = 1:N
        for j = 1:N
            iShift = i - halfN;
            jShift = j - halfN;

            % If radius is smaller than M we want it to be transparent
            if (sqrt(iShift^2 + jShift^2) < M)
                x = iShift * deltaX;
                aperture(i, j) = 0.5 + 0.5 * cos(2 * pi * x / 10);
            end
        end
    end
end

function aperture = getAperture4(N, M, deltaX)
    % Create completely opaque aperture
    aperture = zeros(N);

    % Make center transparent
    halfN = floor(N/2);
    for i = 1:N
        for j = 1:N
            iShift = i - halfN;
            jShift = j - halfN;

            % If radius is smaller than M we want it to be transparent
            if (sqrt(iShift^2 + jShift^2) < M)
                x = iShift * deltaX;
                aperture(i, j) = exp(1j * pi / 2) * sign(cos(2 * pi * x / 10));
            end
        end
    end
end
