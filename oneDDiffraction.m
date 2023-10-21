function [x, y, uz, Iz] = oneDDiffraction(aperture, totalLength, wavelength, distance)
    % Calculate the discrete fourier transform of the aperture and shift it so
    % the DC part is in the middle. We could also shift the transfer function,
    % however, this approach was chosen to be consistent with the transfer
    % function formula used in Goodman's Intro to Fourier Optics.
    angSpec = fft2(aperture);
    angSpecShifted = fftshift(angSpec);
    
    % Create our discrete frequency steps for fx (k). This is a 1D example so
    % we only need fx.
    samplesX = size(aperture, 2);
    samplesY = size(aperture, 1);
    k = 1:samplesX;
    p = (1:samplesY).';

    % Create our transfer function as described in Goodman's Intro to Fourier
    % Optics
    transferFunction = exp(1j * ((2 * pi * distance)/wavelength) * sqrt(1 - (wavelength / totalLength)^2 * ((k - samplesX / 2).^2 + (p - samplesY / 2).^2)));

    % Create a mask. We only need the transfer function if the radius of
    % the frequency is smaller than 1/wavelength
    deltaF = 1 / totalLength;
    mask = ones(size(transferFunction));
    for i = 1:size(mask, 1)
        for j = 1: size(mask, 2)
            freq = sqrt(deltaF * ((k(j) - (samplesX/2))^2 + (p(i) - (samplesY/2))^2));
            if abs(freq) > 1 / wavelength
                mask(i, j) = 0;
            end
        end
    end

    % Multiply element wise with our shifted angular spectrum
    convolutedShifted = angSpecShifted .* (transferFunction .* mask);
    
    % Shift our angular spectrum back so DC is at the beginning and do an
    % inverse fourier transform
    convoluted = ifftshift(convolutedShifted);
    uz = ifft2(convoluted);
    
    Iz = uz .* conj(uz);

    deltaX = totalLength / samplesX;
    x = 0:deltaX:totalLength;
    x = x(1:samplesX);

    deltaY = totalLength / samplesY;
    y = 0:deltaY:totalLength;
    y = y(1:samplesY);
end