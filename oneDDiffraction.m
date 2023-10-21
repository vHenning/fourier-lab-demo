function [x, uz, Iz] = oneDDiffraction(aperture, totalLength, wavelength, distance)
    % Calculate the discrete fourier transform of the aperture and shift it so
    % the DC part is in the middle. We could also shift the transfer function,
    % however, this approach was chosen to be consistent with the transfer
    % function formula used in Goodman's Intro to Fourier Optics.
    angSpec = fft(aperture);
    angSpecShifted = fftshift(angSpec);
    
    % Create our discrete frequency steps for fx (k). This is a 1D example so
    % we only need fx.
    samples = size(aperture, 2);
    k = 1:samples;
    
    % Create our transfer function as described in Goodman's Intro to Fourier
    % Optics
    transferFunction = exp(1j * ((2 * pi * distance)/wavelength) * sqrt(1 - (wavelength / totalLength)^2 * ((k - samples / 2)).^2));

    % Multiply element wise with our shifted angular spectrum
    convolutedShifted = angSpecShifted .* transferFunction;
    
    % Shift our angular spectrum back so DC is at the beginning and do an
    % inverse fourier transform
    convoluted = ifftshift(convolutedShifted);
    uz = ifft(convoluted);
    
    Iz = uz .* conj(uz);

    deltaX = totalLength / samples;
    x = 0:deltaX:totalLength;
    x = x(1:samples);
end