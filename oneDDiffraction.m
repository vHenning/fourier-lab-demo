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

    % Create a mask. We only need the transfer function if the radius of
    % the frequency is smaller than 1/wavelength
    deltaF = 1 / totalLength;
    mask = ones(1, size(transferFunction, 2));
    for i = 1:size(mask, 2)
        freq = (k(i) - (samples/2)) * deltaF;
        if abs(freq) > 1 / wavelength
            mask(i) = 0;
        end
    end

    % Multiply element wise with our shifted angular spectrum
    convolutedShifted = angSpecShifted .* (transferFunction .* mask);
    
    % Shift our angular spectrum back so DC is at the beginning and do an
    % inverse fourier transform
    convoluted = ifftshift(convolutedShifted);
    uz = ifft(convoluted);
    
    Iz = uz .* conj(uz);

    deltaX = totalLength / samples;
    x = 0:deltaX:totalLength;
    x = x(1:samples);
end