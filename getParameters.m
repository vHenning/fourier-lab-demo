function [N, M, Q, L, l] = getParameters(apertureLength, wavelength, distance)
    % Find adequate sampling values
    % We use figure 5.5 (b) from Goodman's Intro to Fourier Optics for this.
    % We assume log(Q(log(NF)) and log(M(log(NF)) are linear.
    fresnelNumber = (apertureLength / 2)^2 / (wavelength * distance);

    slopeQ = -0.6;
    offsetQ = 2.1;
    
    slopeM = 0.4;
    offsetM = 2.8;
    
    % log(Q) = slopeQ * log(fresnelNumber) + offsetQ gives us:
    Q = power(10, offsetQ) * power(fresnelNumber, slopeQ);
    % And the same for M
    M = power(10, offsetM) * power(fresnelNumber, slopeM);
    
    Q = floor(Q);                                      % Ratio of total number of samples and samples in aperture [-]
    M = floor(M);                                      % Number of samples in aperture [-]
    N = Q * M;         % total number of samples [-]
    L = apertureLength * N / M; % total length of aperture (including zero padding) [mm]
    l = apertureLength;
end
