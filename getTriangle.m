% Get a triangle shaped aperture.
% l: triangle side length
% L: total length including zero padding
% N: total number of samples per 
function aperture = getTriangle(L, l, N)
    aperture = zeros(N);
    
    pixelWidth = floor(N * l / L);

    pixelHeight = floor(N * sqrt(3/4) * l / L);

    pixelOffset = pixelHeight / 2;
    pixelSlope = pixelHeight / (pixelWidth / 2);

    for it = floor((N/2) - pixelWidth / 2):floor(N/2 + pixelWidth/2)
        x = it - floor(N/2);
        if x < 0
            mask = pixelSlope * x + pixelOffset;
        else
            mask = -pixelSlope * x + pixelOffset;
        end
        for jt = floor(N/2 - pixelHeight/2):floor(N/2 + pixelHeight/2)
            y = jt - floor(N/2);
            if y < mask
                aperture(it, jt) = 1;
            end
        end
    end
end