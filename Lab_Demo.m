clear all;
close all;
%% Set Aperture Order
n = 3;
    %1: Verticle half plane 
    %2: Horizontal half plane
    %3: Triangle
    %4: Rectangle
    %5: Pentagon
    %6: Hexagon


%% Create transfer function

transferFcn = @(distance,f_x,f_y) ...
    exp(1i*(2*pi*distance) * sqrt(1 - f_x.^2-f_y.^2));

%% Number of Points
N = [9,13];
Q = [10;100]; %vector for Q
NN = 2.^N;
X = Q * 100; %total size of array in the x direction
Y = Q * 100; %total size of array in the y direction

%% Distance Vector
dist = [100;10000];
%% Apertures
switch n
    case 1
        left1 = zeros(NN(1), NN(1)/2);
        left2 = zeros(NN(2), NN(2)/2);
        right1 = ones(NN(1), NN(1)/2);
        right2 = ones(NN(2), NN(2)/2);
        grating_aperture1 = [right1, left1];
        grating_aperture2 = [right2, left2];
    case 2
        up1 = zeros(NN(1)/2, NN(1));
        down1 = ones(NN(1)/2, NN(1));
        up2 = zeros(NN(2)/2, NN(2));
        down2 = ones(NN(2)/2, NN(2));
        grating_aperture1 = [up1; down1];
        grating_aperture2 = [up2; down2];
    case 3

        w = 100;
        h = round(sqrt(3)/2 * w);
        X = Q * w;
        Y = Q * w;

        A = zeros(NN(1),NN(1));
        A(:,(NN(1)/2 - w/2):(NN(1)/2 + w/2)) = 1;
        A(NN(1)/2 + h/2 + 1:end,:) = 0;
        
        triangle1  = zeros(NN(1),NN(1));
        intercept1 = -2*(NN(1)/2 + w/2) + (NN(1)/2 + h/2);
        for row = 1:NN(1)
            colStart =  round(row*2 + intercept1);
            if colStart > 0
                triangle1((colStart:end),row) = 1;
            end
        end
        
        triangle2  = zeros(NN(1),NN(1));
        intercept2 = 2*(NN(1)/2 - w/2) + (NN(1)/2 + h/2);
        for row = 1:NN(1)
            colStart =  round(-row*2 + intercept2);
            if colStart > 0
                triangle2((colStart:end),row) = 1;
            end
        end
        
        grating_aperture1 = A & triangle2 & triangle1;

        AA = zeros(NN(2),NN(2));
        AA(:,(NN(2)/2 - w/2):(NN(2)/2 + w/2)) = 1;
        AA(NN(2)/2 + h/2 + 1:end,:) = 0;
        
        triangle3  = zeros(NN(2),NN(2));
        intercept3 = -2*(NN(2)/2 + w/2) + (NN(2)/2 + h/2);
        for row = 1:NN(2)
            colStart =  round(row*2 + intercept3);
            if colStart > 0
                triangle3((colStart:end),row) = 1;
            end
        end
        
        triangle4  = zeros(NN(2),NN(2));
        intercept4 = 2*(NN(2)/2 - w/2) + (NN(2)/2 + h/2);
        for row = 1:NN(2)
            colStart =  round(-row*2 + intercept4);
            if colStart > 0
                triangle4((colStart:end),row) = 1;
            end
        end
        grating_aperture2 = AA & triangle3 & triangle4;
    case 4
        width = 100;
        height = 50;

        A = zeros(NN(1),NN(1));
        A((NN(1)/2-width/2):(NN(1)/2+width/2), (NN(1)/2-height/2):(NN(1)/2+height/2))=1;
        grating_aperture1  =  A;
        AA = zeros(NN(2),NN(2));
        AA((NN(2)/2-width/2):(NN(2)/2+width/2), (NN(2)/2-height/2):(NN(2)/2+height/2))=1;
        grating_aperture2 = AA;
    case 5 
        radius = 100;
        numSides = 5;
        theta = [linspace(0, 2*pi, numSides + 1)];
        x = radius * cos(theta);
        x1 = x + NN(1)/2;
        x2 = x + NN(2)/2;
        y = radius * sin(theta);
        y1 = y + NN(1)/2;
        y2 = y + NN(2)/2;
        grating_aperture1 = poly2mask(x1, y1, NN(1), NN(1));
        grating_aperture2 = poly2mask(x2, y2, NN(2), NN(2));
    case 6 
        m1 = NN(1)/2;
        r1 = 100;
        x1 = [m1+r1, m1+r1*cosd(60), m1 + r1 * cosd(120), m1+r1*cosd(180), m1+r1*cosd(240), m1+r1*cosd(300)];
        y1 = [m1, m1+r1*sind(60), m1+r1*sind(120), m1+r1*sind(180), m1+r1*sind(240), m1+r1*sind(300)];
        grating_aperture1 = poly2mask(x1, y1, NN(1), NN(1));

        m2 = NN(2)/2;
        r2 = 100;
        x2 = [m2+r2, m2+r2*cosd(60), m2 + r2 * cosd(120), m2+r2*cosd(180), m2+r2*cosd(240), m2+r2*cosd(300)];
        y2 = [m2, m2+r2*sind(60), m2+r2*sind(120), m2+r2*sind(180), m2+r2*sind(240), m2+r2*sind(300)];
        grating_aperture2 = poly2mask(x2, y2, NN(2), NN(2));

end
%% Take the FFT
angularSpectrum1 = fft2(grating_aperture1);
angularSpectrum2 = fft2(grating_aperture2);

%% create Frequency Mask 
fx_single_row1=[0:2^(N(1)-1)-1,-2^(N(1)-1):-1]./X(1,:);
fy_single_row1=[0:2^(N(1)-1)-1,-2^(N(1)-1):-1]./Y(1,:);
fx_single_row2=[0:2^(N(2)-1)-1,-2^(N(2)-1):-1]./X(2,:);
fy_single_row2=[0:2^(N(2)-1)-1,-2^(N(2)-1):-1]./Y(2,:);
[fx1,fy1] = meshgrid(fx_single_row1, fy_single_row1);
[fx2,fy2] = meshgrid(fx_single_row2, fy_single_row2);
propfreq1 = find(abs(fx1).^2+abs(fy1).^2 <1);
propmask1 = zeros(size(fx1));
propmask1(propfreq1) = ones(size(propfreq1));
propfreq2 = find(abs(fx2).^2+abs(fy2).^2 <1);
propmask2 = zeros(size(fx2));
propmask2(propfreq2) = ones(size(propfreq2));

%% Find the Angular Spectrum at Z
angularSpectrumAtZ1 = angularSpectrum1 .* transferFcn(dist(1,:),fx1,fy1) .*propmask1;
angularSpectrumAtZ2 = angularSpectrum2 .* transferFcn(dist(2,:),fx2,fy2) .*propmask2;

%% IFFT to get back to complex fields and intensities 
Uz1 = ifft2(angularSpectrumAtZ1);
Iz1 = Uz1 .* conj(Uz1);

Uz2 = ifft2(angularSpectrumAtZ2);
Iz2 = Uz2 .* conj(Uz2);
%% plot
figure(1);
clf;
imagesc(real(grating_aperture1));
axis image; colorbar;
figure(2);
clf;
imagesc(Iz1)
colorbar; clim([0,0.75]);
figure(3);
clf;
imagesc(Iz2);
colorbar; clim([0,0.75]);
xlim([3700,4500]);
ylim([3700,4500]);















