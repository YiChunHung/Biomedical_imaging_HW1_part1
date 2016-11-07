% Introduction to Biomedical Imaging,   Fall 2016
% EE 441000 Yi-Chun Hung 103061145 HW1 Part1 11/05/2016
%                                      

% ----------- Quiz 1 ------------
% Note that you can use "SI" unit instead of the units I used below
close all;
clear all;
time_offset = 6.48 * 10^(-6); % in usec
fs = 50 * 10^6;    % sampling rate, in MHz
aper_size = 6;  % aperture size, in mm
%focal_pt = ;   % focal point, in mm
dx = 0.050;    % distance between two successive scanning positions, in mm
soundv = 1.54 * 10^6 ;  % speed of sound, in mm/us

% (a)
DR = 100;    % dynamic range, in dB
load('points_rf_data');
rf_data = points_rf_data; 
clear points_rf_data;
[m,n] = size(rf_data);
z_axis = ((0:m-1)/fs+time_offset)*soundv/2;   % z axis
x_axis = (0:n-1)*dx;    % x axis

envolope = hilbert(rf_data);
envolope = abs(envolope);   % envelope detection
envolope_dB = 20*log10(envolope/max(max(envolope)));    % log conversion with respect to the maximum value 


figure
image(x_axis,z_axis,envolope_dB+DR); % or imagesc()?
colormap(gray(DR));
colorbar;
title('Point targets');
xlabel('Lateral position (mm)');
ylabel('Depth (mm)');
axis image;


% (b)
% observe the image in (a) and find the focal length of the transducer
% one of the point targets is located at the focal point

% In my observation, the focal point should be at third point, so
%	focal length = 12.32 mm


% (c)
% check the spatial resolution based on the fundamental definition of the spatial resolution and your eye examination
% for example, in 1D (can be done in 2D) 
% distance = 10; % in terms of sample points. => you have to convert true distance to sample points
% scatterer_x = [ zeros(1,100) 1 zeros(1, distance) 1 zeros(1,100)]
% PSF_x = max(PSF); % lateral beam profile. projection along the depth
% image_x = conv(scatterer_x, PSF_x); % vary the distance, and then try the eye examination

lat_res = 1; %lateral resolution in mm
axial_res = 0.4; %Axial resolution in mm
x_dis = floor(lat_res/dx); %calculate the distance in pixel
z_dis = floor(axial_res * 2/soundv * fs);

%create the new axis
new_x_axis = [(0:n+x_dis-1)*dx];
new_z_axis = [((0:m+z_dis-1)/fs+time_offset)*soundv/2];

%place two image in one image
scatter_xlat1 = ([envolope_dB+DR zeros(m,x_dis)] > 0) .* ([envolope_dB+DR zeros(m,x_dis)]);
scatter_xlat2 = ([zeros(m,x_dis) envolope_dB+DR] > 0) .* ([zeros(m,x_dis) envolope_dB+DR]);
scatter_xlat = scatter_xlat1 + scatter_xlat2 ;
scatter_zlat1 = ([envolope_dB+DR ; zeros(z_dis,n)] > 0) .* ([envolope_dB+DR ; zeros(z_dis,n)]);
scatter_zlat2 = ([zeros(z_dis,n) ; envolope_dB+DR] > 0) .* ([zeros(z_dis,n) ; envolope_dB+DR]);
scatter_zlat = scatter_zlat1 + scatter_zlat2;

scatter_xlat_try = ([envolope zeros(m,x_dis)] + [zeros(m,x_dis) envolope]);
scatter_xlat_try = 20*log10(scatter_xlat_try/max(max(scatter_xlat_try)));

%draw the image
figure
image(new_x_axis, z_axis, scatter_xlat);
colormap(gray(DR));
colorbar;
title('x axis resolution');
xlabel('Lateral position (mm)');
ylabel('Depth (mm)');
axis image;


figure
image(x_axis, new_z_axis, scatter_zlat);
colormap(gray(DR));
colorbar;
title('z axis resolution');
xlabel('Lateral position (mm)');
ylabel('Depth (mm)');
axis image;

%(d)
% check it point target by point target
% e.g., for point target 1, PSF_point1 = envelope_dB(1:190,:);
% then, find the FWHM (-6dB) and - 20dB for the axial and lateral beam plot,
% respectively.
dB_6 = -6;
dB_20 = -20;

%Projection on the laterial axis
PSF_point1 = max(envolope_dB(1:200,:));
PSF_point2 = max(envolope_dB(201:400,:));
PSF_point3 = max(envolope_dB(401:600,:));
PSF_point4 = max(envolope_dB(601:800,:));
PSF_point5 = max(envolope_dB(801:end,:));




%find the bandwidth of -6dB and -20dB
PSF_point1_lat_res6 = find(PSF_point1 >  (max(PSF_point1) + dB_6)) ;
PSF_point2_lat_res6 = find(PSF_point2 >  (max(PSF_point2) + dB_6)) ;
PSF_point3_lat_res6 = find(PSF_point3 >  (max(PSF_point3) + dB_6)) ;
PSF_point4_lat_res6 = find(PSF_point4 >  (max(PSF_point4) + dB_6)) ;
PSF_point5_lat_res6 = find(PSF_point5 >  (max(PSF_point5) + dB_6)) ;
p1_FWHM6 = (PSF_point1_lat_res6(end) - PSF_point1_lat_res6(1)) * dx;
p2_FWHM6 = (PSF_point2_lat_res6(end) - PSF_point2_lat_res6(1)) * dx;
p3_FWHM6 = (PSF_point3_lat_res6(end) - PSF_point3_lat_res6(1)) * dx;
p4_FWHM6 = (PSF_point4_lat_res6(end) - PSF_point4_lat_res6(1)) * dx;
p5_FWHM6 = (PSF_point5_lat_res6(end) - PSF_point5_lat_res6(1)) * dx;
PSF_point1_lat_res20 = find(PSF_point1 >  (max(PSF_point1) + dB_20)) ;
PSF_point2_lat_res20 = find(PSF_point2 >  (max(PSF_point2) + dB_20)) ;
PSF_point3_lat_res20 = find(PSF_point3 >  (max(PSF_point3) + dB_20)) ;
PSF_point4_lat_res20 = find(PSF_point4 >  (max(PSF_point4) + dB_20)) ;
PSF_point5_lat_res20 = find(PSF_point5 >  (max(PSF_point5) + dB_20)) ;
p1_FWHM20 = (PSF_point1_lat_res20(end) - PSF_point1_lat_res20(1)) * dx;
p2_FWHM20 = (PSF_point2_lat_res20(end) - PSF_point2_lat_res20(1)) * dx;
p3_FWHM20 = (PSF_point3_lat_res20(end) - PSF_point3_lat_res20(1)) * dx;
p4_FWHM20 = (PSF_point4_lat_res20(end) - PSF_point4_lat_res20(1)) * dx;
p5_FWHM20 = (PSF_point5_lat_res20(end) - PSF_point5_lat_res20(1)) * dx;



%Project on the axial axis
PSF_point1_z = max(envolope_dB(1:200,:)');
PSF_point2_z = max(envolope_dB(201:400,:)');
PSF_point3_z = max(envolope_dB(401:600,:)');
PSF_point4_z = max(envolope_dB(601:800,:)');
PSF_point5_z = max(envolope_dB(801:end,:)');

%find the bandwidth of -6dB and -20dB
PSF_point1_z_res6 = find(PSF_point1_z >  (max(PSF_point1_z) + dB_6)) ;
PSF_point2_z_res6 = find(PSF_point2_z >  (max(PSF_point2_z) + dB_6)) ;
PSF_point3_z_res6 = find(PSF_point3_z >  (max(PSF_point3_z) + dB_6)) ;
PSF_point4_z_res6 = find(PSF_point4_z >  (max(PSF_point4_z) + dB_6)) ;
PSF_point5_z_res6 = find(PSF_point5_z >  (max(PSF_point5_z) + dB_6)) ;
p1_FWHM6_z = (PSF_point1_z_res6(end) - PSF_point1_z_res6(1)) /fs*soundv/2;
p2_FWHM6_z = (PSF_point2_z_res6(end) - PSF_point2_z_res6(1)) /fs*soundv/2;
p3_FWHM6_z = (PSF_point3_z_res6(end) - PSF_point3_z_res6(1)) /fs*soundv/2;
p4_FWHM6_z = (PSF_point4_z_res6(end) - PSF_point4_z_res6(1)) /fs*soundv/2;
p5_FWHM6_z = (PSF_point5_z_res6(end) - PSF_point5_z_res6(1)) /fs*soundv/2;
PSF_point1_z_res20 = find(PSF_point1_z >  (max(PSF_point1_z) + dB_20)) ;
PSF_point2_z_res20 = find(PSF_point2_z >  (max(PSF_point2_z) + dB_20)) ;
PSF_point3_z_res20 = find(PSF_point3_z >  (max(PSF_point3_z) + dB_20)) ;
PSF_point4_z_res20 = find(PSF_point4_z >  (max(PSF_point4_z) + dB_20)) ;
PSF_point5_z_res20 = find(PSF_point5_z >  (max(PSF_point5_z) + dB_20)) ;
p1_FWHM20_z = (PSF_point1_z_res20(end) - PSF_point1_z_res20(1)) /fs*soundv/2;
p2_FWHM20_z = (PSF_point2_z_res20(end) - PSF_point2_z_res20(1)) /fs*soundv/2;
p3_FWHM20_z = (PSF_point3_z_res20(end) - PSF_point3_z_res20(1)) /fs*soundv/2;
p4_FWHM20_z = (PSF_point4_z_res20(end) - PSF_point4_z_res20(1)) /fs*soundv/2;
p5_FWHM20_z = (PSF_point5_z_res20(end) - PSF_point5_z_res20(1)) /fs*soundv/2;




%{
figure
plot(PSF_point1,'LineStyle' , '-',...
					'Color' , 'r');

hold on
plot(PSF_point2,'LineStyle' , '--',...
					'Color' , 'k');

plot(PSF_point3,'LineStyle' , ':',...
					'Color' , 'k');

plot(PSF_point4,'LineStyle' , '-.',...
					'Color' , 'k');

plot(PSF_point5,'LineStyle' , '-',...
					'Color' , 'b');

hold off
%}

% find the FWHM?
% hint: quick solution by Matlab built in function "find()"
%       but "not accurate enough"
%	for better accuracy, you can try to interpolate the provided data (e.g., by "interp()")
%       you may also try your own codes
% for example,
% idx = find(PSF >= 0.5);
% FWHM_PSF = (idx(end) - idx(1))*?;  % in um


% (e)
% theoretic lateral resolution "at the focal point"
%fc = ?; % center freq. of the transducer, in MHz, try fft
f_axis = linspace(-fs/2,fs/2,m); %create the f_axis in FFT
lambda_all = zeros(1,n);

%take mean of every line's labda
for i = 1:n
	rf_fft_col = abs(fftshift(fft(rf_data(:,i)-mean(rf_data(:,i))))); 
	lambda_temp = find(rf_fft_col == max(rf_fft_col));	
	lambda_all(i) = lambda_temp(2);
end

lambda = soundv/(mean(lambda_all)/m*fs - fs/2) ; % wavelength of the center frequency
[depth,~] = find(envolope_dB == max(max(envolope_dB))); %find the depth of focal point
depth = z_axis(depth);
f_number = depth/aper_size;
lateral_theoretic = f_number * lambda; %calculate the theoretic lateral resolution

% (f)

% (g)

% --------------- Quiz 2 -------------------


%Bonus

%figure
%plot(f_axis,abs(fftshift(fft(rf_data(:,60)-mean(rf_data(:,60))))))
