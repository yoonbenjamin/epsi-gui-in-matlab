%{
This file script creates a GUI on which a read-in proton image and EPSI plot are shown. 
The user has the option to scroll through different proton image slices, contrast 
levels, and EPSI plots.

Author: Shurik Zavriyev
Version: 1.0
_______________________________________________________________________________________

The existing EPSI GUI is maintained with an additional option for the user in which 
they have the option to select a region of interest for a color map to be applied as an 
alternative to superimposing EPSI plots, where the larger the EPSI value the closer to 
yellow and the smaller the closer to blue.

Version: 2.0
Author: Ben Yoon
Date: May 21, 2023
_______________________________________________________________________________________

The color map feature is maintained, but the EPSI GUI is updated such that there now 
exist two panels: one on the top of the user interface (UI) and one on the bottom. The 
image and contrast sliders are place in the lower panel and the EPSI slider, color map 
button, and newly added checkbox are all in the top of the UI. The checkbox starts EPSI 
plotting and the EPSI slider will now have no effect unless the checkbox is checked. 
The user now has the ability to superimpose an EPSI plot ontop of the color map, if 
desired. Also titles are given to sliders and a program title box is displayed.

Version: 2.1
Author: Ben Yoon
Date: May 25, 2023
_______________________________________________________________________________________

Additions to the GUI are made in this update. First, when the color map is displayed, the 
GUI also presents the user with a key to help the user relate the different colors on the 
map to different peak values. Second, the EPSI slider is updated to display which number 
EPSI plot is being shown via a label. 

Version: 2.2
Author: Ben Yoon
Date: Jun 30, 2023
%}

%% Run EPSI GUI

EPSIGUI();

%{
This function sets up the GUI and waits for user input. 

Version: 1.0
Author: Shurik Zavriyev
_______________________________________________________________________________________

Adds the color map feature to the existing EPSI GUI. A button is added that allows the 
user to select a region of interest to which the color map will be applied by calling 
makeMask, instead of superimposing EPSI plots. 

Version: 2.0
Author: Ben Yoon
Date: May 15, 2023
_______________________________________________________________________________________

Updated code to set up the GUI: updates the user interface such that two panels are 
initialized and placed at the top and bottom of the UI. Tools to visualize EPSI are 
located in the top panel (newly added plot EPSI checkbox, color map button, EPSI 
slider). Tools to adjust the proton image are plotted in the bottom panel (image and 
contrast sliders). The new checkbox allows the user to superimpose EPSI plots on top of 
a color map of an ROI. Labels are additionally added to all sliders and a title is 
added to the GUI.

Version 2.1
Author: Ben Yoon
Date: May 22, 2023
_______________________________________________________________________________________

Updated code to create the EPSI data set number label, the text that is to be displayed, 
and to update the action listener for the EPSI slider to give the call to plot_EPSI the 
label and text as parameters to be updated.

Version 2.2
Author: Ben Yoon
Date: Jun 30, 2023
%}
function EPSIGUI()
clear;
close all;

image_parameters.fsems_num = 3;
folder_name = 's_2023041103';
epsi_num = [1];
proton_slice = 9;
image_parameters.xy_shift = [-0.3 -0.4];
xy_shift = [-0.3 -0.4];

image_parameters.scale_by_max_all = true; % Set to true: scale by max spectra val, or 
% to false: scale by max val of each spectrum

image_parameters.moving_avg_window = 1; % Replace with actual value
image_parameters.do_psf_correction = 1;
base_path = [['/Users/benjaminyoon/Desktop/PIGI folder/Projects/Project0 GUI/Mouse ' ...
    'Kidney Data/'] folder_name '/'];
dcm_path = [base_path 'fsems_rat_liver_' num2str(image_parameters.fsems_num.', ...
    '%02d') '.dmc/'];
c13_data_path = [base_path 'epsi_16x12_13c_' num2str(epsi_num.','%02d')];
do_psf_correction = 1;
perform_b1_correction = 0;
study_details.lw_proton = 60;
study_details.centric = 1;
study_details.nimg_to_process = 1;

width_h1 = readprocpar(dcm_path(1:end-5),'lro'); % get image & overlay dimensions
width_h1 = width_h1(2) * 10;
height_h1 = readprocpar(dcm_path(1:end-5),'lpe');
height_h1 = height_h1(2) * 10;
width_c13 = readprocpar(c13_data_path,'lro');
width_c13 = width_c13(2) * 10;
height_c13 = readprocpar(c13_data_path,'lpe');
height_c13 = height_c13(2) * 10;
num_columns_epsi = 16;
num_rows_epsi = 12;

fig = figure(); % Create figure and subplots
fig.Position = [200 100 800 800];

file_list = dir(fullfile(dcm_path, '*.dcm')); % Add message showing the # DICOM files

subfolders = dir(base_path); % Get all subfolders in the parent folder

subfolders = subfolders([subfolders(:).isdir]); % Keep only folders

subfolders = subfolders(3:end); % Remove . and ..

count = 0; % Count the # of subfolders that start with 'epsi'

for i = 1:length(subfolders)
    if and(startsWith(subfolders(i).name, 'epsi'), endsWith(subfolders(i).name, 'fid'))
        count = count + 1;
    end
end

disp(['Number of subfolders starting with "epsi": ' num2str(count)]); % Display result
number_epsis = count;

ha = axes('Parent', fig,'units','normalized', 'position',[0 0 1 1]); % Create figure & 
% axes to plot
uistack(ha,'bottom');
set(ha,'handlevisibility','off', 'visible','off')

epsi_ax = axes('Parent', fig,'units','normalized', 'Position', [((width_h1 - ...
    width_c13) / 2 + xy_shift(1) * width_c13 / num_columns_epsi) / width_h1, 1 - ...
    ((height_h1 - height_c13) / 2 + num_rows_epsi * height_c13 / num_rows_epsi - ...
    xy_shift(2) * height_c13 / num_rows_epsi) / height_h1, width_c13 / width_h1, ...
    height_c13 / height_h1]); % Position subplot in figure
set(epsi_ax,'handlevisibility','off', 'visible','off')

grid_ax = axes('Parent', fig,'units','normalized', 'Position', [((width_h1 - ...
    width_c13) / 2 + xy_shift(1) * width_c13 / num_columns_epsi) / width_h1, ...
    ((height_h1 - height_c13) / 2 + xy_shift(2) * height_c13 / num_rows_epsi) / ...
    height_h1, width_c13 / width_h1, height_c13 / height_h1]); % Pos subplot in fig
set(grid_ax,'handlevisibility','off', 'visible','off')

%{
Updated code to create sliders: updated parent and position to plot in the desired 
location within newly initialized panels for image adjustments and EPSI visualizations 
in the new Version 2.1 GUI (Yoon. 2023. Version 2.1).
%}
panelEPSIOptions = uipanel('Title', 'EPSI Options', 'Position', [0.05 0.845 0.9 0.1]);
panelProtonImgOptions = uipanel('Title', 'Proton Image Options', 'Position', [0.05 ...
    0.045 0.9 0.1]);
slider = uicontrol(panelProtonImgOptions, 'Style', 'slider', 'Units', 'normalized', ...
    'Position', [0.11, 0.6, 0.838, 0.225], 'Value', proton_slice, 'Min', 1, 'Max', ...
    length(file_list), 'SliderStep', [1 / (length(file_list) - 1), ...
    1 / (length(file_list) - 1)], 'FontSize', 16);
contrast_slider = uicontrol(panelProtonImgOptions, 'Style', 'slider', 'Units', ...
    'normalized', 'Position', [0.12, 0.2, 0.83, 0.225], 'Value', .5, 'Min', 0, ...
    'Max', 1, 'SliderStep', [0.05, 0.1], 'FontSize', 16);
EPSI_slider = uicontrol(panelEPSIOptions, 'Style', 'slider', 'Units', 'normalized', ...
    'Position', [0.12, 0.49, 0.83, 0.225], 'Value', epsi_num, 'Min', 1, 'Max', ...
    number_epsis, 'SliderStep', [1 / (number_epsis - 1), 1 / (number_epsis - 1)], ...
    'FontSize', 16);

addlistener(slider, 'Value', 'PostSet', @(~,~)update_plot(ha, file_list, slider, ...
    contrast_slider)); % Listener for the img slider to update plot when moved

addlistener(contrast_slider, 'Value', 'PostSet', @(~,~)update_plot(ha, file_list, ...
    slider, contrast_slider)); % Listener for contrast slider to update plot when moved

%{
Updated code EPSI slider listener: updated the call to update_EPSI to include the 
parameter showPlot added in version 2.1 based on the current value of the newly 
initialized checkbox so that the data in update_EPSI is loaded and shown if checked, 
and only loaded if unchecked. Checking the box also plots the first EPSI plot (Yoon. 
2023. Version 2.1).
*************************************************************************************
Updated code to create a label for the EPSI slider denoting which EPSI data set the program 
is currently looking at. label is placed within the EPSI options panel and is updated 
every time the user interacts with the EPSI slider. update_EPSI is updated to take label 
and labelText as parameters so that they can be updated accordingly (Yoon. 2023. Version
2.2).
%}
checkboxPlotEPSI = uicontrol(panelEPSIOptions, 'Style', 'checkbox', 'String', ...
    'Graph Data Plotting', 'Units', 'normalized', 'Position', [0.05, 0.75, 0.16, ...
    0.225], 'Value', 0, 'Callback', @checkboxCallback);

textEPSIDataSet = 'EPSI Data Set Number:';
labelEPSIDataSet = uicontrol(panelEPSIOptions, 'Style', 'text', 'Units', 'normalized', ...
    'Position', [0.45, 0.75, 0.16, 0.225], 'String', textEPSIDataSet);

addlistener(EPSI_slider, 'Value', 'PostSet', @(~,~)update_EPSI(EPSI_slider, ...
    study_details, epsi_ax, grid_ax, base_path, image_parameters, ...
    get(checkboxPlotEPSI, 'Value'), false, labelEPSIDataSet, textEPSIDataSet));

%{
Adds the color map button to the GUI and a listener to call selectRegion if pressed 
(Yoon. 2023. Version 2.0).
***************************************************************************************
Updated code to create button: updated position to fit in its desired position on the 
top panel in the version 2.1 GUI and listener is made later (Yoon. 2023. Version 2.1).
%}
buttonSelectROI = uicontrol(panelEPSIOptions, 'Style', 'pushbutton', 'String', ...
    'Select Region (Color Map)', 'Units', 'normalized', 'Position', [0.05 0.105 0.9 ...
    0.225], 'Visible', 'on');

%{
Adds titles for the sliders using accordingly positioned text boxes and a large, 
bold, colored (light magenta) program title (Yoon. 2023. Version 2.1).
%}
textTitleSliderEPSI = uicontrol(panelEPSIOptions, 'Style', 'text', 'Units', ...
    'normalized', 'Position', [0.05 0.49 0.06 0.225], 'String', 'Data Set', ...
    'HorizontalAlignment', 'left');
textTitleSliderImg = uicontrol(panelProtonImgOptions, 'Style', 'text', 'Units', ...
    'normalized', 'Position', [0.05 0.6 0.05 0.225], 'String', 'Image', ...
    'HorizontalAlignment', 'left');
textTitleSliderContrast = uicontrol(panelProtonImgOptions, 'Style', 'text', ...
    'Units', 'normalized', 'Position', [0.05 0.2 0.06 0.225], 'String', ...
    'Contrast', 'HorizontalAlignment', 'left');
textTitle = uicontrol('Style', 'text', 'Units', 'normalized', 'Position', [0.25 ...
    0.95 0.5 0.04], 'String', 'Echo-Planar Spectroscopic Imaging', ...
    'HorizontalAlignment', 'center');
set(textTitle, 'BackgroundColor', [0.54, 0.28, 0.54]);
set(textTitle, 'ForegroundColor', 'white');
set(textTitle, 'FontSize', 20);
set(textTitle, 'FontWeight', 'bold');
%{
Nested function, called when the user clicks the checkbox, plots EPSI for the 1st time
%}
    function checkboxCallback (~, ~)
        update_EPSI(EPSI_slider, study_details, epsi_ax, grid_ax, base_path, ...
            image_parameters, true, true, labelEPSIDataSet, textEPSIDataSet);
    end

update_plot(ha, file_list, slider, contrast_slider); % update_plot to plot 1st image

%{
Updated code to call update_EPSI: updated so that update_EPSI returns spectral data 
and stores it in fieldSpectralData, as it is useful outside of update_EPSI as well 
(Yoon. 2023. Version 2.0). 
***************************************************************************************
Updated again to pass false as the value for showPlot so that the EPSI is not plotted 
until the checkbox is checked, but data can still be loaded (Yoon. 2023. Version 2.1).
***************************************************************************************
Updated again to take the EPSI data set number label and the text for the label as parameters 
to updated accordingly (Yoon. 2023. Version 2.2).
%}
fieldSpectralData = update_EPSI(EPSI_slider, study_details, epsi_ax, grid_ax, ...
    base_path, image_parameters, false, false, labelEPSIDataSet, textEPSIDataSet);

%{
Adds the color map button listener last so when the button is clicked it can pass the 
spectral data derived from update_EPSI into makeMask (Yoon. 2023. Version 2.1).
%}
addlistener(buttonSelectROI, 'Value', 'PostSet', @(~,~)makeMask(slider, file_list, ...
    contrast_slider, ha, fieldSpectralData));
end

%% Plot a proton image

%{
This function plots a desired proton image with the desired contrast.

Version 1.0
Author: Shurik Zavriyev
%}
function update_plot(ax, file_list, slider, contrast_slider)
idx = round(get(slider, 'Value')); % Get slider val to find the desired img

img = dicomread(fullfile(file_list(idx).folder, file_list(idx).name)); % Get img from 
% DICOM

contrast_value = get(contrast_slider, 'Value'); % Get desired contrast

low = (1 - contrast_value) / 2; % Adjust contrast of img
high = 1 - low;
contrast_image = imadjust(img, [low, high], []);

imagesc(ax, contrast_image); % Plot the image
colormap(ax, gray);
set(ax,'handlevisibility','off', 'visible','off')

title(ax, sprintf('Image %d of %d', idx, length(file_list)), 'FontSize', 16); % Update 
% title to show curr img index & # images
end

%% Handle EPSI plots and spectral data

%{
This function plots a desired EPSI plot on top of a proton image.

Version 1.0
Author: Shurik Zavriyev
_______________________________________________________________________________________

Updated code to plot EPSI plots: updated to return the initialized spectral data so 
that it can be used in the rest of the script. Also updated to take a boolean showPlot 
as a parameter such that the function will load and plot the EPSI plot if showPlot is 
true, and load but not plot if false.

Version 2.0
Author: Ben Yoon
Date: May 17, 2023
_______________________________________________________________________________________

Updated code to update the EPSI data set number label. Additional parameters labelEPSIData 
and textEPSIData are included which should be the UI label and the text that is displayed 
in the label. This function now will use the current value of the EPSI slider to set the 
text of labelEPSIData appropriately. 

Version 2.2
Author: Ben Yoon
Date: Jun 30, 2023
%}
function dataSpectral = update_EPSI(EPSI_slider, study_details, epsi_ax, grid_ax, ...
    base_path, image_parameters, showPlot, isFirstRun, labelEPSIData, textEPSIData)
epsi_num = get(EPSI_slider, 'Value');
dcm_path = [base_path 'fsems_rat_liver_' num2str(image_parameters.fsems_num.', ...
    '%02d') '.dmc/'];
c13_data_path = [base_path 'epsi_16x12_13c_' num2str(epsi_num.','%02d')];
width_h1 = readprocpar(dcm_path(1:end-5),'lro');
width_h1 = width_h1(2) * 10;
height_h1 = readprocpar(dcm_path(1:end-5),'lpe');
height_h1 = height_h1(2) * 10;
width_c13 = readprocpar(c13_data_path,'lro');
width_c13 = width_c13(2) * 10;
height_c13 = readprocpar(c13_data_path,'lpe');
height_c13 = height_c13(2) * 10;
num_columns_epsi = 16;
num_rows_epsi = 12;

spectral_data = new_epsi_recon(study_details, c13_data_path, 1, ...
    study_details.lw_proton / 4); % Initialize spectral data and scaling option
spectral_data = flip(flip(spectral_data, 1), 2);
%if image_parameters.do_psf_correction
%    spectral_data = psf_correction(spectral_data);
%end

if image_parameters.scale_by_max_all % Scale spectral data as necessary
    max_val = max(spectral_data(:)); % Find max spectra val
    spectral_data = spectral_data ./ max_val; % Scale all spectra by total max val
else
    max_vals = max(spectral_data, [], 3); % Find max for each spectrum
    spectral_data = spectral_data ./ max_vals; % Scale each spectrum by its max val
end

total_data = [];

for i = 1:num_rows_epsi
    row_data = [];
    for j = 1:num_columns_epsi
        if max(spectral_data(i, j, :)) < 0.20 %&& showPlot
            spectral_data(i, j, :) = nan;
        end
        row_data = [squeeze(row_data);
            circshift(squeeze(spectral_data(i, j, :)), 0)];
    end
    total_data = [squeeze(total_data);
        squeeze(row_data + num_rows_epsi - i)];
end

x_vals_for_plotting = repmat(1:length(spectral_data(1, 1, :)) * num_columns_epsi, ...
    1, num_rows_epsi);

for nanpoints = 1:num_rows_epsi-1
    total_data(nanpoints * length(spectral_data(1, 1, :)) * num_columns_epsi) = nan;
end

total_data = movmean(total_data, image_parameters.moving_avg_window, 'Endpoints', ...
    'fill'); % Apply moving average

disp("x_vals_for_plotting")
disp(size(x_vals_for_plotting))

if showPlot
    if epsi_num == 1
        plot(epsi_ax, x_vals_for_plotting, squeeze(total_data), 'Color', 'y', ...
            'LineWidth', 2);
    else        
        plot(epsi_ax, x_vals_for_plotting, squeeze(total_data), 'Color', 'magenta', ...
            'LineWidth', 2);
    end
    ylim([0 num_rows_epsi]);
    xlim([0 num_columns_epsi * length(spectral_data(1, 1, :))]);
    axis off;
    if isFirstRun
        yline(grid_ax, 0:num_rows_epsi, '--w', 'LineWidth', 1, 'Alpha', 0.5);
        xline(grid_ax, 0:num_columns_epsi, '--w', 'LineWidth', 1, 'Alpha', 0.5);
        ylim(grid_ax, [0 num_rows_epsi]);
        xlim(grid_ax, [0 num_columns_epsi]);
        axis off;
    end
end

%{
Updated code: gets the current EPSI value and sets the value of the EPSI data set number 
label accordingly (Yoon. 2023. Version 2.2).
%}
valueEPSISlider = get(EPSI_slider, 'Value');
newTextEPSI = [textEPSIData ' ' num2str(valueEPSISlider)];
set(labelEPSIData, 'String', newTextEPSI);

%{
Updated code: returns the loaded spectral data (Yoon. 2023. Version 2.0).
%}
dataSpectral = spectral_data;
end

%% Select region of interest (ROI), build and show the appropriate color map

%{
This function, called when the user clicks the color map button, replots the proton 
image on top of the EPSI plots to allow the user to select an ROI using drawfreehand, 
which when double clicked is used to create a binary mask matrix (1, ROI, 0, otherwise) 
to which the color map will be applied through calling buildColorMap.

Version 2.0
Author: Ben Yoon
Date: May 18, 2023
%}
function makeMask(image, list, contrast, ax, spectral)
idx = round(get(image, 'Value'));
dicom = dicomread(fullfile(list(idx).folder, list(idx).name));
contrastDesired = get(contrast, 'Value');
low = (1 - contrastDesired) / 2;
hi = 1 - low;
dicomAdjusted = imadjust(dicom, [low, hi], []);
imagesc(ax, dicomAdjusted);
axis(ax, 'image');
colormap(ax, gray);
uistack(gca, 'top');

roi = drawfreehand(ax);
wait(roi);

binaryMask = createMask(roi);

buildColorMap(binaryMask, spectral, dicomAdjusted, ax);
end

%{
This function creates a color map over an ROI, given as a binary mask, over an image. 
The color map is loaded with RGB values computed through calling getColor passing each 
white pixel's coordinates and the EPSI plot. The EPSI plot is passed as a 2 dimensional 
array with coordinate and max values for each plot, built through calling 
build2DArrayEpSIPlot. This loaded color map is then made transparent and superimposed 
over the proton image.

Version 2.0
Author: Ben Yoon
Date: May 19, 2023
_____________________________________________________________________________________

Updated code to display a key for the color map to allow the user to better visualize 
what the colors represent (blue to yellow is interpolated from 0 to 1). Key was designed 
by me and is appropriately sized and positioned in the code below. 

Version 2.2
Author: Ben Yoon
Date: Jun 30, 2023
%}
function buildColorMap(mask, sData, dicom, ax)
colorMap = zeros(size(mask, 1), size(mask, 2), 3); % Initialize color map
[row, col] = find(mask == 1);
arrayEPSIPlot = build2DArrayEPSIPlot(sData);

for i = 1:numel(row)
    x = col(i);
    y = row(i);

    color = getColor(x, y, arrayEPSIPlot);

    redValue = color(1);
    greenValue = color(2);
    blueValue = color(3);

    colorMap(y, x, 1) = redValue;
    colorMap(y, x, 2) = greenValue;
    colorMap(y, x, 3) = blueValue;
end

imagesc(ax, dicom);
axis(ax, 'image');
colormap(ax, gray);
hold on
im = image(colorMap);
im.AlphaData = 0.375;
hold off

%{
Updated code: udpated to read in the mapKey image file (made by me) and to
resize, position, and display it in the desired place in the GUI (Yoon.
2023. Version 2.2).
%}
mapKey = imread(['/Users/benjaminyoon/Desktop/PIGI folder/Projects/Project0 ' ...
    'GUI/Code/Kidney Code/YoonMapKey.png']);
desiredSize = [95 15];
resizedMapKey = imresize(mapKey, desiredSize);
hold on;
shiftX = 225;
shiftY = 80;
imshow(resizedMapKey, 'XData', [shiftX shiftX + size(resizedMapKey, 2) - 1], ...
    'YData', [shiftY shiftY + size(resizedMapKey, 1) - 1]);
hold off;
end

%{
This function initializes, loads, and returns a 192 x 5 matrix that represents an EPSI 
plot. The 192 rows are the EPSI plots and the columns store, in order, the following 
data about each plot: minimum x coordinate, maximum y coordinate, maximum x coordinate, 
minimum y coordinate, maximum EPSI value. Data is loaded based on a ratio of pixels to 
x and y coordinates to find the coordinates of each plot and max is used to find the 
maximum EPSI value. Ratio is found given subplots are 14 x 12 pixels, proton img is 256 
x 256, epsi plot has top left corner coordinate (0.055, 0.8) & coordinates are 
normalized so 256 corresponds to 1.

Version 2.0
Author: Ben Yoon
Date: May 20, 2023
%}
function EPSIPlot = build2DArrayEPSIPlot(data)
plot = zeros(192, 5);

for i = 1:192
    [row, col] = ind2sub([12, 16], i); % Get row & col vals based on order 
    % (arbitrarily: left to right, top to bottom)
    plot(i, 1) = 0.055 + ((14.0 / 256.0) * (col - 1.0));
    plot(i, 2) = 0.2375 + ((12.0 / 256.0) * (row));
    plot(i, 3) = 0.055 + ((14.0 / 256.0) * col);
    plot(i, 4) = 0.2375 + ((12.0 / 256.0) * (row - 1.0));
    plot(i, 5) = max(data(row, col, :));
end

EPSIPlot = plot;
end

%{
This function takes the EPSI plot as a 2 dimensional array and returns the 
corresponding RGB value given as an array where indices 1, 2, 3 are R, G, B. The array 
is derived from given pixel coordinates and max EPSI value. The pixel coordinates are 
checked against each EPSI plot's coordinates until one that overlaps is found. The 
maximum EPSI value in this subplot [0, 1] is interpolated between yellow and blue, 
where yellow is 1 and blue is 0.

Version 2.0
Author: Ben Yoon
Date: May 21, 2023
%}
function arrayRGB = getColor(xCoord, yCoord, data)
maxEPSI = -1;

for i = 1:192
    if xCoord / 256 >= data(i, 1) && yCoord / 256 <= data(i, 2) && ...
        xCoord / 256 <= data(i, 3) && yCoord / 256 >= data(i, 4)
        maxEPSI = data(i, 5);
        break;
    end
end

blue = [0 0 1];
lightBlue = [0 0.5 1];
lighterBlue = [0 1 1];
blueGreen = [0 1 0.5];
green = [0 1 0];
yellowGreen = [0.5 1 0];
yellow = [1 1 0];
colors = [blue; lightBlue; lighterBlue; blueGreen; green; yellowGreen; yellow];
colorValues = [0 0.17 0.33 0.5 0.67 0.83 1];
valRGB = interp1(colorValues, colors, maxEPSI); % Interpolate maxEPSI val b/w blue & 
% yellow
arrayRGB = valRGB;
end