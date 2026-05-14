/**
 * MACRO TO PRE-PROCESS AND THRESHOLD IMMUNOFLUORESCENT CELLS
 * Tim Monko - July 20, 2021 - timmonko@gmail.com - https://github.com/timmonko
 * Helpful resources as a starting point for  https://imagej.net/imaging/principles
 *
 * Scripting parameters (i.e #@ __ designed to create a UI to run process as a batch, see comments and UI notes for further detail)
 * While this is not the ideal way for optional scripting parameters to be used (e.g., the UI scaling from 1 to infinite channels interested to use
 * To add more channels, you would need to copy and paste the channel-specific code, and then add another function call
 *
 * This macro also maintains some legacy parameters (such as Median or Guassian Weighted Median),
 * for compatability with previous iterations of thresholding analysis, but also only maintains reasonable processing (i.e., no longer watershed/watershed irregular features)
 */

#@ File(label = "Input directory", style = "directory") input
#@ File(label = "Output directory", style = "directory") output
#@ String(label = "File suffix", value = ".tif") suffix
#@ Integer(label = "# of Channels", value = 3, min = 1, max = 3, style = "slider") num_channels

#@ Boolean(label = "Image has selected ROI") ROI_boolean
#@ String(label = "Crop to ROI", choices=("No", "Yes, polygon", "Yes, rectangle"), style = "radioButtonHorizontal") ROI_crop

#@ String(label = "Filter Type", choices=("Gaussian Weighted Median", "Median..."), style = "radioButtonHorizontal") filter_type
#@ String(label = "Brightness Normalization", choices = ("Auto", "Manual"), style = "radioButtonHorizontal") normalization_type
#@ String(label = "Threshold Type (Otsu)", choices=("Auto", "Normalized Perc of Max"), style = "radioButtonHorizontal") threshold_type

#@ String(label = "1-Channel Identifier", value = "C1") channel1
#@ String(label = "1-Rename to?", value = "C1") rename1
#@ Integer(label = "1-Background Subtraction Radius", value = 15) bs_radius1
#@ Double(label = "1-Filter Radius", value = 1.5) filter_radius1
#@ Integer(label = "1-Prominence", value = 10) prominence1
#@ Double(label = "1-ThresholdPercent", value = 0.15) threshold1


#@ String(label = "2-Channel Identifier", value = "C2") channel2
#@ String(label = "2-Rename to?", value = "C2") rename2
#@ Integer(label = "2-Background Subtraction Radius", value = 15) bs_radius2
#@ Double(label = "2-Filter Radius", value = 1.5) filter_radius2
#@ Integer(label = "2-Prominence", value = 10) prominence2
#@ Double(label = "2-ThresholdPercent", value = 0.15) threshold2


#@ String(label = "3-Channel Identifier", value = "C3") channel3
#@ String(label = "3-Rename to?", value = "C3") rename3
#@ Integer(label = "3-Background Subtraction Radius", value = 15) bs_radius3
#@ Double(label = "3-Filter Radius", value = 1.5) filter_radius3
#@ Integer(label = "3-Prominence", value = 10) prominence3
#@ Double(label = "3-ThresholdPercent", value = 0.15) threshold3


#@ Boolean(label = "Remove Small Particles?") particles_boolean
#@ Integer(label = "Particle Size") particle_size

// Prevents image windows from opening while the script is running -- unless altered with false, "show", or "hide"
setBatchMode(true);

// For counting the iterations within the function, define the total count of folders, as well as the number of loops
// 	loop starts at 1 instead of zero because it is not being used in a for loop and therefore the i-th index of the looping number starts at 0
file_list = getFileList(input)
count = lengthOf(file_list);
loop = 1;

// This function calls the functions that follow, using the input folder defined by the folder input parameter
processFolder(input);

print("Success");

// Scan folders/subfolders to find files with the correct suffix
function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if (File.isDirectory(input + File.separator + list[i]))
			processFolder(input + File.separator + list[i]);
		if (endsWith(list[i], suffix))
			processFile(input, output, list[i]);
	}
}

function processFile(input, output, file) {

	// Open image using Bio-Formats and then close previous images -- does not seem to open ROIs unless the file is save with Bio-Formats and its inherent ROI management. To preserve backwards compatability, I've chosen to continue with default ImageJ opening, which will show all the ROIs, allowing the Active ROI to be added
	//run("Bio-Formats Importer", "open=[" + input + File.separator + file +"] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
	open(input + File.separator + file);
	close("\\Others");
	title = getTitle();

	// Reset the ROIs in the manager, and then add the active ROI to the manager, to be added back onto the finished image later
	roiManager("reset");
	if (ROI_boolean == 1) {
		roiManager("add");
		run("Select All");
	}
	if (num_channels > 1) {run("Split Channels"); }

	// Function to run code on an individual channel
	processChannel(channel1, bs_radius1, filter_radius1, prominence1, threshold1, rename1, particle_size);
	if (num_channels > 1) {	processChannel(channel2, bs_radius2, filter_radius2, prominence2, threshold2, rename2, particle_size); }
	if (num_channels > 2) {	processChannel(channel3, bs_radius3, filter_radius3, prominence3, threshold3, rename3, particle_size); }

	// Close images that are not named with the title (keep the others), this is because all processed images will have been renamed. Keeping the original 16-bit images with 8-bit binary does not work. Of course, the original file remains.
	close("*" + title, "keep");
	if (num_channels > 1) { run("Images to Stack"); }
	run("8-bit");

	// Rename image (if altered, such as with blinding), then save it with the binned prefix to clarify origin
	rename(title);

	// A convenience function to crop, and then optionally clear outside a selected polygon ROI selection
	if (startsWith(ROI_crop, "Yes")) {
		roiManager("Select", 0);
		run("Crop");
		if(ROI_crop == "Yes, polygon") { run("Clear Outside", "stack");	}
	}

  // Any ROIs present on the image will be sent to the overlay. Then, the 0th indexed ROI will be selected for saving, so next time it is open it will be on the roiManager
	run("To ROI Manager");
	roiManager("Select", 0);
	run("Add Selection...");

  // Save, and the send to log window the iterations of the loop on a file, will only add loop++ if in the processFile of the processFolder function, then printing the image title nicely, as a sort of progress bar
	saveAs(output + File.separator + "thresh_" + title);
	print(loop++ + " out of " + count + " : " + title);
}

function processChannel(channel, bs_radius, filter_radius, prominence_size, threshold_perc, rename_, particle_size) {
	selectWindow(channel + "-" + title);
	run("Duplicate...", " ");
	rename(channel);
	run("Grays");

	// Uses Bio-voxxels convoluted background subtraction https://imagej.net/plugins/biovoxxel-toolbox
	run("Convoluted Background Subtraction", "convolution = Median radius=bs_radius");

	// Use either ImageJ's median filter or Bio-voxxels gaussian-weight median filter. The latter performs better to keep brightness of cells, while effectively preserving the same edges
	run(filter_type, "radius=filter_radius");

	// Detect maxima to create a voronoi segmented particles map -- later add this to the thresholded image to create a prominence based watershed algorithm
	run("Find Maxima...", "prominence=prominence_size output=[Segmented Particles]");
	selectWindow(channel);

	// If desired, manually select areas of signal, while ignoring intense background (such as debris, or off-target signal like autofluorescent blood vessels)
	if (ROI_boolean == 1) { roiManager("Select", 0); }
	if (normalization_type == "Manual") {
		setBatchMode("show");
		waitForUser("Manual Intensity ROI", "Select area with representative highest intensity, avoiding background false positives");
		setBatchMode("hide");
	}

	// Normalize the image to the brightest pixels -- this normalization will not change the algoirthm for Auto Threshold, so results will be the same
	List.setMeasurements;
	max_intensity = List.getValue("Max");
	run("Select All");
	run("32-bit");
	run("Divide...", "value=max_intensity");


	// Threshold settings, will do auto threshold with Otsu's algorithm (ImageJ's "Default" is very similar.
	// Using the "Perc of Max" gives a representative "Percentage" cutoff for pixel brightness after normalization by division of maximum
	setAutoThreshold("Otsu dark");
	if (threshold_type == "Normalized Perc of Max") { setThreshold(threshold_perc, 100); }
	setOption("BlackBackground", true);
	run("Convert to Mask");

	// Combine the Thresholded and Proxima Segmented Particles image to create the Thresholded+Watershed final image
	imageCalculator("AND create", channel, channel + " Segmented");
	close(channel + " Segmented");
	close(channel);
	selectWindow("Result of " + channel);
	rename(rename_);

	// Remove small particles: particles 4 is just up-down-left-right adjacency, not orthogonal, so this is best to be used since we threshold and would otherwise have many 8-connected particles
	if (particles_boolean == 1) { run("Particles4 ", "white show=Particles filter minimum=particle_size maximum=9999999 redirect=None"); }

	// Add the original active ROI back to the image
	if (ROI_boolean == 1) {
	roiManager("Select", 0);
	run("Add Selection...");
	}
}
