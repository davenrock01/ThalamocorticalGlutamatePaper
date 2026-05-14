/* Macro to bin images with either a fixed rectangular bin and/or a polygon selection using ImageJ's ROI tools
 * Tim Monko - July 20, 2021 - timmonko@gmail.com - github.com/timmonko
 * Scripting parameters (i.e #@ __ designed to create a UI to run process as a batch, see comments and UI notes for further detail)
 */

#@ File(label = "Input directory", style = "directory") input
#@ File(label = "Output directory", style = "directory") output 
#@ String(label = "File suffix", value = ".tif") suffix 

#@ Float (label = "Scale px-to-um", value = 1.575) scale
#@ String (value = "1.575 = 10X ; 0.62 = 4X", visibility = "MESSAGE") scale_hint

#@ Boolean (label = "Pre-sized bin") bin_boolean
#@ String (label = "Height", value = 400) bin_height
#@ String (label = "Width", value = 850) bin_width

#@ Boolean (label = "Freehand bin") freehand_boolean

// Create Array to coerce multi-channel images in to a colorblind friendly format. 
//	More colors could be added, but would only cause confusion with readability of the images
LUT_list = newArray("Cyan", "Magenta", "Yellow", "Grays");

// For counting the iterations within the function, define the total count of folders, as well as the number of loops
// 	loop starts at 1 instead of zero because it is not being used in a for loop and therefore the i-th index of the looping number starts at 0

file_list = getFileList(input)
count = lengthOf(file_list); 
loop = 1;

// This function calls the functions that follow, using the input folder defined by the folder input parameter
processFolder(input);

/* Scan folders/subfolders to find files with the correct suffix 
 * Need to allow a recursive function which allows opening of images nested in subfolders
 * which is the only way to call bio-formats options to be used instead of FIJI defaulting to preferred paramaters.
 * This code is sourced from Templates › ImageJ 1.x › Batch › Process Folder
 * The code will generate the list[i] name, which then gets added to the input path for opening -- eventually opening a filename found with the suffix. 
 * If recursively looking into subdirectories, then list[i] subfolder recursively becomes list[i] + subfolder name
 * Thus, the recursive list[i] can be a subfolder which is looked into for the suffix and passed to processFile function
 */ 
 
function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i]))
			processFolder(input + File.separator + list[i]);
		if(endsWith(list[i], suffix)) 
			processFile(input, output, list[i]); 
	}
}

/* Process each file to add ROIs as defined in header intro
 * A series of code could be added which opens all the folders in the file and then combine them
 * See comments within function for more details
 */

function processFile(input, output, file) {
	setBatchMode(true); // prevents image windows from opening while the script is running
	// open image using Bio-Formats
	run("Bio-Formats Importer", "open=[" + input + File.separator + file +"] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
	title = getTitle();
	// Dimensions used for scaling the bin, and placing the bin in a reasonable centered position on the image
	getDimensions(w,h,channels,slices,frames);
	run("Set Scale...", "distance=scale known=1 pixel=1 unit=µm global"); 
	run("Make Composite", "display=Composite");

	// Use the LUT array defined about to iterate over channels
	for (n = 0; n < channels; n++) {
			Stack.setChannel(n+1);
			run(LUT_list[n]);
			resetMinAndMax();
	}

	// Add rectangular bin as defined by script parameters
	if (bin_boolean == 1) {
		setBatchMode("show");
		makeRectangle((w-bin_width*scale)/2, (h-bin_height*scale)/2, bin_width*scale, bin_height*scale);
		waitForUser("Press OK When Finished", "(1) Use 'Selection Rotator' on toolbar \n(2) Click and drag to rotate the bin \n(3) ALT+click or SHFT+click to move the bin");
		run("Add Selection...");
	}

	// Add freehand polygon if defined by script parameter 
	if (freehand_boolean == 1) {
		setBatchMode("show");
		waitForUser("Press OK When Finished", "(1) Use the line tool or polygon tool \n(2) Distance is shown in the FIJI toolbar at the bottom \n(3) Press 'b' to add the line or polybox to the overlay");
		run("Add Selection...");
	}
	setBatchMode("hide"); // Hides image to go back to previously set BatchMode, in this case true

	// Rename image (if altered, such as with blinding), then save it with the binned prefix to clarify origin
	rename(title); 
	saveAs(output + File.separator + "binned_" + title);

	// Send to log window the iterations of the loop on a file, will only add loop++ if in the processFile of the processFolder function, then printing the image title nicely, as a sort of progress bar 
	print(loop++ + " out of " + count + " : " + title);
}


