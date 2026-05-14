/** 
 * MACRO TO COLOCALIZE CELLS
 * Tim Monko - August 10, 2021 - timmonko@gmail.com - https://github.com/timmonko
 *  
 *  Unfortunately, I can not identify a way to accept an array as a script parameter, so you will need to also take care of the arrays below the parameter header
 *  
 * Scripting parameters (i.e #@ __ designed to create a UI to run process as a batch, see comments and UI notes for further detail)
 * While this is not the ideal way for optional scripting parameters to be used (e.g., the UI scaling from 1 to infinite channels interested to use
 * To add more channels, you would need to copy and paste the channel-specific code, and then add another function call
 */

#@ File(label = "Input directory", style = "directory") input
#@ File(label = "Output directory", style = "directory") output 
#@ String(label = "File suffix", value = ".tif") suffix 

#@ Boolean(label = "Check to keep ROI") ROI_boolean
#@ Boolean(label = "Co-occurence?") cooccurrence_boolean
#@ Integer(label = "Overlap Percentage", value = 35) overlap_perc
#@ Boolean(label = "Co-subtraction?") cosubtraction_boolean
#@ Boolean(label = "Remove slices?") sliceremove_boolean
#@ Integer(label = "Start slice") start_slice
#@ Integer(label = "End slice") end_slice

/**
 * Co-occurence variables. 
 * The object is the image that is recreated. 
 * The selector is tested to see if the Selector overlaps the with the Object greater than the overlap percentage.
 */
 

Object 	= newArray(	"ROR",
					"ROR"
					);
Selector= newArray(	"BRN2",
					"CTIP2"
					);

/**
 * Co-subtraction variables. 
 * Both images must have the same original object, since this uses the Image Calculator, not the Binary Feature Extractor 
 * The primary image is kept, while the subtractor image is subtracted from it, usually it will be a co-occurence image.  
 * 
 */

Primary = newArray("BRN2", "ROR", "ROR"
				   );

Subtractor = newArray("ROR_On_BRN2", "ROR_On_BRN2", "ROR_On_CTIP2"
					  );

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
	
	if (sliceremove_boolean == 1) {
		run("Slice Remover", "first=start_slice last=end_slice increment=1");
	}
	run("Stack to Images");
	
	if (cooccurrence_boolean == 1) {
		for (c = 0; c < Object.length; c++) {
			selectWindow(Object[c]);
			run("Duplicate...", " "); 
			rename("Object");
			selectWindow(Selector[c]);
			run("Duplicate...", " "); 
			rename("Selector");
			run("Binary Feature Extractor", "objects=Object selector=Selector object_overlap=overlap_perc");
			rename(Object[c] + "_On_" + Selector[c]);
			close("Object");
			close("Selector");
		}
	}

	if (cosubtraction_boolean == 1) {
	for (s = 0; s < Primary.length; s++) {
		imageCalculator("Subtract create", Primary[s], Subtractor[s]);
		rename(Primary[s] + "_Minus_" + Subtractor[s]);
		
	}
}

	

	run("Images to Stack");
	
	// Add the original active ROI back to the image
	if (ROI_boolean == 1) {
		roiManager("Select", 0);
		run("Add Selection...");
	}

	
	// Rename image (if altered, such as with blinding), then save it with the binned prefix to clarify origin
	rename(title); 

	// Save, and the send to log window the iterations of the loop on a file, will only add loop++ if in the processFile of the processFolder function, then printing the image title nicely, as a sort of progress bar 
	saveAs(output + File.separator + "coloc_" + title);
	print(loop++ + " out of " + count + " : " + title);
}


