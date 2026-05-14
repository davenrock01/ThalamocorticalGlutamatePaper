// === Batch Merge Channels with Single Dialog for C5–C7 ===
// Assigns markers to C5 (cyan), C6 (magenta), C7 (yellow)
// Prompts all three marker inputs in one dialog.
// Author: Daven Rock, July 14, 2025

macro "Batch Merge Channels (One Dialog)" {

    // Folder selection
    inputDir = getDirectory("Choose folder with TIF files");
    if (inputDir == "") exit("No folder selected");

    // Create dialog for markers
    Dialog.create("Merge Channels (C5–C7)");
    Dialog.addString("Marker for Cy2 (cyan):", "", 20);
    Dialog.addString("Marker for Cy3 (magenta):", "", 20);
    Dialog.addString("Marker for Cy5 (yellow):", "", 20);
    Dialog.show();

    c5_marker = Dialog.getString();
    c6_marker = Dialog.getString();
    c7_marker = Dialog.getString();

    // Collect markers and assign channel codes
    markersArr = newArray();
    channelCodes = newArray();
    if (c5_marker != "") {
        markersArr = Array.concat(markersArr, newArray(toUpperCase(trim(c5_marker))));
        channelCodes = Array.concat(channelCodes, newArray("c5"));
    }
    if (c6_marker != "") {
        markersArr = Array.concat(markersArr, newArray(toUpperCase(trim(c6_marker))));
        channelCodes = Array.concat(channelCodes, newArray("c6"));
    }
    if (c7_marker != "") {
        markersArr = Array.concat(markersArr, newArray(toUpperCase(trim(c7_marker))));
        channelCodes = Array.concat(channelCodes, newArray("c7"));
    }

    nMarkers = markersArr.length;
    if (nMarkers < 2) exit("Please specify at least two markers.");

    // Prepare output folder
    outputDir = inputDir + "MERGED" + File.separator;
    File.makeDirectory(outputDir);

    // Gather all tif files
    files = getFileList(inputDir);
    baseNames = newArray();

    // Identify base names (SLx_BRx_SECx)
    for (i = 0; i < files.length; i++) {
        name = files[i];
        lname = toUpperCase(name);
        if (endsWith(lname, ".TIF") || endsWith(lname, ".TIFF")) {
            parts = split(name, "_");
            if (parts.length >= 4) {
                base = parts[0] + "_" + parts[1] + "_" + parts[2];
                found = false;
                for (b = 0; b < baseNames.length; b++)
                    if (baseNames[b] == base) found = true;
                if (!found) baseNames = Array.concat(baseNames, newArray(base));
            }
        }
    }

    // Merge per base name
    for (b = 0; b < baseNames.length; b++) {
        base = baseNames[b];
        presentFiles = newArray();
        skipBase = false;

        // Find each marker file
        for (m = 0; m < nMarkers; m++) {
            marker = markersArr[m];
            matched = "";
            for (f = 0; f < files.length; f++) {
                fName = files[f];
                if (startsWith(fName, base + "_")) {
                    suffix = substring(fName, lengthOf(base) + 1);
                    dot = lastIndexOf(suffix, ".");
                    if (dot > 0) suffixNoExt = substring(suffix, 0, dot);
                    else suffixNoExt = suffix;
                    if (toUpperCase(suffixNoExt) == marker) {
                        matched = inputDir + fName;
                        break;
                    }
                }
            }
            if (matched == "") {
                print("Skipping " + base + " — missing marker " + marker);
                skipBase = true;
                break;
            } else {
                presentFiles = Array.concat(presentFiles, newArray(matched));
            }
        }

        if (skipBase) continue;

        // Open and collect titles
        titles = newArray();
        for (i = 0; i < presentFiles.length; i++) {
            open(presentFiles[i]);
            wait(100);
            titles = Array.concat(titles, newArray(getTitle()));
        }

        // Build merge command (C5–C7 mapping)
        mergeCmd = "";
        for (i = 0; i < nMarkers; i++) {
            mergeCmd += channelCodes[i] + "=" + titles[i] + " ";
        }
        mergeCmd += "create";
        run("Merge Channels...", mergeCmd);

        // Save merged composite
        outName = outputDir + base + "_MERGED.tif";
        saveAs("Tiff", outName);
        print("Saved: " + outName);

        // Close all open images
        while (nImages > 0) {
            selectImage(nImages);
            close();
        }
    }

    print("Merging complete. Output in: " + outputDir);
}
