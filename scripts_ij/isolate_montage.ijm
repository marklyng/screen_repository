/*
 * Macro template to process multiple images in a folder
 */

#@ File (label = "Input directory", style = "directory") input
#@ File (label = "Output directory", style = "directory") output
#@ String (label = "File suffix", value = ".flex") suffix

// See also Process_Folder.py for a version of this code
// in the Python scripting language.

processFolder(input);

// function to scan folders/subfolders/files to find files with correct suffix
function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i]))
			processFolder(input + File.separator + list[i]);
		if(endsWith(list[i], suffix)) {
			processFile(input, output, list[i]);
			dir = File.getParent(input + File.separator + list[i]);
			dir = replace(dir, ".*\\\\", "");
		}
	}
	while (nImages > 0) {
	run("Images to Stack", "name=Stack title=[] use");
	run("Make Montage...", "columns=4 rows=3 scale=0.25");
	run("Scale Bar...", "width=50 height=4 font=14 color=White background=None location=[Lower Right] bold overlay");
//	print(File.getDirectory(list[i]));
	saveAs("jpg", output + File.separator + dir + "_montage");
	close("*");
	}
}

function processFile(input, output, file) {
	// Do the processing here by adding your own code.
	// Leave the print statements until things work, then remove them.
	setBatchMode(true);
run("Bio-Formats Importer", 
		"open=["+input + File.separator + list[i]+"] autoscale color_mode=Composite open_files rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
title = getTitle();
//print(title);

run("Z Project...", "projection=[Max Intensity]");
selectWindow(title);
close();
title = getTitle();

run("Enhance Contrast", "saturated=0.35");
run("Magenta");
run("Next Slice [>]");
run("Enhance Contrast", "saturated=0.35");
run("Green");
run("RGB Color");

selectWindow(title);
close();
}


// Use this function to strip any number of extensions
// off images.
// Returns the title without the extension.
//====================================
function getTitleStripExtension() {
  t = getTitle();
  t = replace(t, ".tif", "");        
  t = replace(t, ".tiff", "");      
  t = replace(t, ".lif", "");      
  t = replace(t, ".lsm", "");    
  t = replace(t, ".czi", "");      
  t = replace(t, ".nd2", "");
  t = replace(t, ".flex", "");        
  return t;
}