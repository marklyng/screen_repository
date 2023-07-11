/*
 * Macro template to process multiple images in a folder
 */

#@ File (label = "Input directory", style = "directory") input
#@ File (label = "Output directory", style = "directory") output
#@ String (label = "File suffix", value = ".tif") suffix

// See also Process_Folder.py for a version of this code
// in the Python scripting language.
setBatchMode(true);
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
		}
	}
	
	saveAs("Results", output + File.separator + "pellicle_intensities.csv");
}


function processFile(input, output, file) {
	// Do the processing here by adding your own code.
	// Leave the print statements until things work, then remove them.

	run("Bio-Formats", "open=[" + input + File.separator + file + "] autoscale color_mode=Default open_files rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");

	image = getTitleStripExtension();
	parent = File.getName(File.getParent(input));
	dir = File.getName(input);
	name = parent + "_" + dir + "_" + image;
	rename(name);
	
	makeOval(78, 78, 747, 747);
	
	run("Set Measurements...", "area mean standard integrated median display redirect=None decimal=3");
	run("Measure");
	run("Next Slice [>]");
	run("Measure");
	close(name);
	
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
  return t;
}