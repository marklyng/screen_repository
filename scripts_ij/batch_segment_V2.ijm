/*
 * Batch segment multiple directories of .flex files and divide into separate binary channels as .tifs
 */

#@ File (label = "Input directory", style = "directory") input
//#@ File (label = "Output directory", style = "directory") output
#@ String (label = "File suffix", value = ".flex") suffix

processFolder(input);

// function to scan folders/subfolders/files to find files with correct suffix
function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i]))
			processFolder(input + File.separator + list[i]);
		if(endsWith(list[i], suffix)) {
			parent_dir = File.getParent(input);
			output = parent_dir + File.separator + "bins";
			File.makeDirectory(output);
			processFile(input, output, list[i]);
		}
	}
}

function processFile(input, output, file) {
	setBatchMode(true);

	run("Bio-Formats Importer", 
		"open=["+input + File.separator + file+"] autoscale color_mode=Default open_files rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT"); //Import files

	run("Convolve...", "text1=[1 1 1 1 1\n1 1 1 1 1\n1 1 1 1 1\n1 1 1 1 1\n1 1 1 1 1\n] normalize stack"); // Noise reduction
	run("Median 3D...", "x=0 y=0 z=2"); // Remove single cells present in only one slice
	
	run("Subtract Background...", "rolling=100 stack"); 

	// Segment with MaxEntropy
	/*setAutoThreshold("MaxEntropy dark");
	setOption("BlackBackground", true);
	run("Convert to Mask", "method=MaxEntropy background=Dark calculate black");
	*/
	run("Auto Threshold", "method=MaxEntropy white stack use_stack_histogram");
	
	title = getTitleStripExtension();
	print("Processing " + input + File.separator + title);

	// Delete already segmented files in a bin-directory
	if(File.exists(output + File.separator + title + "_C0_T0.tif")) {
		print("Segmented file exists, deleting existing file");
		del1 = File.delete(output + File.separator + title + "_C0_T0.tif");
		del2 = File.delete(output + File.separator + title + "_C1_T0.tif");
	}

	// Does not handle "." in the names of files!!
	run("OME-TIFF...", "save=["+output + File.separator + title +".tif] write_each_timepoint write_each_channel use export compression=Uncompressed");
	close("*");
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