#@ File (label = "Input directory", style = "directory") input
#@ String (label = "File suffix", value = ".czi") suffix
//setBatchMode(true);

run("Set Measurements...", "area mean standard integrated display redirect=None decimal=3");
processFolder(input);

// function to scan folders/subfolders/files to find files with correct suffix
function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i]))
			processFolder(input + File.separator + list[i]);
		if(endsWith(list[i], suffix))
			processFile(input, list[i]);
	}
}

saveAs("Results", input + File.separator + "Results.csv");

function processFile(input, file) {
	run("Bio-Formats Importer", 
		"open=["+ input + File.separator + file +"] autoscale color_mode=Default open_all_series rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT"); //Import .czi-files
	
	title = getTitleStripExtension();
	title_ext = getTitle();
	
getDimensions(width, height, channels, slices, frames);
if(channels < 3) {
	print(title + " does not contain 3 channels");
	close();
	}
else {

	// Change the rectangle to include your image
	
	run("Duplicate...", "duplicate channels=1");
	selectWindow(title + "-1.czi");
	run("Gaussian Blur...", "sigma=2 slice");
	
	setAutoThreshold("Otsu dark");
	setOption("BlackBackground", false);
	run("Convert to Mask", "method=Otsu background=Dark calculate");

	
	run("Options...", "iterations=3 count=1 pad do=Close");
//	run("Options...", "iterations=3 count=1 pad do=Dilate");

	run("Analyze Particles...", "size=100000-Infinity include add"); //Remove particles below colony size (if you have some big straggler blobs that "Close" can't remove)
	close(title + "-1.czi");
//	run("Create Selection");
//	roiManager("Add");
	
	selectWindow(title_ext);
	Stack.setChannel(1);
	roiManager("Select", 0);
	run("Measure");
	roiManager("Deselect");

	Stack.setChannel(2);
	roiManager("Select", 0);
	run("Measure");
	
	close(title_ext);
	roiManager("Deselect");
	roiManager("Delete");
	}
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