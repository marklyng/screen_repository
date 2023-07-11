/*
 * Macro template to process multiple images in a folder
 */

#@ File (label = "Input directory", style = "directory") input
//#@ File (label = "Output directory", style = "directory") output
#@ String (label = "File suffix", value = ".czi") suffix

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
		if(endsWith(list[i], suffix))
			processFile(input, list[i]);
	}
}

function processFile(input, file) {
	// Do the processing here by adding your own code.
	// Leave the print statements until things work, then remove them.

	run("Bio-Formats Importer", 
		"open=["+input + File.separator + file+"] autoscale color_mode=Default open_files rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT"); //Import files
	
	title = getTitle();
	title_noext = getTitleStripExtension();
	getPixelSize(unit, pixelWidth, pixelHeight);
	scale = 1/pixelWidth*1000;
	run("Set Scale...", "distance="+scale+" known=1 unit=mm");
	getDimensions(width, height, channels, slices, frames);

	Table.open(input + "/contrast.txt");
	
	mag_min = Table.get("value", 0);
	mag_max = Table.get("value", 1);
	gfp_min = Table.get("value", 2);
	gfp_max = Table.get("value", 3);
	grey_min = Table.get("value", 4);
	grey_max = Table.get("value", 5);

	if(channels <= 1) {
		setMinAndMax(grey_min, grey_max);
		run("Grays");
		run("Scale Bar...", "width=5 height=20 font=42 color=White background=None location=[Lower Right] bold overlay");
		} 
		else {
	
				setMinAndMax(mag_min, mag_max);
				run("Magenta");
				run("Next Slice [>]");
				setMinAndMax(gfp_min, gfp_max);
				run("Fire");
				run("Next Slice [>]");
				setMinAndMax(grey_min, grey_max);
				run("Grays");
				
				run("Stack to Images");
				selectWindow(title_noext + "-0001");
				rename(title_noext + "_mkate");
				
				selectWindow(title_noext + "-0002");
				rename(title_noext + "_gfp");
				run("Calibration Bar...", "location=[Upper Right] fill=None label=White number=5 decimal=0 font=12 zoom=1 bold overlay");
				
				selectWindow(title_noext + "-0003");
				rename(title_noext + "_bf");
				run("Scale Bar...", "width=5 height=20 font=42 color=White background=None location=[Lower Right] bold overlay");				
		}
				
			output_parent = File.getParent(input);
			input_dir = File.getName(input);
			output = output_parent + File.separator + input_dir + "_processed";
			if (!File.exists(output)) {
				File.makeDirectory(output);
			}
		

// Save as .jpeg - can be changed to "png" or "tiff"
    while(nImages > 0){
	title = getTitle();
	run("Export SVG", "filename=" + title + " folder=[" + output + "] exportchannelsseparately=None interpolationrange=0.0 locksensitiverois=false");
	close();
	}

}

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