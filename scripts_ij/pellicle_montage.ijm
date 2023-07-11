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
	run("Set Scale...", "distance=68.4 known=1 unit=mm");
	
	setMinAndMax(300, 1600);
	run("Magenta");
	run("Next Slice [>]");
	setMinAndMax(300, 1400);
	run("Fire");
	run("Next Slice [>]");
	setMinAndMax(550, 6500);
	run("Grays");

	run("Make Montage...", "columns=3 rows=1 scale=0.25");

	selectWindow("Montage");
	rename(title_noext + "_montage");
	run("Scale Bar...", "width=5 height=10 font=24 color=White background=None location=[Lower Right] bold overlay");

	selectWindow(title);
	run("Scale Bar...", "width=5 height=20 font=42 color=White background=None location=[Lower Right] bold overlay");
	run("Split Channels");
	for(i = 0; i < 3; i++) {
		selectWindow("C" + i+1 + "-" + title);
		rename(title_noext + "_C" + i+1);
	}
	
output_parent = File.getParent(input);
input_dir = File.getName(input);
output = output_parent + File.separator + input_dir + "_processed";
if (!File.exists(output)) {
	File.makeDirectory(output);
}

// Save as .jpeg - can be changed to "png" or "tiff"
    while(nImages > 0){
	title = getTitle;
	saveAs("png", output + File.separator + title);
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