/*	StackCellDetection_3Ch_Puncta_V1.IJM
	***********************
	Author: 		Titia+ Julia

 	1) detect cells in channel 4
 	2) measure intensity and other params in the Blue (PV) / Red (mCherry)/ Green (fos) channel
 	3) increase ROI and measure intensity and other params in the Blue Channel= Blue large	
 	
 	**********************************************

		Variable initiation

	**********************************************

*/
//	Arrays
var list 	= newArray();

//	Strings and numbers
var dir = "";						//	directory
var outDir = "";					//	output directory
var index = 0;						//	Index of folder (well)

var minSize = 80, maxSize = 250;			//	min and max size of nuclei for segmentation (in pixels)
var minCirc = 0;			
var enlarge = -1;
var teller=0;
/*
 	**********************************************

		Signal Measurements

	**********************************************
*/

macro "[M] Measure Nuclear Foldings"
{
	setup();
	setBatchMode(true);
	print("Analysis");
	print("************************************************************");
	index = 0;
	for(i=0; i<list.length; i++)
	{

		
		path = dir + list[i];
		if(endsWith(path,'.tif'))
		{
			teller=teller+1;
			run("Collect Garbage");
			open(path);
			id = getImageID;
			title = getTitle;
			print(title);
			prefix = substring(title,0,indexOf(title,".tif"));
			print("Image",i,":",prefix);			
						
			selectImage(id);
			run("Duplicate...", "title=PINKchannel.tif duplicate channels=4");// now gets only cells in channel 4
			PINKid2 = getImageID;			
			
			selectImage(id);
			run("Duplicate...", "title=PINKchannel.tif duplicate channels=4");
			PINKid = getImageID;

			selectImage(PINKid);
			roiManager("reset");
			run("Despeckle","stack");
			run("Gaussian Blur...", "sigma=1 stack");
			setAutoThreshold("Huang"+" dark stack"); // change to filter that you want, that selects the most cells - Li, IJ_isodata or Isodata, Huang
			run("Convert to Mask", "stack");//
			run("Fill Holes", "stack");
			run("Watershed", "stack");
			run("Analyze Particles...", "size="+80+"-"+250+" circularity="+0.5+"-1.00 show=Nothing exclude add stack");// could also add a negative enlarge in here
			n = roiManager("count");
			print(n, "Nuclei");
			roiManager("deselect");
			run("Clear Results");
			run("Set Measurements...", "area mean standard min centroid center bounding shape median stack redirect=None decimal=4");
			selectImage(PINKid2);
			roiManager("Measure");
			selectWindow("Results");
			saveAs("Measurements",dir+prefix+"_Results_PINK.txt");
			roiManager("Save", dir+prefix+"_Roiset_PINK.zip");
			
			//measure intensity values in the Red channel, these are the mCherry +/- cells		
			selectImage(id);
			ns = nSlices;
			run("Duplicate...", "title=Redchannel.tif duplicate channels=3 slices=1-"+ns);
			REDid = getImageID;
			titleRED = getTitle;
			print(titleRED);
			roiManager("deselect");
			run("Clear Results");
			roiManager("Measure");
			saveAs("Measurements",dir+prefix+"_Results_RED.txt");
			selectImage(REDid);
			close;

			// measure intensity values in the Blue channel, these are PV +/- cells
	  		selectImage(id);
			ns = nSlices;
			run("Duplicate...", "title=BLUEchannel.tif duplicate channels=1 slices=1-"+ns);
			BLUEid = getImageID;
			titleBLUE = getTitle;
			print(titleBLUE);
			roiManager("deselect");
			run("Clear Results");
			roiManager("Measure");
			saveAs("Measurements",dir+prefix+"_Results_Blue.txt");			

			// measure intensity values in the Green channel, these are Fos +/- cells
	  		selectImage(id);
			ns = nSlices;
			run("Duplicate...", "title=GREENchannel.tif duplicate channels=2 slices=1-"+ns);
			GREENid = getImageID;
			titleGREEN = getTitle;
			print(titleGREEN);
			roiManager("deselect");
			run("Clear Results");
			roiManager("Measure");
			saveAs("Measurements",dir+prefix+"_Results_GREEN.txt");	
			selectImage(GREENid);
			close;
			
			// Binarize the blue channel and measure signal in the blue channel in the predefined ROIs
			selectImage(id);
			run("Duplicate...", "title=BLUEchannel.tif duplicate channels=1");
			BLUEid2 = getImageID;			
			selectImage(id);
			run("Duplicate...", "title=PINKchannel.tif duplicate channels=1");
			BLUEid3 = getImageID;
			selectImage(BLUEid3);
			run("Make Binary", "method=Triangle background=Dark calculate black");; // change to filter that you want. Check some images
			roiManager("deselect");
			run("Clear Results");
			roiManager("Measure");
			saveAs("Measurements",dir+prefix+"_Results_BlueThreshold.txt");
				

			
		    //enlarge the ROIS for channel 4
			n=roiManager("count");//number of ROIS		
			print(n);
				if (roiManager("count")>0)
				{
				for(r=0; r<n; r++)
					{						
						roiManager("deselect");
						selectImage(PINKid);
						run("Duplicate...","title=PINKchannel.tif");
						PINKid2 = getImageID;
		
						selectImage(PINKid2);
						roiManager("deselect");
						roiManager("select", r);				
						run("Enlarge...", "enlarge=3");// choose size for enlarging)
						roiManager("update");
					}
				}

			//measure the PV intensity for the enlarged ROIs					
			roiManager("deselect");
			run("Clear Results");
			run("Set Measurements...", "area mean standard min centroid center bounding shape median stack redirect=None decimal=4");
			selectImage(BLUEid);
			roiManager("Measure");
			selectWindow("Results");
			saveAs("Measurements",dir+prefix+"_Results_BlueLarge.txt");
			
			//measure the binarized PV pixels for the enlarged ROIs
			roiManager("deselect");
			run("Clear Results");
			run("Set Measurements...", "area mean standard min centroid center bounding shape median stack redirect=None decimal=4");
			selectImage(BLUEid3);
			roiManager("Measure");
			selectWindow("Results");
			saveAs("Measurements",dir+prefix+"_Results_BlueThresholdLarge.txt");

			
//			//threshold the blue channel and measure signal inside enlarged
//			selectImage(BLUEid);
//			run("Duplicate...", "title=BLUEchannel.tif duplicate channels=4");
//			BLUEid3 = getImageID;			
//			selectImage(BLUEid3);
//			roiManager("reset");
//			setAutoThreshold("Triangle"+" dark stack") // change to filter that you want. Check some images
//			roiManager("deselect");
//			run("Clear Results");
//			roiManager("Measure");
//			saveAs("Measurements",dir+prefix+"_Results_BlueThresholdLarge.txt");
//			selectImage(Blueid3);
//			close;
//	

			selectImage(id);
			close;
			selectImage(PINKid);
			close;
			selectImage(PINKid2);
			close;
			selectImage(BLUEid);
			close;	
			selectImage(BLUEid2);
			close;
			selectImage(BLUEid3);
			close;
			
		}
	}
	print("************************************************************");
	setup();

}

/*	
 	**********************************************

		Functions

	**********************************************
*/


function getMoment()
{
     MonthNames = newArray("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec");
     DayNames = newArray("Sun", "Mon","Tue","Wed","Thu","Fri","Sat");
     getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
     TimeString ="Date: "+DayNames[dayOfWeek]+" ";
     if (dayOfMonth<10) {TimeString = TimeString+"0";}
     TimeString = TimeString+dayOfMonth+"-"+MonthNames[month]+"-"+year+"\nTime: ";
     if (hour<10) {TimeString = TimeString+"0";}
     TimeString = TimeString+hour+":";
     if (minute<10) {TimeString = TimeString+"0";}
     TimeString = TimeString+minute+":";
     if (second<10) {TimeString = TimeString+"0";}
     TimeString = TimeString+second;
     return TimeString;
}

function setup()
{
	print("\\Clear");
	run("Close All");
	run("Clear Results");
	roiManager("reset");
	run("Collect Garbage");
	run("Colors...", "foreground=white background=black selection=yellow");	
	setOption("BlackBackground", false);
	run("Set Measurements...", "area mean shape stack redirect=None decimal=4");
	dir = getDirectory("");
	list = getFileList(dir);
	isWin=indexOf(getInfo("os.name"),"Windows")>=0;
	TimeString = getMoment();
	print(TimeString);
	print("************************************************************");
}


function segment(2)
{	
	selectImage(PINKid2);
	roiManager("reset");
	run("Despeckle");
	run("Gaussian Blur...", "sigma=1");
	setAutoThreshold("Li"+" dark");
	run("Convert to Mask");
	run("Fill Holes");
	run("Watershed");
	run("Analyze Particles...", "size="+80+"-"+250+" pixel circularity="+0+"-1.00 show=Nothing exclude clear add");
	n = roiManager("count");
	print(n, "Nuclei");
	return n;
}