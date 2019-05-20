pyCorr Nifti Matching in Python

Vanessa Sochat
Update:
September 5, 2011

Overview

pyCorr is a script that matches user specified images to a template. The applied usage that I created it for was to match group functional networks created with group independent component analysis (gica, via a http://www.fmrib.ox.ac.uk/fsl/melodic/index.html but you could use it for any matching needs that you might have.

Scripts

PyCorr.py Takes a single column list of subjects (–input.txt: a column of complete paths to each subject folder that contains all images that are candidates for being matched) and a single template image to be matched to (–template=) and a list of images that each subject has in their specified folder (–images=) that should also be a single column text file. Lastly, specify an output folder (–output) where a text file with results will be generated. An overview of the method is as follows:

Read in images
Find spots that meet some user specified threshold, convert coordinates to MNI space, call this list “template hotspots”
Look up “template hotspots” in each image you are trying to match after converting MNI coordinates back to raw coordinate space
For shared activation, add one to a shared activation voxel count, and add the value to a shared activation total
Set all “template hotspots” equal to zero and then search for any remaining activation (activation in image but outside template)
Create an equivalent unshared activation voxel count and value total
Calculate an activation difference score
Rank scores in descending order and choose top three
I am still testing this algorithm, and it is subject to change!

# SUBMISSION ON COMMAND LINE
pyCorr.py --subs=input.txt --template=/path/to/template/image/thresh_zstat27.nii.gz --images=images.txt --output=/path/to/output/folder
PyCorrPrint.py Takes the text file output from pyCorr.py and produce a stand alone folder to throw up on a web server with the matching results. the –report tag should point to the output from pyCorr.py, the –template tag should be to some graphical representation of the template image, and the –oname tag should be the name that you want for the folder to contain the web report, which will be created in the present working directory (pwd).

# SUBMISSION ON COMMAND LINE
python pyCorrPrint.py --report=thresh_zstat2.nii_bestcomps.txt --template=../groupmelodic.ica/report/IC_2_thresh.png --oname=thresh2
Modification

The intended usage is with output from the Melodic+ Dual Regression package, and this shouldn't be an issue with pyCorr.py since the input file simply needs paths to folders with images to match. You can make the paths whatever you like, and the path to the template image can be whatever you like as well! The pyCorrPrint.py script, however, is expecting the names of the images to have a number in them (for example, zstat_thresh_1.nii.gz) and this number is used to look up the .png image in the original ica directory. Each .ica folder path in the output file from pyCorr.py should end in .ica, and the script expects to find the .png for use in the web report under /stats/report (this is FSL's standard that I adopted for the ica+ package) within the subject directory. Feel free to modify the script as needed… I threw it together in about two days so it's pretty unsophisticated!
