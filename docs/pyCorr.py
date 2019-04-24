#!/usr/bin/env python2

"""

pyCorr Component --> Template Matching in Python

This python script uses nibabel and numpy to read in a nifti image template, and match
a list of user specified images to the template.  The intended purpose is to match
components (individual functional networks) from an individual ica analysis to the 
most similar component from a group ica (group networks) analysis.  Inputs should be 
as follows:

INPUT:
-h, --help      Print this usage
-s --subs=      Single column text file w/ list of subject folders containing components
-t --template=  The template image to match, such as a group network
-i --images =   Single column text file with a list of component images in folders
-o --output=    Name of output folder.  If not specified, will use pwd

If you input a list of subjects longer than one, keep in mind that each should have the
corresponding component images in the designated folder.  Whether 3D or 4D, the first
timepoint will be used by default to extract data.  If an image's first timepoint is
empty, the script will try the second.  If the second is also empty, it will exit with
error, because there is something wrong with your template or image!

USAGE: python pyCorr.py --subs=sublist.txt --template=/path/to/image.nii.gz --images=imagelist.txt --output=/path/for/outfie

Intended usage is for one template for 1+ subjects with a list of component images.  

OUTPUT: (template_name)_bestcomps.txt w/ top 3 components for each subject

"""

__author__ = "Vanessa Sochat (vsochat@stanford.edu)"
__version__ = "$Revision: 1.0 $"
__date__ = "$Date: 2011/08/20 $"
__license__ = "Python"

import os
import sys
import nibabel as nib
import scitools.numpytools as scinu
import numpy as np
import operator
import getopt

# PYCORR------------------------------------------------------------------------------
class pyCorr:
    def __init__(self,imname):
        self.name = imname  # name of the image, as input
        self.path = None    # Full path to the image        
        self.img = None     # a nibabel Nifti object to hold image
        self.checkFile()

        self.xdim = 0
	self.ydim = 0
        self.zdim = 0       
        self.readDim()

        self.data = []      # the raw Y data 
        self.aff = []       # Affine transformation matrix
        self.readData()
        self.readAff()
                        
        self.XYZ = []       # XYZ coordinates to match raw data
        self.RCP = []       # "raw coordinate points"
        self.readXYZ()

    def __repr__(self):
        return self.name

# CHECK FILE 
    def checkFile(self):
        if os.path.isfile(self.name):
	    self.path = os.path.abspath(self.name)
            try: # Read in template image:
                self.img = nib.load(self.path)
            except:
                 print "Cannot read image data " + self.name + " exiting!"
                 sys.exit()
        else:
            print "Cannot find " + self.name + ". Check the path, and re-run."
            sys.exit()

# READ DATA
    def readDim(self):
        self.xdim = self.img.get_shape()[0]
        self.ydim = self.img.get_shape()[1]
        self.zdim = self.img.get_shape()[2]

# READ AFFINE TRANSFORMATION MATRIX
    def readAff(self):
        self.aff = scinu.mat(self.img.get_affine())
    
# Get affine matrix
    def getAff(self):
        return np.array(self.aff)

# Read raw image data
    def readData(self):
	# Read in all data to a temporary variable
        dataTEMP = self.img.get_data()
	# Check to see if we have a 4D image, represented by a 4th dimension
	if len(self.img.get_shape()) > 3:
            print self.name + " has more than one timepoint, will try using first by default."        		
	    if self.notEmpty(dataTEMP[:,:,:,0:1]):	
	        self.data = dataTEMP[:,:,:,0:1]
	    else:
	        # Once or twice I've seen melodic spit out a component with 2 TPs, the first empty, and the second correct.
	        print "Data at timepoint 1 is an empty or errored image.  Trying next timepoints..."
	        if self.img.get_shape()[3] < 2:
		    print "The template is a 4D file but only has one timepoint that cannot be used.  Exiting!"
		    sys.exit()
	        else:
		    # Here we are checking timepoint 2, which likely has the map.  We could continue checking timepoints
		    # if two is empty, but this probably means something hugely wrong with the image, and we should stop
		    # and alert the user
        	    if self.notEmpty(dataTEMP[:,:,:,1:2]):	
		        print self.name + " has empty first timepoint, using second."
	    	        self.data = dataTEMP[:,:,:,1:2]
	            else:
 		        print self.name + " is empty at both timepoints 1 and 2, and we cannot use it.  Exiting!"
	
        # Otherwise, we have a 3D image and only one set of datapoints to try	
	else:
	    # Make sure that we don't have an empty image
	    if self.notEmpty(dataTEMP):	
	        self.data = dataTEMP
	    else:
	        print self.name + " is empty and cannot be used as a template!  Exiting."
		sys.exit()

# Check if data is empty
    def notEmpty(self,data):
	for p in range(0,np.shape(data)[0]-1):
            for o in range(0,np.shape(data)[1]-1):
                for d in range(0,np.shape(data)[2]-1):
                    if data[p,o,d] != 0:
		        return True
	return False

# Get raw image data
    def getData(self):
        return self.data

# Get XYZ Coordinate Matrix (in MNI space)
    def getXYZ(self):
        return np.array(self.XYZ)

# Get RCP Coordinate from MNI
    def mnitoRCP(self,coord):
        xcor = (coord[0] - self.aff[0,3]) / self.aff[0,0]
	ycor = (coord[1] - self.aff[1,3]) / self.aff[1,1]
	zcor = (coord[2] - self.aff[2,3]) / self.aff[2,2]
	return [xcor,ycor,zcor]

# Get MNI Coordinate from RCP
    def rcptoMNI(self,coord):
        xcor = (coord[0] * self.aff[0,0]) + self.aff[0,3]
	ycor = (coord[1] * self.aff[1,1]) + self.aff[1,3]
	zcor = (coord[2] * self.aff[2,2]) + self.aff[2,3]
	return [xcor,ycor,zcor]

# Get Raw Coordinate Space Matrix
    def getRCP(self):
        return np.array(self.RCP)

# Read XYZ Coordinates 
    def readXYZ(self):
        # I couldn't get this method to work, but will keep to try again...
        # R, C, P = scinu.ndgrid(scinu.seq(1,3),scinu.seq(1,4),scinu.seq(1,5))
        
        # Create coordinate space based on x,y,z dimensions, and multiply by affine matrix
        # Examples shown if xdim = 3, ydim=4, zdim=5
        # Create R row variable [1 2 3 1 2 3...] ydim * zdim times
        Rrow = list(scinu.seq(1,self.xdim)) * (self.ydim*self.zdim)

	# Create C row variable [1 1 1 2 2 2 3 3 3 4 4 4 1 1 1 2 2 2...] zdim X
	Crow = []
	for y in range(1,self.ydim+1):
	  for x in range(0,self.xdim):
	    Crow.append(y)
	Crow = Crow * self.zdim
	
	# Create P row variable [ each of 1:zdim xdim*ydim times ]
	Prow = []
	for z in range(1,self.zdim+1):
	  holder = ([z] * self.xdim*self.ydim)
	  for i in holder:
	    Prow.append(i)

	# Create row of 1s of length zdim*xdim*ydim so we can multiply matrices
	onedim = [1] * self.xdim * self.ydim * self.zdim

	# Stack each row on top of one another
	self.RCP = np.vstack((Rrow,Crow,Prow,onedim))

	# Make it into a matrix
	self.RCP = scinu.mat(self.RCP)

	# Grab the first three rows of the affine transformation matrix (4th is for time)
	affXYZ = self.aff[0:3]

	# Multiply affine transformation matrix by coordinate data to go from coordinate --> MNI space
	self.XYZ = affXYZ * self.RCP

	# self.XYZ[:,0] contains first set of xyz coordinates (in column)
        # It is a 3xn matrix of XYZ locations returned (in mm) 

# RESULT------------------------------------------------------------------------------
class pyCorrRes:
    def __init__(self,output,filename):
        self.output = output      # output folder
        self.file = filename      # filename
        self.name = None
        self.fullpath = None
	self.setPath()
	self.writeHeader()

    def setPath(self):
        base,ext = os.path.splitext(os.path.basename(self.file))
	self.fullpath = self.output + "/" + base + "_bestcomps.txt"
        self.name = base

    def writeHeader(self):
        try:
	    fopen = open(self.fullpath,'w')
	    fopen.write("SubID Match1 Score1 Match2 Score2 Match3 Score3\n")
            fopen.close()
	except:
            print "Cannot write file " + self.fullpath + ". Exiting"
            sys.exit()
        
    def addResult(self,result):
	try:
	    fopen = open(self.fullpath,'a')
	    for entry in result:
		fopen.write(str(entry) + " ",)
	    fopen.write("\n")
            fopen.close()
	except:
            print "Cannot write file " + self.fullpath + ". Exiting"
            sys.exit()


# USAGE ---------------------------------------------------------------------------------
def usage():
    print __doc__

# Reads single column text file, returns list
def readInput(readfile):
    flist = []
    print "Reading input file " + readfile
    try:
        rfile = open(readfile,'r')
	for line in rfile:
	    flist.append(line.rstrip("\n").rstrip())
        rfile.close()
    except:
        print "Cannot open file " + readfile + ". Exiting"
        sys.exit()
    flist
    return flist

# Check that all component files exist for each subject
def checkInput(subinput,compinput):
   print "Checking for all components images for each subject..." 
   for sub in subinput:
       if sub:
           for comp in compinput:
               if not os.path.isfile(sub + "/" + comp):
                   print "Cannot find " + comp + " for " + sub + ". Exiting!"
                   sys.exit()
   print "All components for all subjects have been found!  Continuing analysis..." 
    

# MAIN ----------------------------------------------------------------------------------
def main(argv):
    try:
        opts, args = getopt.getopt(argv, "ht:s:i:o:", ["help","template=","subs=","images=","output="])

    except getopt.GetoptError:
        usage()
        sys.exit(2)
    
    # First cycle through the arguments to collect user variables
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()    
        if opt in ("-t","--template"):
            input1 = arg
        if opt in ("-i", "--images"):
            input2 = arg
	if opt in ("-s","--subs"):
            sublist = arg
        if opt in ("-o","--output"):
            output = arg


    # Get list of subject and component paths
    subfile = readInput(sublist)
    imgfiles = readInput(input2)

    # Check that all components exist for each subject
    checkInput(subfile,imgfiles)
        
    # Read in template image to pyCorr object, and get xyz and raw data
    Template = pyCorr(input1)
    tempXYZ = Template.getXYZ()     # returned as numpy ndarray
    tempDAT = Template.getData()    # returned as numpy ndarray
    tempRCP = Template.getRCP()     # Raw coordinate points
	
    # Prepare output file
    if not output:
        output = os.getcwd()
    Result = pyCorrRes(output,Template.name)        

    # TEMPLATE WORK ------------------------------------------------------------------------
    # Identify voxels that meet criteria.  Since the template is a thresholded network map,
    # for now we will identify all voxels greater than 0.  In the future this could allow for
    # user input of any intensity that matches with an atlas, etc.
    coords = []
    indexes = np.nonzero( tempDAT > 0 )    # indexes[]
					   # indexes[n][0] is x coordinate in mm 
					   # indexes[n][1] is y coordinate in mm
                                           # indexes[n][2] is z coordinate in mm

    # tempDAT(indexes[0][i],indexes[1][0],indexes[2][i],indexes[3][i]) would get the activation value
    # tempXYZ is a 3X(xdim*ydim*zdim) array of X Y Z coordinates, in MNI space.  (tempXYZ[:,column]) would be a set of points
    # tempDAT (the data) is a 3D, X by Y by Z matrix of raw data.  Since we have an array, you reference a point as tempDAT[x][y][z]
    # indexes is a 3 X N list of references to raw coordinates that meet a certain criteria in the data - index value = RCP points
    # To convert the index to MNI space then we need to multiply by the voxel size and add the translation from the affine matrix 
    
    # Cycle through all indexes and save list of MNI coordinates, for use with matching
    # NOTE - coordinates lookup in tempXYZ also tested, results were equivalent to 11th decimal point! 				
    for i in range(0,len(indexes[0])):
	# The index corresponds with raw coordinate space, so we need to convert each to MNI
	coords.append(Template.rcptoMNI([indexes[0][i],indexes[1][i],indexes[2][i]]))

    # COMPONENT IMAGE WORK --------------------------------------------------------------------------   
    # For each subject, compute the similarity score of all components belonging to subject
    # with the template (Adapted from Kaustubh Supekar doTemplateMatching.m, 2008)
    for subject in subfile:
        if subject:
		print "Computing similarity scores for subject " + subject    

		# Create a pyCorr object for each component to match, put in list
		components = []
		for img in imgfiles:
		    img_current = subject + "/" + img
		    try:
		        Match = pyCorr(img_current)
		        components.append(Match)
		    except: 
		        print "Problem with " + img + " for subject " + subject + ". Exiting!"
		        sys.exit()

		# Create dictionaries to hold all results for one template across components
		#activation_difference = {}      # Holds score with direction (+/-)
		activation_differenceabs = {}   # Holds absolute value of score, for ranking
		    
	
		# Cycle through components and...
		for com in components:
		      
		    # Set activation counter variables to zero
		    activation_in_roi = 0
		    activation_in_roiabs = 0			    
		    voxel_in_roi = 0
		    activation_out_roi = 0
		    activation_out_roiabs = 0
		    voxel_out_roi = 0
		    # Get the data in coordinate space
		    data = com.getData()
		      
		    # For each, take the coordinate list (in MNI) and convert to the raw coordinate space
		    coordsRCP = []
		    for point in coords:
		        coordsRCP.append(com.mnitoRCP(point))	      
		     
			# SHARED ACTIVATION
			# For each point, try to look it up.  If we query an index that doesn't exist, this means
			# we don't have data for that point, and we don't use it in our similarity calculation.
		    for point in coordsRCP:
			try:
			    # Get the activation value at the point
			    if data[point[0],point[1],point[2]] != 0:
				# If it isn't zero, then add to shared scoring
				activation_in_roiabs = activation_in_roiabs + abs(data[point[0],point[1],point[2]])
				activation_in_roi = activation_in_roi + data[point[0],point[1],point[2]]
		                voxel_in_roi = voxel_in_roi + 1
			except:
		            print "Coordinate " + str(point) + " is not in " + com.name + " for subject " + subject
		            print "...will not be included in similarity calculation!"

		       # ACTIVATION IN IMAGE NOT IN TEMPLATE
		       # Set the activation of each coordinate that we found to overlap as 0
		    for point in coordsRCP:
			try:
		            data[point[0],point[1],point[2]] = 0
		        except:
			    print "Coordinate " + str(point) + " is not in " + com.name + " for subject " + subject
		            print "This coordinate will not be included in similarity calculation!"
	      
		    # Cycle through the data, and find voxels that still have activation
		    for p in range(0,np.shape(data)[0]-1):
		        for o in range(0,np.shape(data)[1]-1):
		            for d in range(0,np.shape(data)[2]-1):
		                if data[p,o,d] != 0:
				    # Make sure the template includes the point before including it
				    try:
					checkpoint = tempDAT[p,o,d]
		                        activation_out_roiabs = activation_out_roiabs + abs(data[p,o,d])
		                        activation_out_roi = activation_out_roi + data[p,o,d]
		                        voxel_out_roi = voxel_out_roi + 1
				    except:
		                        print "Point " + str(p) + " " + str(o) + " " + str(d) + " has activation, but not present in template ROI.  Will not count!"
		                     
		    # Each subject will have an activation difference score for each component to the template. 
		    # I kept the absolute_difference (without absolute value) variable here in case it is wanted in the future.
		    if (voxel_in_roi == 0) and (voxel_out_roi == 0):
			#activation_difference[com.name] = 0
                        activation_differenceabs[com.name] = 0
		    elif voxel_in_roi == 0:
			#activation_difference[com.name] = (0 - (activation_out_roi/voxel_out_roi))          
			activation_differenceabs[com.name] = (0 - (activation_out_roiabs/voxel_out_roi))
                    elif voxel_out_roi == 0:
                        #activation_difference[com.name] = (activation_in_roi/voxel_in_roi)          
			activation_differenceabs[com.name] = (activation_in_roiabs/voxel_in_roi)
                    else:
			#activation_difference[com.name] = ((activation_in_roi/voxel_in_roi) - (activation_out_roi/voxel_out_roi))          
			activation_differenceabs[com.name] = ((activation_in_roiabs/voxel_in_roi) - (activation_out_roiabs/voxel_out_roi))
		
		# CHOOSING TOP RESULTS ----------------------------------------------------------------------------------	
		# When we finish cycling through the components, we want to find the top three matching (the most similar) components
		# We rank with activation_differenceabs, which is looking at absolute activation values (negative and positive Z ranked equally)
			
		topmatch = list(sorted(activation_differenceabs.iteritems(), key=operator.itemgetter(1))[-1])
		secondmatch = list(sorted(activation_differenceabs.iteritems(), key=operator.itemgetter(1))[-2])
		thirdmatch = list(sorted(activation_differenceabs.iteritems(), key=operator.itemgetter(1))[-3])

		# Add the information about the top three to the final results log, for printing
		resultitem = [subject,os.path.basename(topmatch[0]),topmatch[1],os.path.basename(secondmatch[0]),secondmatch[1],os.path.basename(thirdmatch[0]),thirdmatch[1]]	
		Result.addResult(resultitem)

if __name__ == "__main__":
    main(sys.argv[1:])
