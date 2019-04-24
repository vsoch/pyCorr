#!/usr/bin/env python2

"""

pyCorrPrint Component --> Template Matching Report Generator in Python

This python script takes a report produced by pyCorr.py and prints an HTML
report for visually seeing the matched images.

INPUT:
-h, --help      Print this usage
-r --report=    Text report file
-t --template=  Path to image to display for template
-o --oname=     Name for output folder in PWD

USAGE: python pyCorrPrint.py --report=thresh_zstat1.nii_bestcomps.txt --template=/path/to/image.png -oname=thresh1

OUTPUT: report.html and images in oname directory in pwd

"""

__author__ = "Vanessa Sochat (vsochat@stanford.edu)"
__version__ = "$Revision: 1.0 $"
__date__ = "$Date: 2011/08/25 $"
__license__ = "Python"

import os
import sys
import re
import shutil
import operator
import getopt

# USAGE ---------------------------------------------------------------------------------
def usage():
    print __doc__

# Print HTML report with motion charts for flagged subjects
def printHTML(output,result):
    if not os.path.isfile(output + "/report.html"):
        print "Creating results HTML report in " + output + "..."
        report = open(output + "/report.html",'w')
        report.write("<html>\n<body>\n<h1>pyCorr Result Report</h1>\n")
        report.write("<a href=\"result.txt\">Text Report</a>")
        report.write("<p><strong>Template: </strong>: " + output + "</p>\n")
	
	# Print the template image
        report.write("<img src=\"img/template.png\" />\"\n")
	
        # Cycle through list of results, print name and links to component images:
        report.write("<h1>Top Three Matched Components per Subject</h1>\n<p>")
        for res in result:
	    report.write("<p><strong>" + res[0] + "</strong></p>\n")
	    for i in [1,3,5]:
	        report.write("<img src=\"" + res[i] + "\" width=\"30%\" height=\"30%\" />\"")
	    report.write("<br /><br />\n")
            report.write("Matching Scores: <strong>1)</strong> " + str(res[2]) + " <strong>2)</strong> " + str(res[4]) + "<strong> 3) </strong> " + str(res[6])) 
            report.write("<br /><br />\n")
    

        report.write("</body>\n</html>")
        report.close()

# Reads single column text file, returns list
def readInput(readfile):
    result = []
    print "Reading input file " + readfile
    try:
        rfile = open(readfile,'r')
	for line in rfile:
	    line = line.rstrip("\n").rstrip(" ").rstrip()
	    sub,match1,val1,match2,val2,match3,val3 = line.split(" ")
            if sub not in ("SubID"):
                result.append([sub.rstrip(),match1.rstrip(),val1.rstrip(),match2.rstrip(),val2.rstrip(),match3.rstrip(),val3.rstrip().rstrip("\n")])
        rfile.close()
    except:
        print "Cannot open file " + readfile + ". Exiting"
        sys.exit()
    return result

# Get full paths for components and subject folders
def fullPaths(result):
    # result[n][0] is the subject ID
    # result[n][1] -- first match name, result[n][2] -- first match value
    # result[n][3] -- second match name, result[n][4] -- second match value
    # result[n][5] -- third match name, result[n][6] -- third match value
    count = 0
    
    # First fix the subject folder path
    for res in result:
        # First grab the path up to the subject .ica folder
	ica = re.compile(".ica")
	ica.search(res[0])
	found = ica.search(res[0])

	# Replace this path with the path to the subject report folder        
	result[count][0] = res[0][0:found.start()+4]
        
        # Now extract the thresh_zstat number and match to the report png image
        for i in [1,3,5]:
            znumexp = re.compile('[0-9]{1,2}')
            znum = znumexp.search(res[i])
            znumber = res[i][znum.start():znum.end()]
	    result[count][i] = "IC_" + znumber + "_thresh.png"
        count = count + 1

    # Return the result with all partial paths
    return result

def setupOut(output,tempimg,result,infile):
    # Create output directory, if doesn't exist
    if not os.path.exists(output):
        os.makedirs(output)
    if not os.path.exists(output + "/img"):
        os.makedirs(output + "/img")
    if os.path.isfile(infile):
        shutil.copy(infile,output + "/result.txt")
    if os.path.isfile(tempimg):
        shutil.copy(tempimg,output + "/img/template.png")
    else:
        print "Cannot find " + tempimg + ". Exiting!"
        sys.exit()

    # Copy each subject image into the image folder, number subjects 1 to N
    count = 0
    for res in result:
        for i in [1,3,5]:
            shutil.copy(res[0] + "/stats/report/" + res[i],output + "/img/" + str(count + 1) + res[i]) 
	    result[count][i] = "img/" + str(count + 1) + res[i]
	count = count + 1

    return result

# MAIN ----------------------------------------------------------------------------------
def main(argv):
    try:
        opts, args = getopt.getopt(argv, "hr:o:t:", ["help","report=","oname=","template="])

    except getopt.GetoptError:
        usage()
        sys.exit(2)
    
    # First cycle through the arguments to collect user variables
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()    
        if opt in ("-r","--report"):
            input1 = arg
        if opt in ("-o","--oname"):
            output = arg
        if opt in ("-t","--template"):
            tempimg = arg

    # Read in text file input
    rawres = readInput(input1)

    # Convert image paths to .png paths
    pathres = fullPaths(rawres)

    # Setup output folders and images
    linkres = setupOut(output,tempimg,pathres,input1)

    # Print HTML report
    printHTML(output,linkres)

if __name__ == "__main__":
    main(sys.argv[1:])
