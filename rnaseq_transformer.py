#!/usr/bin/python
'''
rnaseq_transformer.py

Author: Jeffrey Thompson

Purpose: Performs various transformations to an expression pcl file.

naive: log2 transformation
qn: quantile normalization, requires reference file to be specified
tdm: training distribution matching, requires reference file to be specified
relog: inverse log, followed by log2(x) + 1

The required argument to this script is always the expression file to transform.
Optional arguments are -t, for transformation type and -r for reference file.
The transformation type is one of the four values listed above. The reference
file is required for grt and tdm. It should be the expression file used as
training data.

The output of this script is written to stdout.
'''

import math
import argparse
import sys
import math_functions
from subprocess import call

# Specify the command line arguments that are allowed.
parser = argparse.ArgumentParser()
parser.add_argument("expfile", help="path to file with expression values to transform")
parser.add_argument("-t", "--transformtype",
                    help="naive|qn|tdm|relog",
                    default="naive")
parser.add_argument("-r", "--referencefile",
                    help="an expression file that the input will be transformed in reference to")

# Parse command line arguments.
args = parser.parse_args()

# First get the path to the file to log transform and the transformation type.
geneExpressionFileName = args.expfile
transformType = args.transformtype

if transformType == "qn":
    call(["Rscript", "/home/jeff/Repos/tdm2015test2/validation/code/data_utils/quantile_normalizer.R", geneExpressionFileName, args.referencefile])
    sys.exit()
if transformType == "tdm":
    call(["Rscript", "/home/jeff/Repos/tdm2015test2/validation/code/data_utils/tdm_normalizer.R", geneExpressionFileName, args.referencefile])
    sys.exit()

# Open the gene expression data.
geneExpressionFileHandle = open(geneExpressionFileName, 'r')

# Skip the first line (contains header).
# Obviously, this is a bit fragile and could be improved.
geneExpressionFileHandle.readline()

# Variable to hold the max expression value in a file.
fileMaxVal = 0

# Variable to hold all expression values.
expVals = []

# Variables to hold the min and max expression values in ref file and
# the first and third quartiles.
refMinVal = sys.maxint
refMaxVal = 0
refFirstQ = 0
refThirdQ = 0

# If the reference file was specified, then read it to find the min, max,
# first and third quartiles.
if args.referencefile != None:
    # Open the reference expression file.
    refExpressionFileHandle = open(args.referencefile, 'r')

    # Skip the first line.
    refExpressionFileHandle.readline()

    # Traverse the reference expression file and find the max and min values as
    # well as the IQR
    for line in refExpressionFileHandle:
        columns = line.split()
        for i in range(1, len(columns)):
            columns[i] = 0 if columns[i].strip() in ["NA","null"] else float(columns[i])

        lineMinMax = math_functions.getMinMax(columns[1:])
        expVals.extend(columns[1:])
        if lineMinMax[1] > refMaxVal:
            refMaxVal = lineMinMax[1]
        if lineMinMax[0] < refMinVal:
            refMinVal = lineMinMax[0]

    percentiles = math_functions.getPercentiles(expVals)
    refThirdQ = percentiles[2]
    refFirstQ = percentiles[0]
    expVals = []

# Traverse the gene expression file and find the highest expression value.
for line in geneExpressionFileHandle:
    columns = line.split()
    for i in range(1, len(columns)):
        columns[i] = 0 if columns[i].strip() in ["NA","null"] else float(columns[i])

    lineMax = math_functions.getMinMax(columns[1:])[1]
    expVals.extend(columns[1:])
    if(lineMax > fileMaxVal):
        fileMaxVal = lineMax

# Reset file position.
geneExpressionFileHandle.seek(0)

# Get file percentiles.
percentiles = math_functions.getPercentiles(expVals)
fileFirstQ = percentiles[0]
fileThirdQ = percentiles[2]

# Find the inter-quartile range
iqr = fileThirdQ - fileFirstQ

# Open the gene expression data, to perform transformation.
geneExpressionFileHandle = open(geneExpressionFileName, 'r')

# Print the header information to output.
print geneExpressionFileHandle.readline().rstrip()

# Traverse the gene expression file and log transform all non-zero values.
for i, line in enumerate(geneExpressionFileHandle):
    # Split a row into columns. Each row contains expression values for a gene,
    # while each column contains a sample's expression value for that gene.
    columns = line.split('\t')

    # Find the max expression value for this gene.
    for i in range(1, len(columns)):
        columns[i] = 0 if columns[i].strip() in ["NA","null"] else float(columns[i])

    maxVal = math_functions.getMinMax(columns[1:])[1]

    # Traverse the columns of the current row (samples for a single gene).
    for j, col in enumerate(columns):
        # If it is the first column, print it to the output.
        if j == 0:
            print col.rstrip(),
        # Otherwise, perform the transformation.
        else:
            # Seperate columns by tabs.
            print "\t",

            # Do a simple log2 transform of the value in this column.
            if transformType == 'naive':
                print math.log(1 + float(col), 2),

            # Perform TDM transformation.
            elif transformType == 'tdm':
                colVal = float(col)
                topscale = (2**refMaxVal- 2**refThirdQ)/(2**refThirdQ - 2**refFirstQ)
                bottomscale = (2**refFirstQ - 2**refMinVal)/(2**refThirdQ - 2**refFirstQ)

                upandout = fileThirdQ + topscale*iqr
                downandout = fileFirstQ - bottomscale*iqr

                if downandout < refMinVal:
                    downandout = refMinVal

                if colVal > upandout:
                    colVal = upandout
                elif colVal < downandout:
                    colVal = downandout

                print (math.log(
                    ((colVal - downandout) / (upandout - downandout)) * (2**refMaxVal - 2**refMinVal) + 2**refMinVal, 2)),

            # Relog the expression values so that log2(x) becomes log2(1 + x).
            elif transformType == 'relog':
                if col=="NA" or col=="null": col=0
                colVal = float(col)
                colVal = 2**colVal
                colVal = math.log(1 + colVal, 2)
                print colVal,

            elif transformType == 'unlog':
                if col=="NA" or col=="null": col=0
                colVal = float(col)
                colVal = 2**colVal
                print colVal,

            elif transformType == 'round':
                if col=="NA": col=0
                colVal = float(col)
                
                print "{0:.5f}".format(round(colVal,5)),


            else:
                raise Exception("Transformation type not found!")
    print

