#!/usr/bin/env python
import argparse
import os
from svglib.svglib import svg2rlg
from reportlab.graphics import renderPDF, renderPM


def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='rnastructureplot.py')
    parser.add_argument("-i", help='Input file path', required=True)
    parser.add_argument("-o", help='Output file path', required=True)
    args = parser.parse_args()
    outputpathsplit=args.o.split("/")
    print(outputpathsplit)
    filenameextension=outputpathsplit[-1]
    print(filenameextension)
    filenameextensionlist=filenameextension.split(".")
    filename=filenameextensionlist[0]
    first_line=">" + filename + "\n"
    infilehandle = open(args.i, 'r')
    intempfile=args.i + "_tmp"
    intempfilehandle=open(intempfile, 'w')
    count = 0
    for line in infilehandle:
        count += 1
        if count==1:
            intempfilehandle.write(first_line)
        else:
            intempfilehandle.write(line)
    # Closing files
    infilehandle.close()
    intempfilehandle.close()
    cmd1 = "RNAplot --infile " + intempfile + " -t 4 --output-format=svg "
    print(cmd1)
    returned_value = os.system(cmd1)    
    #cmd2 = "inkscape " + filename + "_ss.svg --export-png=" + filename + ".png"
    drawing = svg2rlg(filename + "_ss.svg")
    renderPM.drawToFile(drawing, filename + ".png", fmt="PNG")    
    #print(cmd2)
    #returned_value = os.system(cmd2)
    cmd3 = "mv " + filename + ".png " + args.o
    print(cmd3)
    returned_value = os.system(cmd3)
    
if __name__ == '__main__':
    main()
