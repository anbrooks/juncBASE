#!/lab/64/bin/python
# buildVirtualReference_median.py
# Author: Angela Brooks
# Program Completion Date:
# Modification Date(s):
# Copyright (c) 2011, Angela Brooks. anbrooks@gmail.com
# All rights reserved.

"""<One line description>.
"""

import sys
import optparse 
import pdb

import rpy2.robjects as robjects

r = robjects.r
#############
# CONSTANTS #
#############
DEF_THRESH = 25
NA = "NA"
#################
# END CONSTANTS #
#################


###########
# CLASSES #
###########
class OptionParser(optparse.OptionParser):
    """
    Adding a method for required arguments.
    Taken from:
    http://www.python.org/doc/2.3/lib/optparse-extending-examples.html
    """
    def check_required(self, opt):
        option = self.get_option(opt)

        # Assumes the option's 'default' is set to None!
        if getattr(self.values, option.dest) is None:
            print "%s option not supplied" % option
            self.print_help()
            sys.exit(1)


###############
# END CLASSES #
###############
 
########
# MAIN #	
########
def main():
	
    opt_parser = OptionParser()
   
    # Add Options. Required options should have default=None
    opt_parser.add_option("--input",
                          dest="input_file",
                          type="string",
                          help="""Input file containing exclusion and inclusion
                                  counts for each event.""",
                          default=None)
#   opt_parser.add_option("--output_median",
#                         dest="output_median",
#                         type="string",
#                         help="""Output file containing exclusion and inclusion
#                                 counts for each event plus an extra column
#                                 for a virtual reference.  This virtual
#                                 reference will contain the median exclusion
#                                 count and the median inclusion count.""",
#                         default=None)
    opt_parser.add_option("--output_median_ratio",
                          dest="output_median_ratio",
                          type="string",
                          help="""Output file containing exclusion and inclusion
                                  counts for each event plus an extra column
                                  for a virtual reference.  This virtual
                                  reference uses the median expression level
                                  for the event and uses the median PSI to
                                  decide how to distribute the exclusion,
                                  inclusion counts.""",
                          default=None)
    opt_parser.add_option("--output_median_psi",
                          dest="output_median_psi",
                          type="string",
                          help="""Additional, optional output file containing 
                                  PSI values for all samples with an extra
                                  column for the median PSI of the samples.
                                  Used as the virtual reference PSI.""",
                          default=None)
    opt_parser.add_option("--total_thresh",
                          dest="total_thresh",
                          type="int",
                          help="Threshold for total number of counts.",
                          default=DEF_THRESH)
    opt_parser.add_option("--remove_from_median",
                          dest="remove_from_median",
                          type="string",
                          help="""Comma separated list of samples that will not
                                  be used when calculating the median.""",
                          default=None)
    opt_parser.add_option("--weights",
                          dest="weights",
                          type="string",
                          help="""Comma separated list of weights given in the
                                  order of the samples in the table. Weights are
                                  used to create a weighted median. Default is
                                  equal weight for all samples.""",
                          default=None)

    (options, args) = opt_parser.parse_args()
	
    # validate the command line arguments
    opt_parser.check_required("--input")
#    opt_parser.check_required("--output_median")
    opt_parser.check_required("--output_median_ratio")

    input_file = open(options.input_file)

#    output_median = open(options.output_median, "w")
    output_median_ratio = open(options.output_median_ratio, "w")

    output_median_psi = None
    if options.output_median_psi:
        output_median_psi = open(options.output_median_psi, "w")

    total_thresh = options.total_thresh

    remove_from_median = []
    if options.remove_from_median:
        remove_from_median = options.remove_from_median.split(",")

    weights = None
    if options.weights:
        weights = map(float, options.weights.split(","))

        # Use R limma package
        try:
            r.library("limma")
        except:
            print """In order to use weighted median, please install the limma package from Bioconductor: 
                     http://www.bioconductor.org/packages/release/bioc/html/limma.html"""
            print """In R:\nsource("http://bioconductor.org/biocLite.R")\nbiocLite("limma")"""

    remove_from_median_idx = []

    weights = None
    if options.weights:
        weights = map(float, options.weights.split(","))

    for line in input_file:

        line = formatLine(line)

        line_list = line.split("\t")

        if line.startswith("#"):
            # Print header
            out_header_list = list(line_list[:11])
            out_header_list.append("virtual")
            out_header_list.extend(line_list[11:])

            for sample in remove_from_median:
                idx = None
                try:
                    idx = line_list[11:].index(sample)
                except ValueError:       
                    print "%s is not a sample.  Check --remove_from_median option." % sample
                    continue

                if idx is not None: 
                    remove_from_median_idx.append(idx)

#            output_median.write("\t".join(out_header_list) + "\n")
            output_median_ratio.write("\t".join(out_header_list) + "\n")

            if output_median_psi:
                output_median_psi.write("\t".join(out_header_list) + "\n")
            continue

            if weights:
                if len(weights) != len(line_list[11:])-1:
                    print "Weights for every sample needs to be given"
                    opt_parser.print_help()


        (median_virtual,
         median_ratio_virtual,
         median_psi) = getMedianVirtualReferences(line_list[11:],
                                                  total_thresh,
                                                  weights,
                                                  remove_from_median_idx)

#        outline_list_median = list(line_list[:11])
        outline_list_median_ratio = list(line_list[:11])
        outline_median_psi = list(line_list[:11])

#        outline_list_median.append(median_virtual)
        outline_list_median_ratio.append(median_ratio_virtual)
        outline_median_psi.append(median_psi)

#        outline_list_median.extend(line_list[11:])
        outline_list_median_ratio.extend(line_list[11:])

        for excl_incl_count in line_list[11:]:
            outline_median_psi.append(getPSI(excl_incl_count, total_thresh))

#        output_median.write("\t".join(outline_list_median) + "\n")
        output_median_ratio.write("\t".join(outline_list_median_ratio) + "\n")
        if output_median_psi:
            output_median_psi.write("\t".join(outline_median_psi) + "\n")
        
    input_file.close()
#    output_median.close()
    output_median_ratio.close()
    if output_median_psi:
        output_median_psi.close()
        
			
    sys.exit(0)

############
# END_MAIN #
############

#############
# FUNCTIONS #
#############
def formatLine(line):
    line = line.replace("\r","")
    line = line.replace("\n","")
    return line

def getMedianVirtualReferences(excl_incl_vals, total_thresh,
                               all_weights,
                               remove_from_median_idx):
    """
    Returns two different virtual counts:
    median_virtual - median(exclusion);median(inclusion)
    median_ratio = median total expression * (1.00-median(PSI)); median total expression * median(PSI)
    Also outputs medianPSI
    """
    excl_counts = []
    incl_counts = []
    total_expressions = []
    psis = []
    weights = [] 

    for i in range(len(excl_incl_vals)):
        
        if i in remove_from_median_idx:
            continue

        elems = excl_incl_vals[i]

        excl_count, incl_count = map(int, elems.split(";"))

        total_count = excl_count + incl_count
    
        if total_count == 0:
            continue

        # Only consider median for events that are above a threshold for
        # expression
        if total_count < total_thresh:
            continue

        total_expressions.append(total_count)

        excl_counts.append(excl_count)
        incl_counts.append(incl_count)

        psis.append(float(incl_count)/total_count)

        if all_weights:
            weights.append(all_weights[i])

    if excl_counts == []:
        return "NA", "NA", "NA"

    if all_weights:
        virtual_median = "%d;%d" % (int(round(r['weighted.median'](robjects.IntVector(excl_counts),
                                                                   robjects.FloatVector(weights))[0])), 
                                    int(round(r['weighted.median'](robjects.IntVector(incl_counts),
                                                                   robjects.FloatVector(weights))[0])))
        median_total = int(round(r['weighted.median'](robjects.IntVector(total_expressions),
                                                      robjects.FloatVector(weights))[0]))
        median_psi = r['weighted.median'](robjects.FloatVector(psis),
                                          robjects.FloatVector(weights))[0]

    else:
        virtual_median = "%d;%d" % (int(round(r['median'](robjects.IntVector(excl_counts))[0])), 
                                    int(round(r['median'](robjects.IntVector(incl_counts))[0])))
        median_total = int(round(r['median'](robjects.IntVector(total_expressions))[0]))

        median_psi = r['median'](robjects.FloatVector(psis))[0]

    virtual_median_ratio = "%d;%d" % (int(round(median_total*(1.00-median_psi))),
                                      int(round(median_total*(median_psi))))

    return virtual_median, virtual_median_ratio, "%.2f" % (median_psi * 100)
    
   
def getPSI(excl_incl_ct_str, total_thresh):

    excl_str, incl_str = excl_incl_ct_str.split(";")

    excl = float(excl_str)
    incl = float(incl_str)

    if excl + incl == 0:
        return NA

    if excl + incl < total_thresh:
        return NA

    psi = (incl/(incl + excl)) * 100

    psi_str = "%.2f" % psi

    return psi_str
#################
# END FUNCTIONS #	
#################	
if __name__ == "__main__": main()
