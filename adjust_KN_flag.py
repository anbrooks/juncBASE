#!/lab/64/bin/python
# <Program Title>
# Author: Angela Brooks
# Program Completion Date:
# Modification Date(s):
"""<One line description>.
"""

import sys
import optparse 
import pdb
#############
# CONSTANTS #
#############

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
                          dest="input_table",
                          type="string",
                          help="""Input table. First column is the K (known) or N
                                  (novel) flag""",
                          default=None)
    opt_parser.add_option("--junctions",
                          dest="ref_junctions",
                          type="string",
                          help="""List of junctions that will be considered the
                                  \"known\" set. Format should be \"chr start end\" """,
                          default=None)
    opt_parser.add_option("--output",
                          dest="output_table",
                          type="string",
                          help="New table which contains the proper flag.",
                          default=None)

    (options, args) = opt_parser.parse_args()
	
    # validate the command line arguments
    opt_parser.check_required("--input")
    opt_parser.check_required("--junctions")
    opt_parser.check_required("--output")

    # Get set of junctions
    jcn_set = parseJcns(options.ref_junctions)

    input_table = open(options.input_table)
    output_table = open(options.output_table, "w")

    for line in input_table:
        line = formatLine(line)

        if line.startswith("#"):
            output_table.write(line + "\n")
            continue

        these_jcns = []

        line_list = line.split("\t")

        these_jcns.extend(parseJcnStr(line_list[5])) 
        these_jcns.extend(parseJcnStr(line_list[6])) 

        n_or_k = "K"
        for jcn in these_jcns:
            if jcn not in jcn_set:
                n_or_k = "N"
                break

        line_list[0] = n_or_k

        output_table.write("\t".join(line_list) + "\n")
        
    
			
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

def parseJcns(junction_file_name):

    jcn_set = set([])

    jcn_file = open(junction_file_name)

    for line in jcn_file:
        line = formatLine(line)

        line_list = line.split("\t")
        jcn_str = "%s:%s-%s" % (line_list[0],
                                line_list[1],
                                line_list[2])

        jcn_set.add(jcn_str)

    return jcn_set

def parseJcnStr(jcn_str):

    list1 = jcn_str.split(";")

    list2 = []

    for comp in list1:
        for comp2 in comp.split(","):
            if comp2 == "":
                continue
            list2.append(comp2)

    return list2

#################
# END FUNCTIONS #	
#################	
if __name__ == "__main__": main()
