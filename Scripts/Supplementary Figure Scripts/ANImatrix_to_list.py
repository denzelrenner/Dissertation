import os
import logging
import argparse

vars = argparse.ArgumentParser(description='Interactive viewer of scoary output')

vars.add_argument('--output_directory',dest='output_dir', type=str, required=False, help='Specify absolute path to output directory')

vars.add_argument('-i',dest= 'ani_file', type=str,help='File containing data matrix of pairwise ANI values')

vars.add_argument('--cutoff',dest='cutoff',type=float,default=999,help='ANI values greater than or equal this value will not be plotted. Default is 999 to include all ANI values')

# vars.add_argument('--fastANI_file', type=str,help='File containing data matrix of pairwise ANI values using fastANI')

args = vars.parse_args()

# open the ANI file
with open(args.ani_file,'r') as ani_input,open(f"{args.output_dir}/ANI_values.txt",'w') as out_file:

    # write header for outfile
    out_file.write('ANI_Values\n')

    # store all lines in a variable. ignore the first line because it just has values
    ani_data = ani_input.readlines()[1:]

    print(len(ani_data[0].split('\t')))

    print(len(ani_data))

    # go through each line, split the line and ignore the first value in that list because it will be a genome name
    for line in ani_data:

        # store ani values in a variable
        ani_values = line.split('\t')[1:]

        #print(len(ani_values))

        # now loop through each value and check if it is greater than the cutoff
        for value in ani_values:

            if float(value.rstrip()) <= args.cutoff:

                # write to our file
                out_file.write(f"{float(value.rstrip())*100}\n")

            #else:
            #    print(f"Error for value {value}")


    # # fast ani header
    # fastaniout_file.write('ANI_Values\n')

    # # store data
    # fastani_data = fastani_input.readlines()

    # # loop through fast ani data
    # for line in fastani_data:

    #     ANI = float(line.split('\t')[2])

    #     fastaniout_file.write(f"{ANI/100}\n")



