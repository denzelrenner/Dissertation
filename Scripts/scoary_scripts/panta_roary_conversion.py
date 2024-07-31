import pandas as pd
import os
import argparse
import logging

vars = argparse.ArgumentParser(description='Convert Panta Rtab file to a true Roary csv file')

vars.add_argument('-r','--rtab',dest='rtab_file', type=str, required=True, help='give the absolute path to the gene presence Rtab file produced by panta')

vars.add_argument('-od','--out_dir',dest='output_dir', type=str, required=True, help='give the absolute path to an output directory to put all the output files into')

args = vars.parse_args()

# # create log file directory if it doesnt exist
# if not os.path.isdir(f"{os.path.expanduser('~')}/log_files"):
#     os.makedirs(f"{os.path.expanduser('~')}/log_files")

# # # configure the logger we will be using
# logging.basicConfig(level=logging.DEBUG,format='%(asctime)s %(levelname)s %(message)s',handlers=[logging.FileHandler(f"{os.path.expanduser('~')}/log_files/panta_to_roary.log",mode='w')],datefmt='%H:%M:%S')

# # create a logger
# logger = logging.getLogger(__name__)

#  function to apply lambda transformation to numeric values only

# create output directory if it doesnt exist
if not os.path.isdir(args.output_dir):
    os.makedirs(args.output_dir)


# open input and output files
with open(args.rtab_file,'r') as rtab_file, open(f"{args.output_dir}/gene_presence_absence.Rtab.csv",'w') as out_file:

    # store all data in a varibale
    rtab_data = rtab_file.readlines()

    # write the first line of data to the output file
    out_file.writelines(rtab_data[0])

    # go through each line of data
    for line in rtab_data[1:]:

        # split up all data in the line
        columns = line.split('\t')

        # get the first column
        gene_col = columns[0]

        # put the rest of the data into a variable
        value_cols = columns[1:]

        # create variable for new value column
        new_value_cols = []

        # go through the values cols and change all values greater than 1 to 1
        for value in value_cols:

            # check if the value is greater than one or not
            if int(value.strip()) >= 1:
                new_value_cols.append('1')

            # if it is lower than one then just return itself
            elif int(value.strip()) < 1:
                new_value_cols.append('0')

        # add the line back together by first creating a new list with the gene name at the beginning
        combined_line = [gene_col] + new_value_cols

        # create a line to write. This will separate all values in our combined line into tabs, and then add a new line at the end as well
        line_to_write = '\t'.join(combined_line) + '\n'

        # write to output file
        out_file.write(line_to_write)

        if new_value_cols.count('1') >= 1162:
            print(f"OK WE FOUND CORE GENES {new_value_cols.count('1')}")

# convert to a csv and add columns for annotations and non unique gene name

# load in new Rtab file
df = pd.read_csv(filepath_or_buffer=f"{args.output_dir}/gene_presence_absence.Rtab.csv",sep='\t')

# check if the data is actually separated by lines
print(df.iloc[1])

# # Renaming the second column to 'Non unique gene names'
# df.columns.values[1] = 'Non-unique Gene name'

# # renaming the third column to Annotations
# df.columns.values[2] = 'Annotation'

# get the index of the first column
first_column_index = df.columns.get_loc('Gene')

# Create two new columns with values from original_column
df.insert(first_column_index + 1, 'Non-unique Gene name', df['Gene'])
df.insert(first_column_index + 2, 'Annotation', df['Gene'])

# change every value that is above 1 to the value 1
# df.apply(lambda x: 1 if (isinstance(x, int) and x > 1) else x)
# df.apply(transform_numeric)

# check if change was good
print(df.iloc[1])

# Convert it to a csv by just writing to output
df.to_csv(path_or_buf=f"{args.output_dir}/final_gene_presence_absence.Rtab.csv",index=False)
