# This code was a combination of editing from scripts on the Roary github page and the matplot lib tutorial page.

# links to sites that helped:
# roary: https://github.com/sanger-pathogens/Roary/blob/master/contrib/roary_plots/roary_plots.py
# matplotlib: https://matplotlib.org/stable/gallery/pie_and_polar_charts/pie_and_donut_labels.html#sphx-glr-gallery-pie-and-polar-charts-pie-and-donut-labels-py

import matplotlib.pyplot as plt
import numpy as np

# my_vars = argparse.ArgumentParser(description='Interactive viewer of scoary output')

# my_vars.add_argument('-c','--csv',dest='csv_file', type=str, required=False, help='give the absolute path to the file containing the data matrix')

# my_vars.add_argument('-N', '--skipped-columns', action='store',type=int,default=8,help='First N columns of panta\'s output to exclude [Default: 8]')

# my_vars.add_argument('--format',choices=('png','tiff','pdf','svg'),default='png',help='Output format [Default: png]')

# args = my_vars.parse_args()

# # function to take the different gene families and
# # Load panta
# panta = pd.read_csv(args.csv_file, low_memory=False)

# # Set index (group name)
# panta.set_index('Gene', inplace=True)

# # Drop the other info columns
# panta.drop(list(panta.columns[:args.skipped_columns-1]), axis=1, inplace=True)

# # Transform it in a presence/absence matrix (1/0)
# panta.replace('.{2,100}', 1, regex=True, inplace=True)
# panta.replace(np.nan, 0, regex=True, inplace=True)

# # Sort the matrix by the sum of strains presence
# idx = panta.sum(axis=1).sort_values(ascending=False).index
# panta_sorted = panta.loc[idx]

# core     = panta[(panta.sum(axis=1) >= panta.shape[1]*0.99) & (panta.sum(axis=1) <= panta.shape[1]     )].shape[0]
# softcore = panta[(panta.sum(axis=1) >= panta.shape[1]*0.95) & (panta.sum(axis=1) <  panta.shape[1]*0.99)].shape[0]
# shell    = panta[(panta.sum(axis=1) >= panta.shape[1]*0.15) & (panta.sum(axis=1) <  panta.shape[1]*0.95)].shape[0]
# cloud    = panta[panta.sum(axis=1)  < panta.shape[1]*0.15].shape[0]

fig, ax = plt.subplots(figsize=(10, 10), subplot_kw=dict(aspect="equal"))

#recipe = ["338 genes", "138 genes", "3000 genes", "200000 genes"]
#data = [338, 138, 3000, 200000]
# added from roary script
# Define custom labels for the legend
# custom_labels = ['core\n(99% <= strains <= 100%)',
#                 'soft-core\n(95% <= strains < 99%)',
#                 'shell\n(15% <= strains < 95%)',
#                 'cloud\n(strains < 15%)']

custom_labels = ['core',
                'soft-core',
                'shell',
                'cloud']

# pan-genome distirbution
data = [338, 183, 3793, 266906]

# due to time limitations we are manually adding the strain values, but keep in mind we used the roary code above to get how many of our strains contiributed to the different
# parts of the pan-genome. THeres just not enough time to incorporate it from the preliminary script
recipe = [f"338 genes\n(1161 <= strains <= 1173)", f"183 genes\n(1161 <= strains <= 1114))", f"3793 genes\n(175 <= strains <= 1114) ", f"266906 genes\n(strains < 175)"]

# choosing colorblind friendly hex codes
colors = ['#1f77b4', '#ff7f0e', '#6a3d9a','#6baed6']

# creating pie chart. Removed wedgeprops=dict(width=0.5) so it becomes a circle and not a donut
wedges, _ = ax.pie(data, labels=None, colors=colors,startangle=-40)

bbox_props = dict(boxstyle="round,pad=0.3", fc="w", ec="k", lw=0.72)
kw = dict(arrowprops=dict(arrowstyle="-"),
          bbox=bbox_props, zorder=0, fontsize=12)

# # setting the label positions for better layout
# label_positions = [
#     (-10, -100),  # for 338 core genes - slightly to the left of 138 genes
#     (140, -100),  # for 138 soft core genes - current position
#     (260, -100),  # for 3793 shell
#     (-100, 100), # 200000 cloud genes
# ]

# i printed the actual position the script was using, and then i adjusted it here so now the labels are being placed where i want
label_positions = [
    (1.35, -0.9956969511766536),  # for 338 core genes - slightly to the left of 138 genes
    (1.65, -0.7391873280956749),  # for 138 soft core genes - current position
    (1.48, -0.4384601830247148),  # for 3793 shell
    (-1.35, 0.8452108916217387), # 200000 cloud genes
]

for i, p in enumerate(wedges):
    ang = (p.theta2 - p.theta1) / 2. + p.theta1
    y = np.sin(np.deg2rad(ang))
    x = np.cos(np.deg2rad(ang))

    # find out the actual values being used for xytext
    pos = (1.35*np.sign(x), 1.4*y)
    print(pos)

    horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]

    # used to draw the different names for the pies. like the 338 genes and stuff. Also specifies where the boxes are placed, specifically the xytext=label positions.
    # thats why we created the label positions so we can control where they get placed. ,textcoords='offset points',xytext=label_positions[i]
    ax.annotate(recipe[i], xy=(x, y),
                xytext=label_positions[i], va='center',
                ha='center', **kw)

# Create legend based on custom labels and colors at top right - added by me
ax.legend(custom_labels,
          title="Gene Categories",
          loc="upper right",
          bbox_to_anchor=(1.3, 1))

ax.set_title("Pan-genome Composition")

# display the pie chart and from there you can do what you wish
plt.show()
