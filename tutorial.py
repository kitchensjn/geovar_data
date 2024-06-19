import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pkg_resources
from geovar import *

data_path = pkg_resources.resource_filename("geovar", "data")

# Filepath to the VCF File
vcf_file = "{}/new_1kg_nygc.chr22.biallelic_snps.filt.n5000.vcf.gz".format(data_path)

# Filepath to the population panel file
population_panel = "{}/integrated_call_samples_v3.20130502.1kg_superpops.panel".format(data_path)

# Reading the population dataframe
pop_df = read_pop_panel(population_panel)

# Writing out VCF to a Frequency Table
af_df = vcf_to_freq_table(vcf_file, pop_df=pop_df, outfile="{}/test.freq.csv".format(data_path), minor_allele=True)

# Creating the GeoVar Object
geovar_test = GeoVar()

# Adding in the frequency file (all of it)
geovar_test.add_freq_mat(freq_mat_file="{}/test.freq.csv".format(data_path))

# Generate a geovar binning with the binning we used in our paper
geovar_test.geovar_binning()

geovar_plot = GeoVarPlot()

# Adding data directly from the geovar object itself
geovar_plot.add_data_geovar(geovar_test)

# Filter to remove very rare categories (only to speed up plotting)
geovar_plot.filter_data()

# Adding in a colormap (see code for alternative ideas beyond default)
geovar_plot.add_cmap()

# Shifting order of regional groupings to reflect OOA-demography
geovar_plot.reorder_pops(np.array(['AFR','EUR','SAS','EAS','AMR']))

fig, ax = plt.subplots(1,1,figsize=(3,6))
# The full plotting routine
geovar_plot.plot_geovar(ax);
ax.set_xticklabels(geovar_plot.poplist)

plt.show()