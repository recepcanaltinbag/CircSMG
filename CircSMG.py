from configparser import ConfigParser
import pandas as pd
from pathlib import Path
from metaDataCreate import *

#Recep Can Altınbağ, 2021

parser = ConfigParser()
parser.read("config.conf")

#INPUTS (files)
genome_features = parser.get("inputs", "genome_features")
pathways_genes = parser.get("inputs", "pathways_genes")
pathways_names = parser.get("inputs", "pathways_names")

#SETTINGS


colors = parser.get("settings", "colors").split(',')
UCSC_COLORs = parser.get("settings", "UCSC_COLORs")
multiplication_factor = parser.getint("settings", "multiplication_factor")
organism_id = parser.get("settings", "organism_id")
genome_start = parser.get("settings", "genome_start")
genome_end = parser.get("settings", "genome_end")
genome_color = parser.get("settings", "genome_color")
karyotype_file_name = parser.get("settings", "karyotype_file_name")


#OUTPUTS

pt_class_color_dict = parser.get("outputs", "pt_class_color_dict")
pt_class_patric_color_dict = parser.get("outputs", "pt_class_patric_color_dict")



output_folder = parser.get("outputs", "output_folder")
Path(output_folder).mkdir(parents=True, exist_ok=True)



#PATRIC FILES JOIN
df_genome_feature = pd.DataFrame(pd.read_csv(genome_features))[['PATRIC ID', 'Start', 'End', 'Strand', 'Product', 'GO']]
df_pat_genes = pd.DataFrame(pd.read_csv(pathways_genes))[['PATRIC ID', 'Pathway Name', 'Product']]
df_path_names = pd.DataFrame(pd.read_csv(pathways_names))[['Pathway ID', 'Pathway Name',  'Pathway Class']]

left_join = pd.merge(df_pat_genes, df_genome_feature, on='PATRIC ID', how='left').drop_duplicates(subset=['PATRIC ID', 'Start','End'])
all_join = pd.merge(left_join, df_path_names, on='Pathway Name', how='left').drop_duplicates(subset=['PATRIC ID', 'Start','End'])
all_join = all_join.sort_values(by=['Pathway Class', 'Pathway Name'])

print(all_join)
count_pathway = all_join['Pathway Class'].value_counts()
total_metabolic = count_pathway.sum()
print(total_metabolic)
count_pathway = count_pathway.to_dict()
print(count_pathway)

pathway_classes = all_join['Pathway Class'].unique().tolist()
#sorting the classes due to their length
pathway_classes = sorted(pathway_classes, key=lambda k: count_pathway[k], reverse=True)
print(pathway_classes)
color_dict = dict(zip(pathway_classes, colors))


#KARYOTYPE FILE WRITING
karyotype_str = 'chr\t-\t' + organism_id + '\t' + organism_id + '\t' + str(genome_start) + '\t' + str(genome_end) + '\t' + genome_color + '\n'
karyotype_writing(pathway_classes, multiplication_factor, count_pathway, color_dict, karyotype_str, karyotype_file_name)

#SPECIAL GENE LISTS
special_ones = parser.get("inputs", "special_genes").split(',')
protein_list_file = parser.get("inputs", "protein_list_file")


id_special_list = []

for element in special_ones:
    id_special_list.append(SpecialList(element.split('.')[0].split('/')[1], take_ids_from_file(element, '\t')))

gene_writing(id_special_list, protein_list_file, output_folder, organism_id)

strain_specific = parser.get("inputs", "strain_specific")

start_list = []
end_list = []
gene_start_end_add(strain_specific, start_list, end_list)
print('start_list',len(start_list))

col_name_for_pt = parser.get("settings", "col_name_for_pt")
link_out_name = parser.get("outputs", "link_out_name")
link_in_circ_conf = parser.get("outputs", "link_in_circ_conf")
battery_tile_out_file = parser.get("outputs", "battery_tile_out_file")

Path(battery_tile_out_file).mkdir(parents=True, exist_ok=True)
Path(link_out_name).mkdir(parents=True, exist_ok=True)

link_in_circ_conf_str = '<links>\n'
battery_tile_str_half = ''
battery_tile_str_all = ''
battery_plot_str = ''
battery_text_label_str = ''

link_create(pathway_classes, all_join, link_out_name, col_name_for_pt, start_list, multiplication_factor, organism_id, battery_tile_out_file, color_dict, battery_plot_str, battery_text_label_str, link_in_circ_conf_str, link_in_circ_conf)
############################# LINKS WRITED ########################################################################


#####################BATTERY CODESS ###############################
####TEXT





###SSGRSs
clustering_factor_bp = parser.getint("settings", "clustering_factor_bp")
clustering_factor_size = parser.getint("settings", "clustering_factor_size")
special_genes_file = strain_specific

SSGR = parser.get("outputs", "SSGR")
SSGRText = parser.get("outputs", "SSGRText")
SSGR_find(clustering_factor_bp, clustering_factor_size, strain_specific, organism_id, SSGR, SSGRText)


#CIRCOS MAIN FILES
###CONF FILE
out_conf = parser.get("outputs", "out_circ_conf")
conf_file_create(ticks=True, ideogram=True, link=True, plot=True, out_conf=out_conf)
###TICKS FILE
out_tick = parser.get("outputs", "out_tick")
tick_create(organism_id, out_tick)
###IDEOGRAM FILE
out_ideo = parser.get("outputs", "out_ideo")
idegram_create(out_ideo)


##### PLOTTING FILE
phages = parser.get("inputs", "phages")
gc = parser.get("inputs", "gc")
gene_clusters = parser.get("inputs", "gene_clusters").split(',')
plt_out = parser.get("outputs", "plt_out")

strainspecific = parser.get("inputs", "strain_specific")
strain_color = parser.get("inputs", "strain_specific_color")
species = parser.get("inputs", "species")
species_color = parser.get("inputs", "species_color")

genes_d = output_folder + "/genes_d.txt"
genes_c = output_folder + "/genes_c.txt"

plot_create(True, SSGR, SSGRText, strainspecific, strain_color, species, species_color, genes_d, genes_c, phages, gene_clusters, gc, plt_out)

print('END')
print('Ready For Circos')
