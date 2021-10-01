
import pandas as pd


# COLOR CONVERSION
def rgb2hex(r,g,b):
    return "#{:02x}{:02x}{:02x}".format(r, g, b)


def convert_color(UCSC_COLORs, colors, color_dict, pt_class_color_dict, pt_class_patric_color_dict):
    file1 = open(UCSC_COLORs, 'r')
    lines = file1.readlines()
    hex_color_dict = {}
    patric_color_dict = {}
    for line in lines:
        my_str = line.strip().replace(' ', '')
        my_name_d = my_str.split('=')[0]
        rgb = my_str.split('=')[1].split(',')
        r_c = int(rgb[0])
        g_c = int(rgb[1])
        b_c = int(rgb[2])
        new_color = rgb2hex(r_c, g_c, b_c)
        if my_name_d in colors:
            hex_color_dict[(list(color_dict.keys())[list(color_dict.values()).index(my_name_d)])] = new_color
            patric_color_dict[(list(color_dict.keys())[list(color_dict.values()).index(my_name_d)])] = my_name_d

    (pd.DataFrame.from_dict(data=hex_color_dict, orient='index').to_csv(pt_class_color_dict, header=False))
    (pd.DataFrame.from_dict(data=patric_color_dict, orient='index').to_csv(pt_class_patric_color_dict, header=False))


def karyotype_writing(pathway_classes, multiplication_factor, count_pathway, color_dict, karyotype_str, karyotype_file_name):
    for element in pathway_classes:
        ID = "".join([word[0] for word in element.split()])
        start = 0
        end = multiplication_factor * int(count_pathway[element])
        temp_str = 'chr\t-\t' + ID + '\t' + ID + '\t' + str(start) + '\t' + str(end) + '\t' + color_dict[element] + '\n'
        karyotype_str += temp_str
    f_karyotype = open(karyotype_file_name, 'w')
    f_karyotype.write(karyotype_str)
    f_karyotype.close()
    print('KARYOTYPE WRITTEN to', karyotype_file_name)


def gene_start_end_add(file_name, start_list, end_list):
    with open(file_name, 'r') as alex_file:
        alex_lines = alex_file.readlines()
    for line in alex_lines:
        protein_info_list = line.split(' ')
        my_start = protein_info_list[1]
        my_end = protein_info_list[2]
        start_list.append(my_start)
        end_list.append(my_end)


def take_ids_from_file(file_name, delimiter):
    with open(file_name,'r') as my_protein_file:
        protein_lines = my_protein_file.readlines()
    prot_id_list = []
    for protein in protein_lines:
        protein_info_list = protein.split(delimiter)
        prot_id = protein_info_list[0]
        prot_id_list.append(prot_id)
    return prot_id_list


def ref_seq_id_found(protein):
    for protein_info in protein.split('\t'):
        for element in protein_info.split(' '):
            if 'genomic.gbff:' in element or 'sequence:' in element:
                return element.split(':')[-1]


def control_id_with_list(prot_id, prot_id_list):
    if prot_id in prot_id_list:
        return True
    else:
        return False


def gene_writing(id_special_list, protein_list_file, output_folder, the_organism):
    with open(protein_list_file,'r') as my_protein_file:
        protein_lines = my_protein_file.readlines()
    gene_string = ''
    d_gene_string = ''
    c_gene_string = ''
    spec_gene_names_str = ''

    for protein in protein_lines:
        protein_info_list = protein.split('\t')
        start = protein_info_list[0]
        end = protein_info_list[1]
        strand = protein_info_list[2]
        type_CDS = protein_info_list[3]
        the_prot_id = protein_info_list[4]
        gene_name = ''

        if protein_info_list[5] != 'NA':
            gene_name = protein_info_list[5]
        ref_seq_id = ref_seq_id_found(protein)
        #print(start, end, strand, type_CDS, the_prot_id, gene_name)
        #print(ref_seq_id)
        gene_string += the_organism + ' ' + start + ' ' + end + ' ' + type_CDS + '\n'
        if strand == 'D':
            d_gene_string += the_organism + ' ' + start + ' ' + end + ' ' + type_CDS + '\n'
        elif strand == 'C':
            c_gene_string += the_organism + ' ' + start + ' ' + end + ' ' + type_CDS + '\n'

        for element in id_special_list:
            if control_id_with_list(the_prot_id, element.list):
                element.mystr += the_organism + ' ' + start + ' ' + end + ' ' + gene_name + '\n'

        if gene_name != '':
            spec_gene_names_str += the_organism + ' ' + start + ' ' + end + ' ' + gene_name + '\n'

    open(output_folder + '/genes.txt', 'w').write(gene_string)
    open(output_folder + '/genes_d.txt', 'w').write(d_gene_string)
    open(output_folder + '/genes_c.txt', 'w').write(c_gene_string)

    for element in id_special_list:
        open(output_folder + '/' + element.name + '.txt', 'w').write(element.mystr)

    open(output_folder + '/genes_names.txt', 'w').write(spec_gene_names_str)
    print('Genes were written to the file in', output_folder)


class SpecialList:
    def __init__(self, str_name, mlist):
        self.name = str_name
        self.mystr = ''
        self.list = mlist


def plot_text_str():
    my_str = '<plot>\ntype = text\ncolor = black\nfile = data/Links_Tile/Text.txt\nr1 = 0.90r\nr0 = 0.85r\nshow_links = no\n'
    my_str += 'link_dims = 4p,4p,8p,4p,4p\nlink_thickness = 1p\nlink_color = black\nlabel_size = 12p\n'
    my_str += 'label_font = condensed\npadding  = 0p\nrpadding = 0p\n</plot>\n\n'
    return my_str


def plot_str(name, color):
    my_str = '<plot>\nfile = data/Links_Tile/' + str(name) + '.txt\nr1 = 0.85r\nr0 = 0.80r\ncolor = '+ color + '\nthickness = 50p'
    my_str += '\nshow_label = yes\nlayers_overflow=collapse\norientation = out\nstroke_thickness = 0.5\nstroke_color = black\n</plot>\n\n'
    return my_str


def link_create(pathway_classes, all_join, link_out_name, col_name_for_pt, start_list, multiplication_factor, organism_id, battery_tile_out_file, color_dict, battery_plot_str, battery_text_label_str, link_in_circ_conf_str, link_in_circ_conf):
    for element in pathway_classes:
        ID = "".join([word[0] for word in element.split()])
        selected_rows = all_join.loc[all_join[col_name_for_pt] == element]
        link_file_name_in = link_out_name + '/' + ID + 'in.txt'
        link_file_name_out = link_out_name + '/' + ID + 'out.txt'

        print(selected_rows)
        counter = 0
        link_str_in = ''
        link_str_out = ''

        print('starta_list', len(start_list))


        my_boom = selected_rows[selected_rows['Start'].isin(start_list)]
        my_boom_nego = selected_rows[~selected_rows['Start'].isin(start_list)]

        in_list = []
        out_list = []

        for index, row in my_boom.iterrows():
            link_str_in += ID + '\t' + str(counter) + '\t' + str(counter + multiplication_factor) + '\t'
            link_str_in += organism_id + '\t' + str(row['Start']) + '\t' + str(row['End']) + '\n'
            in_list.append((counter, counter + multiplication_factor))
            counter += multiplication_factor
        for index, row in my_boom_nego.iterrows():
            link_str_out += ID + '\t' + str(counter) + '\t' + str(counter + multiplication_factor) + '\t'
            link_str_out += organism_id + '\t' + str(row['Start']) + '\t' + str(row['End']) + '\n'
            out_list.append((counter, counter + multiplication_factor))
            counter += multiplication_factor

        battery_half_start = 0
        if len(in_list) == 0:
            battery_half_end = 1
        else:
            battery_half_end = in_list[-1][1]
        battery_all_start = 0
        battery_all_end = out_list[-1][1]

        battery_tile_out_filename_half = battery_tile_out_file + '/' + ID + 'half.txt'
        battery_tile_out_filename_all = battery_tile_out_file + '/' +  ID + 'all.txt'

        battery_tile_str_half = ID + ' ' + str(battery_half_start) + ' ' + str(battery_half_end) + '\n'
        battery_tile_str_all = ID + ' ' + str(battery_all_start) + ' ' + str(battery_all_end) + '\n'

        btr_file = open(battery_tile_out_filename_half, 'w')
        btr_file.write(battery_tile_str_half)
        btr_file.close()
        btr_file = open(battery_tile_out_filename_all, 'w')
        btr_file.write(battery_tile_str_all)
        btr_file.close()

        battery_plot_str += plot_str(ID + 'all', color_dict[element] + '_a2')
        battery_plot_str += plot_str(ID + 'half', color_dict[element])

        the_value_percent = 0
        if len(in_list) > 0:
            the_value_percent = int((battery_half_end / battery_all_end) * 100)

        battery_text_label_str += ID + ' ' + str(battery_half_start) + ' ' + str(battery_half_end) + ' %' + str(
            the_value_percent) + '\n'

        link_file = open(link_file_name_in, 'w')
        link_file.write(link_str_in)
        link_file.close()

        link_file = open(link_file_name_out, 'w')
        link_file.write(link_str_out)
        link_file.close()

        thickness = 1
        z_value = 1
        link_in_circ_conf_str += '<link>\n'
        link_in_circ_conf_str += 'file = ' + link_file_name_out + '\nribon = yes\nthickness = ' + str(
            thickness) + '\ncolor = ' + \
                                 color_dict[element] + '_a2'
        link_in_circ_conf_str += '\nradius = 0.80r' + '\nz = ' + str(z_value) + '\n</link>\n'

        thickness = 2
        z_value = 40
        link_in_circ_conf_str += '<link>\n'
        link_in_circ_conf_str += 'file = ' + link_file_name_in + '\nthickness = ' + str(thickness) + '\ncolor = ' + \
                                 color_dict[element]
        link_in_circ_conf_str += '\nradius = 0.80r' + '\nz = ' + str(z_value) + '\n</link>\n'

    link_in_circ_conf_str += '</links>\n'
    link_file_c = open(link_in_circ_conf, 'w')
    link_file_c.write(link_in_circ_conf_str)
    link_file_c.close()

    link_text = open(battery_tile_out_file + '/Text.txt', 'w')
    link_text.write(battery_text_label_str)
    link_text.close()

    text_plotting_str = plot_text_str()

    ####PLOT
    link_plot_file_end = open(battery_tile_out_file + '/plot.txt', 'w')
    link_plot_file_end.write(battery_plot_str + '\n' + text_plotting_str)
    link_plot_file_end.close()


def SSGR_find(clustering_factor_bp, clustering_factor_size, strain_specific, organism_id, SSGR, SSGRText):
    special_pd = pd.read_csv(strain_specific, delimiter='\s+', usecols=[1, 2],
                                    header=None)
    list_special = special_pd.values.tolist()
    cluster_size = 1
    list_of_clusters = []
    temp_end = list_special[0][1]
    temp_start = list_special[0][0]
    list_of_start_end = []
    pre_Name = 'SSGR-'
    cluster_counter = 1

    for element in range(1, len(list_special)):
        start = list_special[element][0]
        end = list_special[element][1]
        if abs(int(start) - int(temp_end)) < clustering_factor_bp or int(start) < int(temp_end):
            cluster_size += 1
            list_of_start_end.append((start, end))
        else:
            if cluster_size > clustering_factor_size:
                print(temp_end)
                list_of_clusters.append(((pre_Name + str(cluster_counter)), list_of_start_end))
                cluster_counter += 1
            list_of_id = []
            list_of_start_end = []
            cluster_size = 1
            list_of_start_end.append((start, end))

        temp_end = end

    print(list_of_clusters)
    print('Len: ', len(list_of_clusters))

    circos_file_list = []
    for cluster in list_of_clusters:
        circos_file_list.append((cluster[0], cluster[1][0][0], cluster[1][-1][1]))

    print(circos_file_list)

    cluster_circos_str = ''
    cluster_circos_str_unlabeled = ''

    for cluster in circos_file_list:
        cluster_circos_str += organism_id + ' ' + str(cluster[1]) + ' ' + str(cluster[2]) + ' ' + cluster[0] + '\n'
        cluster_circos_str_unlabeled += organism_id + ' ' + str(cluster[1]) + ' ' + str(
            cluster[2]) + ' dblue' + '\n'

    open(SSGR, 'w').write(cluster_circos_str_unlabeled)
    open(SSGRText, 'w').write(cluster_circos_str)


def tick_create(organism_id, out_tick):
    tick_str = ""
    tick_str += "show_ticks = yes\nshow_tick_labels = yes\n<ticks>\nchromosomes_display_default = no\nchromosomes = " + organism_id
    tick_str += "radius = dims(ideogram,radius_outer) + 5p\nlabel_offset = 5p\norientation = out\nlabel_multiplier = 1e-6\ncolor = black\n"
    tick_str += "<tick>\nspacing  = 1u\nsize  = 12p\ncolor = black\nthickness = 2p\nshow_label = yes\nlabel_size = 15p\nlabel_color = black\n"
    tick_str += """label_offset  = 3p\nformat   = %d\nsuffix = " Mbp"\n</tick>\n<tick>\nspacing  = 0.1u\nsize  = 12p\n"""
    tick_str += "color = black\nthickness   = 2p\nshow_label  = yes\nlabel_size  = 15p\nlabel_color  = black\nlabel_offset = 3p\n"
    tick_str += "format   = %.1f\n</tick>\n\n<tick>\nspacing  = 0.01u\nsize = 6p\ncolor  = dgrey\nthickness = 1p\nshow_label  = no\n</tick>\n</ticks>\n"
    open(out_tick, 'w').write(tick_str)
    print('Tick File was Created')

def idegram_create(out_ideo):
    ideo_str = ''
    ideo_str += "<ideogram>\n<spacing>\ndefault = 0.002r\n</spacing>\n"
    ideo_str += "radius = 0.85r\nthickness = 20p\nfill = yes\nstroke_thickness = 1\nstroke_color = black\n"
    ideo_str += "show_label = yes\nlabel_font = default\nlabel_radius = 1.035r\nlabel_size = 15\n"
    ideo_str += "label_parallel = no\n</ideogram>"
    open(out_ideo, 'w').write(ideo_str)
    print('Tick File was Created')


def conf_file_create(ticks, ideogram, link, plot, out_conf):
    conf_str = ''
    conf_str += "karyotype = data/reference.karyotype.txt\n"
    conf_str += "chromosomes_units = 1000000\n"
    if ticks == True:
        conf_str += "<<include ticks.conf>>\n"
    if ideogram == True:
        conf_str += "<<include ideogram.conf>>\n"
    if link == True:
        conf_str += "<<include data/Links/link.conf>>\n"
    if plot == True:
        conf_str += "<<include plots.conf>>\n"

    conf_str += "<image>\n<<include etc/image.conf>>\n</image>\n"
    conf_str += "<<include etc/colors_fonts_patterns.conf>>\n"
    conf_str += "<<include etc/housekeeping.conf>>"
    open(out_conf, 'w').write(conf_str)
    print('Conf File was Created')

def plot_template_plot(r1, r0, file, color, thickness, stroke_th, stroke_co, out):
    str_p = ""
    str_p += "<plot>\nfile = " + file + "\nr1 = " + r1 + "\nr0 = " + r0 + "\ncolor = " + color + "\n"
    str_p += "thickness = " + thickness + "\n" + "layers=1\nlayers_overflow=collapse\nshow_label=yes\norientation="+out+"\n"
    str_p += "stroke_thickness = " + stroke_th + "\nstroke_color=" + stroke_co + "\n</plot>\n\n\n"
    return str_p

def plot_template_text(r1, r0, file, color, size):
    str_p = ""
    str_p += "<plot>\ntype=text\nfile = " + file + "\nr1 = " + r1 + "\nr0 = " + r0 + "\ncolor = " + color + "\n"
    str_p += "show_links=no\nlink_dims = 4p,4p,8p,4p,4p\nlink_thickness = 1p\nlink_color = black\nlabel_size = " + size + "\n"
    str_p += "padding  = 0p\nrpadding = 0p\n</plot>\n\n\n"
    return str_p

def plot_create(link, SSGR, SSGRText, strainspecific, strain_color, species, species_color, genes_d, genes_c, phages, gene_clusters, gc, plt_out):
    plot_str = ""
    plot_str += "<plots>\ntype = tile\n"
    if link:
        plot_str += "<<include data/Links_Tile/plot.txt>>\n\n"

    #SSGR
    plot_str += plot_template_plot("1.09r", "1.07r", SSGR, "dblue", "25p", "0.5", "black", "out")

    plot_str += plot_template_text("1.15r", "1.10r", SSGRText, "black", "11p")



    plot_str += plot_template_plot("0.92r", " 0.88r", strainspecific, strain_color, "50p", "0", "grey", "out")

    plot_str += plot_template_plot("0.92r", " 0.88r", species, species_color, "50p", "0", "grey", "out")

    plot_str += plot_template_plot("0.88r", " 0.85r", genes_d, "green", "25p", "0", "black", "out")
    plot_str += plot_template_plot("0.85r", " 0.82r", genes_c, "green", "25p", "0", "black", "in")

    plot_str += plot_template_plot("0.96r", " 0.92r", phages, "purple", "25p", "0.1", "black", "center")

    for cluster in gene_clusters:
        plot_str += plot_template_plot("0.96r", " 0.92r", cluster, "orange", "25p", "0.1", "black", "center")
        plot_str += plot_template_text("1.10r", "1.00r", cluster, "black", "12p")

    plot_str += "<plot>\ntype = heatmap\nfile = " + gc + "\ncolor = spectral-5-div\nr1 = 0.98r\nr0 = 0.96r\n</plot>\n"
    plot_str += "<plot>\ntype = line\nfile = " + gc + "\ncolor = black\nr1 = 0.98r\nr0 = 0.96r\n</plot>\n"

    plot_str += "</plots>\n"

    open(plt_out, 'w').write(plot_str)
    print('Plot File was Created')
