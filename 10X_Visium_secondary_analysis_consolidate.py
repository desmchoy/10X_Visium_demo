#!/usr/bin/env python3

import os
import re
import argparse
import glob
import pandas as pd
import collections
import xlsxwriter
from xlsxwriter.utility import xl_rowcol_to_cell


workdir = os.getcwd()
project = os.path.basename(workdir)

help_text = '10X_Visium_secondary_analysis_consolidate.py: Consolidating 10X Visium secondary analysis outputs'

parser = argparse.ArgumentParser(description=help_text)

parser.add_argument(
    dest='samplename',
    action='store',
    type=str,
    help='Sample Name (required)',
    metavar='samplename'
    )

arguments = parser.parse_args()
samplename = arguments.samplename

output = pd.ExcelWriter(samplename+'_secondary_analysis.xlsx', engine='xlsxwriter')
workbook = output.book
center = workbook.add_format({'align': 'center'})
center_bold_format = workbook.add_format({'align': 'center', 'bold': 1, 'valign': 'vcenter'})
bold_format = workbook.add_format({'bold': 1, 'valign': 'vcenter'})
italic_center_format = workbook.add_format({'italic': 1, 'align': 'center', 'valign': 'vcenter'})
header_format = workbook.add_format({'bold': 1, 'align': 'center', 'valign': 'vcenter'})
header_format.set_text_wrap()
integer_format = workbook.add_format({'num_format': '#,##0', 'align': 'center'})
float_format = workbook.add_format({'num_format': '#,##0.00', 'align': 'center'})
percentage_format = workbook.add_format({'num_format': '0.00%', 'align': 'center'})

results_dir = workdir+'/results_secondary_analysis/'+samplename


worksheet2 = workbook.add_worksheet('Image')
scan_all_clusters = results_dir+'/plots/'+samplename+'_all_clusters.png'
worksheet2.insert_image('B2', scan_all_clusters)


### RNA Clusters
counter = 0
while counter < (len(glob.glob(results_dir+'/diff_exp_RNA/'+samplename+'_RNA_cluster_1vs*_positive.csv')) + 1) :

    RNA_cluster_RNA_diff_pos_file = results_dir+'/diff_exp_RNA/'+samplename+'_RNA_cluster_'+str(counter)+'_positive.csv'
    RNA_cluster_RNA_diff_pos_results = pd.read_csv(RNA_cluster_RNA_diff_pos_file, delimiter=',')
    RNA_cluster_RNA_diff_pos_results.to_excel(output, sheet_name='RNA', startrow=62, startcol=(10 * counter + 0), index=False, header=False)
    RNA_cluster_RNA_diff_neg_file = results_dir+'/diff_exp_RNA/'+samplename+'_RNA_cluster_'+str(counter)+'_negative.csv'
    RNA_cluster_RNA_diff_neg_results = pd.read_csv(RNA_cluster_RNA_diff_neg_file, delimiter=',')
    RNA_cluster_RNA_diff_neg_results.to_excel(output, sheet_name='RNA', startrow=62+len(RNA_cluster_RNA_diff_pos_results)+5, startcol=(10 * counter + 0), index=False, header=False)

    counter += 1


worksheet = output.sheets['RNA']

final_data_file = results_dir+'/tables/'+samplename+'_remaining_final.csv'
final_data = pd.read_csv(final_data_file, delimiter=',', index_col=0)
worksheet.write('E2', 'No. of Spots', header_format)
worksheet.write('E3', len(final_data), center)
worksheet.write('E5', 'No. of Clusters', header_format)
SNN_RNA_column = [ x for x in list(final_data.columns) if re.compile('^SNN_SCT').search(x) ][0]
worksheet.write('E6', len(final_data[SNN_RNA_column].unique()), center)

TSNE_RNA = results_dir+'/plots/'+samplename+'_TSNE_RNA_RNA_clusters.png'
worksheet.insert_image('A1', TSNE_RNA)
UMAP_RNA = results_dir+'/plots/'+samplename+'_UMAP_RNA_RNA_clusters.png'
worksheet.insert_image('F1', UMAP_RNA)
features_RNA = results_dir+'/plots/'+samplename+'_feature_RNA_clusters.png'
worksheet.insert_image('L1', features_RNA)
mito_RNA = results_dir+'/plots/'+samplename+'_mito_RNA_clusters.png'
worksheet.insert_image('P1', mito_RNA)
ribo_RNA = results_dir+'/plots/'+samplename+'_ribo_RNA_clusters.png'
worksheet.insert_image('V1', ribo_RNA)
POLR_RNA = results_dir+'/plots/'+samplename+'_POLR_RNA_clusters.png'
worksheet.insert_image('Z1', POLR_RNA)
HSP_RNA = results_dir+'/plots/'+samplename+'_HSP_RNA_clusters.png'
worksheet.insert_image('AF1', HSP_RNA)
HIST_RNA = results_dir+'/plots/'+samplename+'_HIST_RNA_clusters.png'
worksheet.insert_image('AJ1', HIST_RNA)
#clonotypes_RNA = results_dir+'/plots/'+samplename+'_unique_clonotypes_RNA_clusters.png'
#worksheet.insert_image('AP1', clonotypes_RNA)

counter = 0
while counter < (len(glob.glob(results_dir+'/diff_exp_RNA/'+samplename+'_RNA_cluster_1vs*_positive.csv')) + 1) :

    RNA_cluster_RNA_diff_pos_file = results_dir+'/diff_exp_RNA/'+samplename+'_RNA_cluster_'+str(counter)+'_positive.csv'
    RNA_df_for_testing_pos = pd.read_csv(RNA_cluster_RNA_diff_pos_file, delimiter=',')
    RNA_cluster_RNA_diff_neg_file = results_dir+'/diff_exp_RNA/'+samplename+'_RNA_cluster_'+str(counter)+'_negative.csv'
    RNA_df_for_testing_neg = pd.read_csv(RNA_cluster_RNA_diff_neg_file, delimiter=',')
    scan_spots = results_dir+'/plots/'+samplename+'_cluster_'+str(counter)+'.png'
    worksheet.insert_image(xl_rowcol_to_cell(27, 10 * counter + 3), scan_spots)


    cluster_file = results_dir+'/clustering_SCT/'+samplename+'_SCT_cluster_'+str(counter)+'.csv'
    cluster_info = pd.read_csv(cluster_file, delimiter=',')
    worksheet.write(xl_rowcol_to_cell(27, 10 * counter + 0), 'Cluster '+str(counter), bold_format)
    worksheet.write(xl_rowcol_to_cell(28, 10 * counter + 0), 'No. of Spots')
    worksheet.write_number(xl_rowcol_to_cell(28, 10 * counter + 1), len(cluster_info), integer_format)
    worksheet.write(xl_rowcol_to_cell(31, 10 * counter + 0), 'Gene Count Mean')
    worksheet.write_number(xl_rowcol_to_cell(31, 10 * counter + 1), cluster_info.mean()['nCount_Spatial'], integer_format)
    worksheet.write(xl_rowcol_to_cell(32, 10 * counter + 0), 'Gene Count SD')
    worksheet.write_number(xl_rowcol_to_cell(32, 10 * counter + 1), cluster_info.std()['nCount_Spatial'], integer_format)
    worksheet.write(xl_rowcol_to_cell(34, 10 * counter + 0), 'Unique Features Mean')
    worksheet.write_number(xl_rowcol_to_cell(34, 10 * counter + 1), cluster_info.mean()['nFeature_Spatial'], integer_format)
    worksheet.write(xl_rowcol_to_cell(35, 10 * counter + 0), 'Unique Features SD')
    worksheet.write_number(xl_rowcol_to_cell(35, 10 * counter + 1), cluster_info.std()['nFeature_Spatial'], integer_format)
    worksheet.write(xl_rowcol_to_cell(37, 10 * counter + 0), '% MT Mean')
    worksheet.write_number(xl_rowcol_to_cell(37, 10 * counter + 1), cluster_info.mean()['percent.mt'], float_format)
    worksheet.write(xl_rowcol_to_cell(38, 10 * counter + 0), '% MT SD')
    worksheet.write_number(xl_rowcol_to_cell(38, 10 * counter + 1), cluster_info.std()['percent.mt'], float_format)

    if len(cluster_info) > 1 :

        worksheet.write(xl_rowcol_to_cell(29, 10 * counter + 0), 'Most likely marker')
        worksheet.write_formula(xl_rowcol_to_cell(29, 10 * counter + 1), '='+xl_rowcol_to_cell(62, 10 * counter + 0), center)

    worksheet.write(xl_rowcol_to_cell(58, 10 * counter + 0), 'Differential Expression (against other clusters)', bold_format)
    worksheet.write(xl_rowcol_to_cell(59, 10 * counter + 0), '(avg_logFC > 0)')
    worksheet.write(xl_rowcol_to_cell(60, 10 * counter + 0), 'Genes', bold_format)
    worksheet.write(xl_rowcol_to_cell(61, 10 * counter + 0), 'DE Gene', bold_format)
    worksheet.write(xl_rowcol_to_cell(61, 10 * counter + 1), 'p-Value', bold_format)
    worksheet.write(xl_rowcol_to_cell(61, 10 * counter + 2), 'Average Log FC', bold_format)
    worksheet.write(xl_rowcol_to_cell(61, 10 * counter + 3), 'Percentage in Cluster', bold_format)
    worksheet.write(xl_rowcol_to_cell(61, 10 * counter + 4), 'Percentage in Rest of Spots', bold_format)
    worksheet.write(xl_rowcol_to_cell(61, 10 * counter + 5), 'Adjusted p-Value', bold_format)
    buffer1 = len(RNA_df_for_testing_pos)+5
    worksheet.write(xl_rowcol_to_cell(59+buffer1, 10 * counter + 0), '(avg_logFC < 0)', bold_format)
    worksheet.write(xl_rowcol_to_cell(60+buffer1, 10 * counter + 0), 'Genes', bold_format)
    worksheet.write(xl_rowcol_to_cell(61+buffer1, 10 * counter + 0), 'DE Gene', bold_format)
    worksheet.write(xl_rowcol_to_cell(61+buffer1, 10 * counter + 1), 'p-Value', bold_format)
    worksheet.write(xl_rowcol_to_cell(61+buffer1, 10 * counter + 2), 'Average Log FC', bold_format)
    worksheet.write(xl_rowcol_to_cell(61+buffer1, 10 * counter + 3), 'Percentage in Cluster', bold_format)
    worksheet.write(xl_rowcol_to_cell(61+buffer1, 10 * counter + 4), 'Percentage in Rest of Spots', bold_format)
    worksheet.write(xl_rowcol_to_cell(61+buffer1, 10 * counter + 5), 'Adjusted p-Value', bold_format)

    worksheet.set_column(10 * counter + 0, 10 * counter + 0, 21)
    worksheet.set_column(10 * counter + 1, 10 * counter + 1, 12)
    worksheet.set_column(10 * counter + 2, 10 * counter + 2, 14)
    worksheet.set_column(10 * counter + 3, 10 * counter + 3, 19.5)
    worksheet.set_column(10 * counter + 4, 10 * counter + 4, 24.5)
    worksheet.set_column(10 * counter + 5, 10 * counter + 5, 16)

    counter += 1


counter = 0
to_count = len(glob.glob(results_dir+'/diff_exp_RNA/'+samplename+'_RNA_cluster_1vs*_positive.csv')) + 1
while counter < to_count :

    counter2 = 0
    while counter2 < to_count :

        if counter2 != counter :

            RNA_cluster_RNA_diff_pos_file = results_dir+'/diff_exp_RNA/'+samplename+'_RNA_cluster_'+str(counter)+'vs'+str(counter2)+'_positive.csv'
            RNA_cluster_RNA_diff_pos_results = pd.read_csv(RNA_cluster_RNA_diff_pos_file, delimiter=',')
            RNA_cluster_RNA_diff_pos_results.to_excel(output, sheet_name='RNA_'+str(counter), startrow=4, startcol=(10 * counter2 + 0), index=False, header=False)
            RNA_cluster_RNA_diff_neg_file = results_dir+'/diff_exp_RNA/'+samplename+'_RNA_cluster_'+str(counter)+'vs'+str(counter2)+'_negative.csv'
            RNA_cluster_RNA_diff_neg_results = pd.read_csv(RNA_cluster_RNA_diff_neg_file, delimiter=',')
            RNA_cluster_RNA_diff_neg_results.to_excel(output, sheet_name='RNA_'+str(counter), startrow=4+len(RNA_cluster_RNA_diff_pos_results)+5, startcol=(10 * counter2 + 0), index=False, header=False)

            worksheet = output.sheets['RNA_'+str(counter)]
           
            worksheet.write(xl_rowcol_to_cell(0, 10 * counter2 + 0), 'DE against cluster '+str(counter2), bold_format)
            worksheet.write(xl_rowcol_to_cell(1, 10 * counter2 + 0), '(avg_logFC > 0)')
            worksheet.write(xl_rowcol_to_cell(2, 10 * counter2 + 0), 'Genes', bold_format)
            worksheet.write(xl_rowcol_to_cell(3, 10 * counter2 + 0), 'DE Gene', bold_format)
            worksheet.write(xl_rowcol_to_cell(3, 10 * counter2 + 1), 'p-Value', bold_format)
            worksheet.write(xl_rowcol_to_cell(3, 10 * counter2 + 2), 'Average Log FC', bold_format)
            worksheet.write(xl_rowcol_to_cell(3, 10 * counter2 + 3), 'Percentage in Cluster', bold_format)
            worksheet.write(xl_rowcol_to_cell(3, 10 * counter2 + 4), 'Percentage in Rest of Spots', bold_format)
            worksheet.write(xl_rowcol_to_cell(3, 10 * counter2 + 5), 'Adjusted p-Value', bold_format)
            buffer1 = len(RNA_cluster_RNA_diff_pos_results)+5
            worksheet.write(xl_rowcol_to_cell(1+buffer1, 10 * counter2 + 0), '(avg_logFC < 0)', bold_format)
            worksheet.write(xl_rowcol_to_cell(2+buffer1, 10 * counter2 + 0), 'Genes', bold_format)
            worksheet.write(xl_rowcol_to_cell(3+buffer1, 10 * counter2 + 0), 'DE Gene', bold_format)
            worksheet.write(xl_rowcol_to_cell(3+buffer1, 10 * counter2 + 1), 'p-Value', bold_format)
            worksheet.write(xl_rowcol_to_cell(3+buffer1, 10 * counter2 + 2), 'Average Log FC', bold_format)
            worksheet.write(xl_rowcol_to_cell(3+buffer1, 10 * counter2 + 3), 'Percentage in Cluster', bold_format)
            worksheet.write(xl_rowcol_to_cell(3+buffer1, 10 * counter2 + 4), 'Percentage in Rest of Spots', bold_format)
            worksheet.write(xl_rowcol_to_cell(3+buffer1, 10 * counter2 + 5), 'Adjusted p-Value', bold_format)

            worksheet.set_column(10 * counter2 + 0, 10 * counter2 + 0, 21)
            worksheet.set_column(10 * counter2 + 1, 10 * counter2 + 1, 12)
            worksheet.set_column(10 * counter2 + 2, 10 * counter2 + 2, 14)
            worksheet.set_column(10 * counter2 + 3, 10 * counter2 + 3, 19.5)
            worksheet.set_column(10 * counter2 + 4, 10 * counter2 + 4, 24.5)
            worksheet.set_column(10 * counter2 + 5, 10 * counter2 + 5, 16)

        counter2 += 1

    counter += 1



workbook.close()
output.save()
