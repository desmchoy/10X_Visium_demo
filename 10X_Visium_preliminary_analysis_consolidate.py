#!/usr/bin/env python3

import os
import re
import argparse
import pandas as pd
import xlsxwriter
from xlsxwriter.utility import xl_rowcol_to_cell
import glob


workdir = os.getcwd()
project = os.path.basename(workdir)

help_text = '10X_Visium_preliminary_analysis_consolidate.py: Consolidating 10X Visium preliminary analysis outputs'

parser = argparse.ArgumentParser(description=help_text)

parser.add_argument(
    dest='samplelist',
    action='store',
    type=str,
    help='Sample list (required)',
    metavar='samplelist'
    )

arguments = parser.parse_args()
samplelist = arguments.samplelist

output = pd.ExcelWriter(project+'_preliminary_analysis.xlsx', engine='xlsxwriter')
workbook = output.book
center = workbook.add_format({'align': 'center'})
bold_format = workbook.add_format({'bold': 1, 'valign': 'vcenter'})
italic_center_format = workbook.add_format({'italic': 1, 'align': 'center', 'valign': 'vcenter'})
header_format = workbook.add_format({'bold': 1, 'align': 'center', 'valign': 'vcenter'})
header_format.set_text_wrap()
integer_format = workbook.add_format({'num_format': '#,##0', 'align': 'center'})
float_format = workbook.add_format({'num_format': '#,##0.00', 'align': 'center'})
percentage_format = workbook.add_format({'num_format': '0.00%', 'align': 'center'})

filelist = open(samplelist, 'r')
samples = filelist.read().splitlines()
samples = [x.strip(' ') for x in samples]
filelist.close()

for samplename in samples :

    if len(samplename) > 30 :
        sheetname = re.sub('CytAssist_11mm_FFPE_Human_', '', samplename)
    else :
        sheetname = samplename


    results_dir = workdir+'/results_preliminary_analysis/'+samplename

    file1 = results_dir+'/tables/'+samplename+'_feature_summary.csv'
    feature_stats = pd.read_csv(file1, delimiter=',')
    file2 = results_dir+'/tables/'+samplename+'_mito_summary.csv'
    mito_stats = pd.read_csv(file2, delimiter=',')
    file3 = results_dir+'/tables/'+samplename+'_IG_summary.csv'
    IG_stats = pd.read_csv(file3, delimiter=',')
    file4 = results_dir+'/tables/'+samplename+'_TR_summary.csv'
    TR_stats = pd.read_csv(file4, delimiter=',')

    file5 = results_dir+'/tables/'+samplename+'_top_100_variable_genes.csv'
    top_genes = pd.read_csv(file5, delimiter=',')
    file6 = results_dir+'/tables/'+samplename+'_top_100_variable_genes_no_immune.csv'
    top_genes_no_VDJ = pd.read_csv(file6, delimiter=',')


    feature_stats.to_excel(output, sheet_name=sheetname, startrow=5, startcol=11, index=False, header=False)
    mito_stats.to_excel(output, sheet_name=sheetname, startrow=43, startcol=11, index=False, header=False)
    IG_stats.to_excel(output, sheet_name=sheetname, startrow=81, startcol=11, index=False, header=False)
    TR_stats.to_excel(output, sheet_name=sheetname, startrow=119, startcol=11, index=False, header=False)
    top_genes.to_excel(output, sheet_name=sheetname, startrow=197, startcol=10, index=False, header=False)
    top_genes_no_VDJ.to_excel(output, sheet_name=sheetname, startrow=197, startcol=23, index=False, header=False)

    worksheet = output.sheets[sheetname]

    worksheet.set_column(0, 0, 22, bold_format)
    worksheet.set_column(10, 10, 12.5, center)
    worksheet.set_column(11, 11, 12.5, center)
    worksheet.set_column(12, 12, 15, center)
    worksheet.set_column(14, 14, 30, center)
    worksheet.set_column(23, 23, 12.5, center)

    worksheet.write('A1', 'Total Number of Spots')
    worksheet.write('A4', 'Feature Counts')
    worksheet.write('A42', 'Mitochondrial Content')
    worksheet.write('A80', 'IG-expressing Spots')
    worksheet.write('A118', 'TR-expressing Spots')
    worksheet.write('A156', 'Selected Gene Counts')
    worksheet.write('A197', 'Top 100 Variable Genes')
    worksheet.write('O197', 'Top 100 Variable Genes (no VDJ)', bold_format)

    worksheet.write('M5', 'Number of Spots', bold_format)
    worksheet.write('M43', 'Number of Spots', bold_format)
    worksheet.write('M81', 'Number of Spots', bold_format)
    worksheet.write('M119', 'Number of Spots', bold_format)

    worksheet.write_formula('B1', '=SUM(M6:M13)')

    worksheet.insert_image('A5', results_dir+'/plots/'+samplename+'_feature.png')
    worksheet.insert_image('O5', results_dir+'/plots/'+samplename+'_violin_feature.png')
    worksheet.insert_image('A43', results_dir+'/plots/'+samplename+'_mito.png')
    worksheet.insert_image('O43', results_dir+'/plots/'+samplename+'_violin_mito.png')
    worksheet.insert_image('A81', results_dir+'/plots/'+samplename+'_IG.png')
    worksheet.insert_image('O81', results_dir+'/plots/'+samplename+'_violin_IG.png')
    worksheet.insert_image('A119', results_dir+'/plots/'+samplename+'_TR.png')
    worksheet.insert_image('O119', results_dir+'/plots/'+samplename+'_violin_TR.png')
    worksheet.insert_image('A157', results_dir+'/plots/'+samplename+'_selected_gene_counts.png')
    worksheet.insert_image('A198', results_dir+'/plots/'+samplename+'_variable_genes.png')
    worksheet.insert_image('O198', results_dir+'/plots/'+samplename+'_variable_genes_no_immune.png')


workbook.close()
output.save()
