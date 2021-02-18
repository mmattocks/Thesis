import os
import numpy as np
import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib import cm
from matplotlib import colors

from PIL import Image
from io import BytesIO

def get_chr_props(genome):
    genome_dict = SeqIO.index(genome, &quot;fasta&quot;)
    length_dict = {}
    prop_dict = {}
    
    total_length = 0
    nonchr_length = 0
    for entry in genome_dict:
        length=len(genome_dict[entry].seq)
        if entry.isdigit():
            length_dict[entry] = length
            total_length += length
        else:
            nonchr_length += length
            total_length += length
    length_dict['nonchr']=nonchr_length
    
    for entry in length_dict:
        prop_dict[entry] = length_dict[entry]/total_length

    return prop_dict

def build_position_dict(df):
    observed_dict={}
    observed_nonchr_positions = 0
    
    for chr, chrgrp_data in df.groupby('chr'):
        length=len(chrgrp_data)
        if chr.isdigit():
            observed_dict[chr]=length
        else:
            observed_nonchr_positions+=length
    observed_dict['nonchr']=observed_nonchr_positions
    
    return observed_dict

def build_signal_dict(df):
    observed_dict={}
    observed_nonchr_signal = 0
    
    for chr, chrgrp_data in df.groupby('chr'):
        signal=chrgrp_data['smt_value'].sum()
        if chr.isdigit():
            observed_dict[chr]=signal
        else:
            observed_nonchr_signal+=signal
    observed_dict['nonchr']=observed_nonchr_signal
    
    return observed_dict

#relevant file locations
reference_genome = flow_variables['reference_genome_input']
sib_positions = os.path.join(input_table_1['PathToDanposResults'][0],'pooled',input_table_2['DanposGroupNames'][0]+'.Fnor.smooth.positions.xls')
rys_positions = os.path.join(input_table_1['PathToDanposResults'][0],'pooled',input_table_2['DanposGroupNames'][1]+'.Fnor.smooth.positions.xls')

#read in dataframes
sib_pos_df = pd.read_table(sib_positions, dtype={'chr':'category', 'start':'uint32', 'end':'uint32', 'smt_pos':'uint32', 'smt_value':'float32', 'fuzziness_score':'float32'})
rys_pos_df = pd.read_table(rys_positions, dtype={'chr':'category', 'start':'uint32', 'end':'uint32', 'smt_pos':'uint32', 'smt_value':'float32', 'fuzziness_score':'float32'})

#get proportional lengths of extrachromosomal material from reference genome
prop_dict=get_chr_props(reference_genome)

#get total positions and signals from sib and rys dfs
sib_total_positions=len(sib_pos_df)
sib_total_signal=sib_pos_df['smt_value'].sum()
rys_total_positions=len(rys_pos_df)
rys_total_signal=rys_pos_df['smt_value'].sum()

#build dict of expected position numbers and signal sums for sibs based on percentage of genome comprised by relevant chr
sib_expected_pos_d={}
sib_expected_sig_d={}
for entry in prop_dict:
    sib_expected_pos_d[entry]=sib_total_positions * prop_dict[entry]
    sib_expected_sig_d[entry]=sib_total_signal * prop_dict[entry]

#build dicts of observed position numbers and signal summit sums by chromosome
sib_posn_d=build_position_dict(sib_pos_df)
sib_sig_d=build_signal_dict(sib_pos_df)
rys_posn_d=build_position_dict(rys_pos_df)
rys_sig_d=build_signal_dict(rys_pos_df)

#build arrays of observed positions, signals, and proportional differences
sib_observed_positions=[]
sib_observed_signal=[]
rys_observed_positions=[]
rys_observed_signal=[]

rys_expected_posn_from_sib=[]
rys_expected_sig_from_sib=[]

sib_posn_dif_from_expected=[]
sib_sig_dif_from_expected=[]
rys_posn_dif_from_sib=[]
rys_sig_dif_from_sib=[]

#for chr1-25
for i in range (1,26):
    #build observed position &amp; signal arrays
    sib_observed_positions.append(sib_posn_d[str(i)])
    sib_observed_signal.append(sib_sig_d[str(i)])
    rys_observed_positions.append(rys_posn_d[str(i)])
    rys_observed_signal.append(rys_sig_d[str(i)])
    
    #build rys expected arrays
    rys_expected_posn_from_sib.append(rys_total_positions*(sib_posn_d[str(i)]/sib_total_positions))
    rys_expected_sig_from_sib.append(rys_total_signal*(sib_sig_d[str(i)]/sib_total_signal))
    
    #build proportional difference arrays
    sib_posn_dif_from_expected.append((sib_posn_d[str(i)]-sib_expected_pos_d[str(i)])/sib_expected_pos_d[str(i)])
    sib_sig_dif_from_expected.append((sib_sig_d[str(i)]-sib_expected_sig_d[str(i)])/sib_expected_sig_d[str(i)])
    rys_posn_dif_from_sib.append((rys_posn_d[str(i)]-rys_expected_posn_from_sib[i-1])/rys_expected_posn_from_sib[i-1])
    rys_sig_dif_from_sib.append((rys_sig_d[str(i)]-rys_expected_sig_from_sib[i-1])/rys_expected_sig_from_sib[i-1])

#lastly append nonchrosomal material to all arrays
sib_observed_positions.append(sib_posn_d['nonchr'])
sib_observed_signal.append(sib_sig_d['nonchr'])
rys_observed_positions.append(rys_posn_d['nonchr'])
rys_observed_signal.append(rys_sig_d['nonchr'])

rys_expected_posn_from_sib.append(rys_total_positions*(sib_posn_d['nonchr']/sib_total_positions))
rys_expected_sig_from_sib.append(rys_total_signal*(sib_sig_d['nonchr']/sib_total_signal))

sib_posn_dif_from_expected.append((sib_posn_d['nonchr']-sib_expected_pos_d['nonchr'])/sib_expected_pos_d['nonchr'])
sib_sig_dif_from_expected.append((sib_sig_d['nonchr']-sib_expected_sig_d['nonchr'])/sib_expected_sig_d['nonchr'])
rys_posn_dif_from_sib.append((rys_posn_d['nonchr']-rys_expected_posn_from_sib[25])/rys_expected_posn_from_sib[25])
rys_sig_dif_from_sib.append((rys_sig_d['nonchr']-rys_expected_sig_from_sib[25])/rys_expected_sig_from_sib[25])

#!!sib nonchromosomal signal is blown out relative to other scaffolds. store value and annotate, color separately
NC_sig_dif = sib_sig_dif_from_expected[25]
sib_sig_dif_from_expected[25] = 0

#build proportional difference color arrays
sib_norm = colors.Normalize(vmin=-abs(max(sib_posn_dif_from_expected+sib_sig_dif_from_expected, key=abs)), vmax=abs(max(sib_posn_dif_from_expected+sib_sig_dif_from_expected, key=abs)))
rys_norm = colors.Normalize(vmin=-abs(max(rys_posn_dif_from_sib+rys_sig_dif_from_sib, key=abs)), vmax=abs(max(rys_posn_dif_from_sib+rys_sig_dif_from_sib, key=abs)))

sib_pos_colors = cm.PiYG(sib_norm(sib_posn_dif_from_expected))
sib_sig_colors = cm.PiYG(sib_norm(sib_sig_dif_from_expected))
sib_sig_colors[25] = [0, .1, 1, 1]
rys_pos_colors = cm.PRGn(rys_norm(rys_posn_dif_from_sib))
rys_sig_colors = cm.PRGn(rys_norm(rys_sig_dif_from_sib))

#build figure
fig = plt.figure(figsize=(6,6.5))

gs1 = GridSpec (2,2, wspace=.1, hspace=.05, left=.05, right=.95, top=.95, bottom=.15)
pos_sib = plt.subplot(gs1[0, 0])
pos_rys = plt.subplot(gs1[0, 1])
sig_sib = plt.subplot(gs1[1, 0])
sig_rys = plt.subplot(gs1[1, 1])

pos_sib.set_ylabel('Positions', fontsize=12, weight='bold')
sig_sib.set_ylabel('Occupancy', fontsize=12, weight='bold')

pos_sib.text(.01,.95, 'A', transform=pos_sib.transAxes, fontweight='bold', va='top')
pos_rys.text(.01,.95, 'B', transform=pos_rys.transAxes, fontweight='bold', va='top')
sig_sib.text(.01,.95, 'C', transform=sig_sib.transAxes, fontweight='bold', va='top')
sig_rys.text(.01,.95, 'D', transform=sig_rys.transAxes, fontweight='bold', va='top')

gs2 = GridSpec(4, 2, wspace=.1, hspace=.05, left=.05, right=.95, top=.15, bottom=.05)
sib_cbax = plt.subplot(gs2[3, 0])
rys_cbax = plt.subplot(gs2[3, 1])

#colorbar stuff
ssm = cm.ScalarMappable(cmap=cm.PiYG, norm=sib_norm)
ssm._A =[]
sib_cb=plt.colorbar(ssm, cax=sib_cbax, orientation='horizontal', extend='max', ticks=[-abs(max(sib_posn_dif_from_expected+sib_sig_dif_from_expected, key=abs)), -.3, -.2, -.1, 0, .1, .2, .3])
sib_cb.ax.set_xticklabels(['{0:.2f}'.format(-abs(max(sib_posn_dif_from_expected+sib_sig_dif_from_expected, key=abs))),'-0.3','-0.2','-0.1','0.0','0.1','0.2','0.3'])
sib_cbax.set_xlabel('Fold difference from per-kb uniformity')
sib_cbax.xaxis.set_label_position('top')

rsm = cm.ScalarMappable(cmap=cm.PRGn, norm=rys_norm)
rsm._A =[]
rsm_cb=plt.colorbar(rsm, cax=rys_cbax, orientation='horizontal', ticks=[-abs(max(rys_posn_dif_from_sib+rys_sig_dif_from_sib, key=abs)), -.05, -.025, 0, .025, .05, abs(max(rys_posn_dif_from_sib+rys_sig_dif_from_sib, key=abs))])
rsm_cb.ax.set_xticklabels(['{0:.2f}'.format(-abs(max(rys_posn_dif_from_sib+rys_sig_dif_from_sib, key=abs))),'-0.05','-0.025','0.0','0.025','0.05','{0:.2f}'.format(abs(max(rys_posn_dif_from_sib+rys_sig_dif_from_sib, key=abs)))])
rys_cbax.set_xlabel('Fold difference from sib proportions')
rys_cbax.xaxis.set_label_position('top')

#labels for chr wedges
labels=list(range(1,26))
labels.append('NC')

#make the donut plots
pos_sib.pie(sib_observed_positions, colors=sib_pos_colors, radius=1, startangle=90, counterclock=False, wedgeprops=dict(width=.7,edgecolor='w'), labels=labels, textprops=dict(fontsize=6))
pos_sib.text(0.5, 0.5, 'sib', horizontalalignment='center', verticalalignment='center', transform=pos_sib.transAxes, weight='bold', fontsize='16')
pos_rys.pie(rys_observed_positions, colors=rys_pos_colors, radius=1, startangle=90, counterclock=False, wedgeprops=dict(width=.7,edgecolor='w'), labels=labels, textprops=dict(fontsize=6))
pos_rys.text(0.5, 0.5, 'rys', horizontalalignment='center', verticalalignment='center', transform=pos_rys.transAxes, weight='bold', fontsize='16')

wedges, texts = sig_sib.pie(sib_observed_signal, colors=sib_sig_colors, radius=1, startangle=90, counterclock=False, wedgeprops=dict(width=.7,edgecolor='w'), labels=labels, textprops=dict(fontsize=6))
sig_sib.text(0.5, 0.5, 'sib', horizontalalignment='center', verticalalignment='center', transform=sig_sib.transAxes, weight='bold', fontsize='16')
sig_rys.pie(rys_observed_signal, colors=rys_sig_colors, radius=1, startangle=90, counterclock=False, wedgeprops=dict(width=.7,edgecolor='w'), labels=labels, textprops=dict(fontsize=6))
sig_rys.text(0.5, 0.5, 'rys', horizontalalignment='center', verticalalignment='center', transform=sig_rys.transAxes, weight='bold', fontsize='16')

#annotate NC wedge
wedges[25].set_hatch('//')
ang = (wedges[25].theta2 - wedges[25].theta1)/2. + wedges[25].theta1
y = np.sin(np.deg2rad(ang))
x = np.cos(np.deg2rad(ang))
horizontalalignment = {-1: &quot;right&quot;, 1: &quot;left&quot;}[int(np.sign(x))]
connectionstyle = &quot;angle,angleA=0,angleB={}&quot;.format(ang)
kw = dict(xycoords='data', textcoords='data', arrowprops=dict(arrowstyle=&quot;-&quot;),
            zorder=0, va=&quot;center&quot;, fontsize=9)
kw[&quot;arrowprops&quot;].update({&quot;connectionstyle&quot;: connectionstyle})
sig_sib.annotate('{0:.2f}'.format(NC_sig_dif)+'-fold', xy=(x, y), xytext=(.8*np.sign(x), 1.30*y),
    horizontalalignment=horizontalalignment, **kw)

fig.savefig(os.path.join(input_table_1['PathToDanposResults'][0],'proportional_chromosome_occupancy.png'),dpi=600)
png_memory = BytesIO()
fig.savefig(png_memory, format='png', dpi=600)
PILpng = Image.open(png_memory)
PILpng.save(os.path.join(input_table_1['PathToDanposResults'][0],'proportional_chromosome_occupancy.tiff'))
png_memory.close()

# Copy input to output
output_table_1 = input_table_1.copy()
output_table_2 = input_table_2.copy()
