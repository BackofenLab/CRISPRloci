"""
    Casboundary
    Copyright (C) 2020 Victor Alexandre Padilha <victorpadilha@usp.br>,
                       Omer Salem Alkhnbashi <alkhanbo@informatik.uni-freiburg.de>,
                       Van Dinh Tran <dinh@informatik.uni-freiburg.de>,
                       Shiraz Ali Shah <shiraz.shah@dbac.dk>,
                       André Carlos Ponce de Leon Ferreira de Carvalho <andre@icmc.usp.br>,
                       Rolf Backofen <backofen@informatik.uni-freiburg.de>
    
    This file is part of Casboundary.
    
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import sys
import glob
import pandas as pd
import matplotlib.pyplot as plt

def draw_protein_arrow(ax, x, y, start, end, width, cas_name, strand, font_size_small):
    if cas_name == "cas1" or cas_name == "cas2" or cas_name == "cas4":
        cascolor = 'gold'
    elif cas_name == "cas6":
        cascolor = 'red'
    elif cas_name == "unknown":
        cascolor = 'grey'
    else:
        cascolor = 'cornflowerblue' 
    if(strand == 1):
        ax.arrow(x, y, width ,0, width=width, head_width=2.8, head_length=0.3, fc=cascolor,ec='black',alpha=0.7)
    else:
        ax.arrow(x + width, y,- width,0,width=width, head_width=2.5, head_length=0.3, fc=cascolor,ec='black', alpha=0.7)
    
    ax.text((2 * x + width)/2, (2 * y)/2,cas_name, color='black', ha='center', va='center', fontsize = font_size_small)
    ax.text(x,y + width - 0.6,str(start), color= 'black',ha='left', fontsize = font_size_small)
    ax.text(x + width,y - width - 0.6,str(end), color= 'black',ha='right', fontsize = font_size_small)

def draw_CRISPR_proteins(filename):
    dataframe = pd.read_csv(filename, header=0, index_col=0)
    n = dataframe.shape[0]
    axis = 3.2*n
    x = 0.7
    y = axis/2
    width = 2
    end_last = 0
    strand_last = 0
    font_size_small = int(600/axis)
    fig = plt.figure(figsize=(40,10))
    fig.suptitle(filename, fontsize=20, fontweight='bold')
    ax = fig.add_subplot(111)
    fig.subplots_adjust(top=0.85)
    ax.set_xlabel('xlabel')
    ax.set_ylabel('ylabel')

    for idx, (_, row) in enumerate(dataframe.iterrows()):
        if idx > 0:
            cas_name = row.iloc[4]
            strand = int(row.iloc[2])
            start = int(row.iloc[0])
            end = int(row.iloc[1])

            if idx == 1:
                end_last = end 
                strand_last = strand
            elif start < end_last:
                x = x + width
            else:
                if (strand == 1 and strand_last == 1) or (strand == -1 and strand_last == -1) :
                    x = x + width + 0.7
                elif strand == -1 and strand_last == 1 :
                    x = x + width + 1
                elif strand == 1 and strand_last == -1 :   
                    x = x + width + 0.4
            
            draw_protein_arrow(ax,x,y,start,end,width,cas_name,strand,font_size_small)
            end_last = end
            strand_last = strand
    
    ax.axis([0,axis, 0, axis])
    plt.axis('off')
    plt.grid(b=None)
    if(n < 20):
        dpi = 100
    elif(n > 20 and n < 40):
        dpi = 200
    else:
        dpi = 500
    
    plt.tight_layout()
    out_filename = filename.replace('.csv', '.png')
    print('Saving', out_filename)
    fig.savefig(out_filename, dpi = dpi)
    plt.close(fig)
