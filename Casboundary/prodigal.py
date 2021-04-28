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

import subprocess as sp
import os
import pandas as pd
import re

from collections import defaultdict
from Bio import SeqIO

def run_prodigal(fasta_file, completeness, output_dir, prodigal_path='prodigal'):
    meta = ' -p meta' if completeness == 'partial' else ''
    fasta_file_preffix = fasta_file.split('/')[-1].rsplit('.', 1)[0]
    output_fasta_file = output_dir + '/' + fasta_file_preffix + '_proteins.fasta'
    log_file = output_dir + '/' + fasta_file_preffix + '.prodigal.log'
    prodigal = f'{prodigal_path} -i {fasta_file} -c -m -g 11 -a {output_fasta_file} -q {meta}'
    
    with open(log_file, 'w') as logf:
        sp.call(prodigal.split(), stdout=logf)
    
    return output_fasta_file

def prodigal_fasta_to_genome_info(prodigal_fasta, output_dir, sequence_type):
    org_id = prodigal_fasta.split('/')[-1].replace('.fasta', '')
    info = defaultdict(list)
    sequences = []

    with open(prodigal_fasta, 'r') as f:
        for s in SeqIO.parse(f, 'fasta'):
            if sequence_type == 'dna':
                description = s.description.split('#')
                prot_id_first_part = description[0].strip()
                prot_id_second_part = description[-1].strip().split(';')[0]
                prot_id = prot_id_first_part + '_' + prot_id_second_part
                start = description[1].strip()
                end = description[2].strip()
                strand = description[3].strip()
            else:
                prot_id = s.description.split()[0]
                start = 'none'
                end = 'none'
                strand = 'none'

            info['ProteinId'].append(prot_id)
            info['Start'].append(start)
            info['End'].append(end)
            info['Strand'].append(strand)

            s.id = prot_id
            s.description = ''
            sequences.append(s)
    
    if sequence_type == 'dna':
        info = pd.DataFrame(info).astype({'Start' : 'int'}).sort_values(by='Start')
    else:
        info = pd.DataFrame(info)
    
    info.set_index('ProteinId', inplace=True)
    out_path = output_dir + '/' + prodigal_fasta.split('/')[-1]

    with open(out_path, 'w') as f:
        for s in sequences:
            f.write(s.format('fasta'))
    
    return info
