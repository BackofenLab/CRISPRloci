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
import glob
import multiprocessing as mp

from pathlib import Path

def hmmsearch(log_file_path, output_file_path, cutoff, hmm_file_path, fasta_file_path, n_cpus, hmmsearch_cmd='hmmsearch'):
    with open(log_file_path, 'w') as log_file:
            sp.call([hmmsearch_cmd, '--tblout', output_file_path, '-E', str(cutoff),
                     '--cpu', str(n_cpus), hmm_file_path, fasta_file_path],
                     stdout=log_file, stderr=log_file)

def run_hmmsearch(fasta_file, hmm_dir, output_dir, n_cpus, cutoff, hmmsearch_cmd='hmmsearch', use_mp=False):
    if not os.path.exists(output_dir):
        Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    hmm_models = glob.glob(hmm_dir + '/*.hmm')

    if use_mp:
        params = []

    for hmm in hmm_models:
        output_file_path = os.path.join(output_dir, hmm.rsplit('/', 1)[1].replace('.hmm', '.tab'))
        log_file_path = os.path.join(output_dir, hmm.rsplit('/', 1)[1].replace('.hmm', '.log'))

        if use_mp:
            p = [log_file_path, output_file_path, cutoff, hmm, fasta_file, n_cpus, hmmsearch_cmd]
            params.append(p)
        else:
            hmmsearch(log_file_path, output_file_path, cutoff, hmm, fasta_file, n_cpus, hmmsearch_cmd)
    
    if use_mp:
        chunk = len(params) // n_cpus if len(params) % n_cpus == 0 else len(params) // n_cpus + 1
        pool = mp.Pool(processes=n_cpus)
        pool.starmap(hmmsearch, params, chunksize=chunk)
        pool.close()
        pool.join()
