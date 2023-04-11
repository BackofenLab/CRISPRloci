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

import os
import glob
import numpy as np
import pandas as pd
import joblib
import time
import warnings
import shutil
warnings.simplefilter('ignore')

from pathlib import Path
from Bio import SeqIO
from scipy import sparse

# Project imports
from prodigal import run_prodigal, prodigal_fasta_to_genome_info
from hmmer import run_hmmsearch
from utils import extract_targz
from protein_features import get_protein_features
from wrapper import ClassifierWrapper
from draw import draw_CRISPR_proteins

CAS_HMM_TAR_PATH = 'Cas_HMM.tar.gz'
SIG_HMM_TAR_PATH = 'Sig_HMM.tar.gz'
GEN_HMM_TAR_PATH = 'Gen_HMM.tar.gz'
ML_MODELS_TAR_PATH = 'ML_models.tar.gz'

CAS_HMM_DIR = 'Cas_HMM'
SIG_HMM_DIR = 'Sig_HMM'
GEN_HMM_DIR = 'Gen_HMM'
ML_MODELS_DIR = 'ML_models'

GEN_HMM_NAMES_FILE = 'gen_hmm_names.txt'
CAS_HMM_NAMES_FILE = 'cas_hmm_names.txt'

# minimum number of genes to define a cassette
MIN_GENES = 3

def find_potential_regions(fasta_file, sequence_completeness, output_dir, hmmsearch_output_dir, n_cpus, sequence_type, offset=50):
    if sequence_type == 'dna':
        print('Running Prodigal on input Fasta file.')
        proteins_fasta_file = run_prodigal(fasta_file, sequence_completeness, output_dir)
    else:
        fasta_file_cp = output_dir + '/' + fasta_file.split('/')[-1]
        shutil.copyfile(fasta_file, fasta_file_cp)
        proteins_fasta_file = fasta_file_cp

    print('Generating Genome DataFrame')
    info_df = prodigal_fasta_to_genome_info(proteins_fasta_file, output_dir, sequence_type)
    info_df.to_csv(output_dir + '/Genome.csv')

    print('Searching for potential signature proteins')
    sig_output_dir = hmmsearch_output_dir + '/' + SIG_HMM_DIR
    run_hmmsearch(proteins_fasta_file, SIG_HMM_DIR, sig_output_dir, n_cpus, 1, use_mp=True)
    signatures = find_signatures(sig_output_dir)

    print('Extracting potential regions.')
    output_dir += '/potential_regions/'

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    
    regions = []
    for j, sig in enumerate(signatures):
        i = np.where(info_df.index == sig)[0][0]
        start = max(0, i - offset)
        end = min(len(info_df), i + offset) + 1
        r = info_df.iloc[start:end]
        r = r.assign(RefSignature=['no'] * r.shape[0])
        r.at[sig, 'RefSignature'] = 'yes'
        regions.append(r)
        r.to_csv(output_dir + f'/region_{j+1}.csv')
    
    return proteins_fasta_file, regions, info_df

def find_signatures(hmmsearch_output_dir):
    files = glob.glob(hmmsearch_output_dir + '/*.tab')
    signatures = set()

    for f in files:
        hmm_results = np.loadtxt(f, usecols=[0, 2, 4], dtype=np.str)

        if len(hmm_results) and len(hmm_results.shape) == 1:
            hmm_results = np.expand_dims(hmm_results, axis=0)
        
        signatures.update(prot for prot, _, _ in hmm_results)
    
    return sorted(signatures, key=lambda s : (len(s), s))

def extract_regions_sequences(proteins_fasta_file, regions, output_dir):
    output_dir = output_dir + '/potential_regions'

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    
    with open(proteins_fasta_file, 'r') as file_:
        sequences = SeqIO.to_dict(SeqIO.parse(file_, 'fasta'))
    
    regions_prot_ids = [set(r.index) for r in regions]
    try:
        proteins_ids_regions = set.union(*regions_prot_ids)
    except TypeError:
        proteins_ids_regions = []
    sequences = [sequences[prot_id] for prot_id in sorted(proteins_ids_regions, key=lambda s : (len(s), s))]
    regions_fasta_file = output_dir + f'/regions_genes.fasta'
    
    with open(regions_fasta_file, 'w') as out_file:
        for s in sequences:
            out_file.write(s.format('fasta'))
    
    return regions_fasta_file

def build_hmm_matrix(region, region_name, matrix_output_dir, hmmsearch_output_dirs, hmm_to_index, cas=False):
    prot_to_index = dict(zip(region.index, np.arange(region.shape[0])))
    hmm_results_files = sum((glob.glob(ho_dir + '/*.tab') for ho_dir in hmmsearch_output_dirs), [])
    X_hmm = sparse.dok_matrix((region.shape[0], len(hmm_to_index)), dtype=np.double)

    for hmm_f in hmm_results_files:
        hmm = np.loadtxt(hmm_f, usecols=[0, 2, 5], dtype=np.str)
        
        if len(hmm) and len(hmm.shape) == 1:
            hmm = np.expand_dims(hmm, axis=0)

        for prot, model, bitscore in hmm:
            if prot in prot_to_index:
                i = prot_to_index[prot]

                if cas:
                    name = hmm_f.rsplit('/', 1)[1].replace('.tab', '')
                else:
                    name = model
                
                j = hmm_to_index[name]
                bitscore = float(bitscore)
                X_hmm[i, j] = max(X_hmm[i, j], bitscore)
        
    X_hmm = sparse.csr_matrix(X_hmm)
    sparse.save_npz(matrix_output_dir + '/' + region_name + '.npz', X_hmm, compressed=True)
    
    return X_hmm

def build_protein_properties_matrix(region, regions_fasta, region_name, matrix_output_dir):
    with open(regions_fasta, 'r') as file_:
        sequences = SeqIO.to_dict(SeqIO.parse(file_, 'fasta'))
    
    X_prop = []
    for prot_id in region.index:
        seq_record = sequences[prot_id]
        names, prot_features = get_protein_features(str(seq_record.seq))
        X_prop.append(prot_features)
    
    X_prop = np.array(X_prop)
    np.savetxt(matrix_output_dir + '/' + region_name + '.prop.txt', X_prop, header=' '.join(names))
    
    return X_prop

def find_boundaries(y_pred, max_gap, i_sig):
    last_left = i_sig
    left = i_sig - 1
    gap = 0

    while left >= 0 and gap <= max_gap:
        if y_pred[left] == 1:
            gap = 0
            last_left = left
        else:
            gap += 1
        
        left -= 1
    
    last_right = i_sig
    right = i_sig + 1
    gap = 0

    while right < len(y_pred) and gap <= max_gap:
        if y_pred[right] == 1:
            gap = 0
            last_right = right
        else:
            gap += 1
        
        right += 1
    
    return last_left, last_right

def filter_and_merge_cassettes(cassettes_list, info_df, max_gap, protein_sequences, output_dir, min_genes=3, max_unk=0.5, max_overlap=0.0001):
    sorted_indices = sorted(range(len(cassettes_list)), key=lambda i : -cassettes_list[i].shape[0])
    pre_filtered_cassettes_list = []

    for i in sorted_indices:
        c = cassettes_list[i]
        
        if (c['CasType'] == 'unknown').sum() <= int(c.shape[0] * max_unk):
            pre_filtered_cassettes_list.append(c)
    
    cassettes_list = pre_filtered_cassettes_list

    output_dir += '/raw_cassettes/'
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    for i, cassette_pred in enumerate(cassettes_list):
        cassette_sequences = [protein_sequences[idx] for idx in cassette_pred.index]
        cassette_pred.to_csv(output_dir + f'cassette_{i+1}.csv')

        out_path = output_dir + f'cassette_{i+1}.fasta'
        with open(out_path, 'w') as out_file:
            for s in cassette_sequences:
                out_file.write(s.format('fasta'))

    boundaries_list = []
    info_index = list(map(str, list(info_df.index)))
    # print(info_df.index)
    for c in cassettes_list:
        start = np.where(info_df.index.values == c.index[0])[0][0]
        end = np.where(info_df.index.values == c.index[-1])[0][0]
        boundaries_list.append((start, end))
    #     print((c.index[0], c.index[-1]), (start, end))
    #     print(info_df.index[start], info_df.index[end])
    #     print()
    
    # input()

    merged = True
    while merged:
        merged = False
        merged_cassettes = []
        merged_boundaries = []

        while cassettes_list:
            c0 = cassettes_list.pop()
            b0 = boundaries_list.pop()
            indices_to_merge = []

            for i in reversed(range(len(cassettes_list))):
                intersection = np.intersect1d(c0.index.values, cassettes_list[i].index.values)
                intersection = len(intersection) / min(len(c0.index.values), len(cassettes_list[i].index.values))

                if intersection >= max_overlap or \
                   (b0[0] > boundaries_list[i][1] and b0[0] <= boundaries_list[i][1] + max_gap + 1) or \
                   (b0[1] < boundaries_list[i][0] and b0[1] + max_gap + 1 >= boundaries_list[i][0]):
                    indices_to_merge.append(i)
            
            if indices_to_merge:
                for i in indices_to_merge:
                    c1 = cassettes_list.pop(i)
                    c0 = pd.concat([c0, c1]).sort_values(by=['Start', 'RefSignature'], ascending=[True, False])
                    c0 = c0.loc[~c0.index.duplicated(keep='first')]

                    b1 = boundaries_list.pop(i)
                    b0 = min(b0[0], b1[0]), max(b0[1], b1[1])
                
                merged = True
                cassettes_list.append(c0)
                boundaries_list.append(b0)
            
            elif c0.shape[0] >= min_genes and not np.all(c0['CasType'] == 'unknown'):
                merged_cassettes.append(c0)
                merged_boundaries.append(b0)
        
        cassettes_list = merged_cassettes
        boundaries_list = merged_boundaries
        assert len(cassettes_list) == len(boundaries_list)
    
    return [c.sort_values(by='Start') for c in cassettes_list]

def label_best_hmm_hits(region, hmmsearch_output_dirs):
    labels = {idx : ('unknown', 0.0) for idx in region.index}
    
    for hmm_dir in hmmsearch_output_dirs:
        files = glob.glob(f'{hmm_dir}/*.tab')

        for hmm_f in files:
            hmm = np.loadtxt(hmm_f, usecols=[0, 2, 5], dtype=np.str)
            
            if len(hmm) and len(hmm.shape) == 1:
                hmm = np.expand_dims(hmm, axis=0)

            for prot, model, bitscore in hmm:
                if prot in labels:
                    bitscore = float(bitscore)
                    if bitscore > labels[prot][1]:
                        cas = hmm_f.split('/')[-1].split('.')[0].split('_')[0].split('gr')[0].lower()
                        labels[prot] = (cas, bitscore)
    
    return labels

def decompose_into_modules(predictions_df_list):
    for i, predictions in enumerate(predictions_df_list):
        modules = np.repeat('interference', predictions.shape[0]).astype('<U13')
        
        adaptation = np.where((predictions['CasType'] == 'cas1') | \
                              (predictions['CasType'] == 'cas2') | \
                              (predictions['CasType'] == 'cas4'))[0]
        modules[adaptation] = 'adaptation'

        processing = np.where(predictions['CasType'] == 'cas6')[0]
        modules[processing] = 'processing'

        potential_cas = np.where(predictions['CasType'] == 'unknown')[0]
        modules[potential_cas] = 'potential_cas'

        predictions_df_list[i] = predictions.assign(Module=modules)

def write_final_predictions(cassette_pred_list, info_df, protein_sequences, output_dir, max_gap, min_genes=3):
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    
    cassette_pred_list = filter_and_merge_cassettes(cassette_pred_list, info_df, max_gap, protein_sequences, output_dir)

    for i, cassette_pred in enumerate(cassette_pred_list):
        cassette_sequences = [protein_sequences[idx] for idx in cassette_pred.index]
        cassette_pred.to_csv(output_dir + f'/cassette_{i+1}.csv')

        out_path = output_dir + f'/cassette_{i+1}.fasta'
        with open(out_path, 'w') as out_file:
            for s in cassette_sequences:
                out_file.write(s.format('fasta'))

if __name__ == '__main__':
    start_time = time.time()

    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument('-f', '--fasta', dest='fasta_file', help='Organism DNA Fasta file.', metavar='/path/to/organism.fa', type=str)
    parser.add_argument('-c', '--sequence-completeness', nargs='?', dest='sequence_completeness', help='Sequence completeness of DNA Fasta file. Available options: complete or partial.', metavar='seq_comp', choices=['complete', 'partial'], type=str)
    parser.add_argument('-o', '--output-directory', nargs='?', dest='output_dir', help='Output directory path.', metavar='output_dir', default='./casboundary_output', type=str)
    parser.add_argument('-n', '--number-of-cpus', nargs='?', dest='n_cpus', help='Number of CPUs to use (default: 8).', default=8, type=int)
    parser.add_argument('-g', '--maximum-gap', nargs='?', dest='max_gap', help='Maximum number of contiguous gaps allowed in a cassette. Available options: 0 <= gap <= 3 (default: 2).', type=int, choices=range(4), default=2)
    parser.add_argument('-m', '--model', nargs='?', dest='model', help='Which ML model will be used. Available obtions: ert and dnn (default: ert).', choices=['ert', 'dnn', 'ERT', 'DNN'], default='ert')
    parser.add_argument('-ho', '--hmmsearch-output-dir', nargs='?', dest='hmmsearch_output_dir', help='Hmmsearch output directory path.', metavar='hmmsearch_output_dir', default='./hmmsearch_output', type=str)
    parser.add_argument('-d', '--draw-cassettes', dest='draw_cassettes', action='store_true', help='Whether to draw the found cassettes.')
    parser.add_argument('-st', '--sequence-type', nargs='?', dest='sequence_type', default='dna', help='Sequence type. Available options: dna or protein (default: dna).', metavar='seq_type', choices=['dna', 'protein'])
    args = parser.parse_args()

    if not os.path.exists(args.fasta_file):
        raise FileNotFoundError('No such file {}'.format(args.fasta_file))

    if not os.path.exists(CAS_HMM_DIR):
        extract_targz(CAS_HMM_TAR_PATH)
    
    if not os.path.exists(SIG_HMM_DIR):
        extract_targz(SIG_HMM_TAR_PATH)

    if not os.path.exists(GEN_HMM_DIR):
        extract_targz(GEN_HMM_TAR_PATH)
    
    if not os.path.exists(ML_MODELS_DIR):
        extract_targz(ML_MODELS_TAR_PATH)

    if not os.path.exists(args.output_dir):
        Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    
    if not os.path.exists(args.hmmsearch_output_dir):
        Path(args.hmmsearch_output_dir).mkdir(parents=True, exist_ok=True)
    
    proteins_fasta_file, regions_dataframes, info_df = find_potential_regions(args.fasta_file, args.sequence_completeness,
                                                          args.output_dir, args.hmmsearch_output_dir,
                                                          args.n_cpus, args.sequence_type)
    
    regions_fasta_file = extract_regions_sequences(proteins_fasta_file, regions_dataframes, args.output_dir)

    gen_output_dir = args.hmmsearch_output_dir + '/' + GEN_HMM_DIR
    gen_hmm_names = np.loadtxt(GEN_HMM_NAMES_FILE, dtype=np.str)
    gen_hmm_to_index = dict(zip(gen_hmm_names, np.arange(len(gen_hmm_names))))
    
    sig_output_dir = args.hmmsearch_output_dir + '/' + SIG_HMM_DIR
    cas_output_dir = args.hmmsearch_output_dir + '/' + CAS_HMM_DIR
    cas_hmm_names = np.loadtxt(CAS_HMM_NAMES_FILE, dtype=np.str)
    cas_hmm_to_index = dict(zip(cas_hmm_names, np.arange(len(cas_hmm_names))))

    print('Running hmmsearches.')
    run_hmmsearch(regions_fasta_file, GEN_HMM_DIR, gen_output_dir, args.n_cpus, 1000, use_mp=True)
    run_hmmsearch(regions_fasta_file, CAS_HMM_DIR, cas_output_dir, args.n_cpus, 1, use_mp=True)
    
    regions_gen_hmm_list = []
    regions_cas_hmm_list = []
    regions_prop_list = []
    matrix_output_dir = args.output_dir + '/potential_regions'
    
    print('Extracting HMM features and protein properties features.')
    for i, r in enumerate(regions_dataframes):
        r_name = f'region_{i+1}'

        X_gen_hmm = build_hmm_matrix(r, r_name + '.gen', matrix_output_dir, [gen_output_dir], gen_hmm_to_index)
        regions_gen_hmm_list.append(X_gen_hmm)

        X_cas_hmm = build_hmm_matrix(r, r_name + '.cas', matrix_output_dir,
                                     [sig_output_dir, cas_output_dir],
                                     cas_hmm_to_index, cas=True)
        
        regions_cas_hmm_list.append(X_cas_hmm)

        X_prop = build_protein_properties_matrix(r, regions_fasta_file, r_name, matrix_output_dir)
        regions_prop_list.append(X_prop)
    
    assert len(regions_dataframes) == len(regions_gen_hmm_list) == \
           len(regions_cas_hmm_list) == len(regions_prop_list)
    
    print('Finding cassette boundaries.')

    if args.model == 'ert':
        model_scaler = joblib.load(ML_MODELS_DIR + '/ert_bound/MaxAbsScaler_hmm.joblib')
        model_svd = joblib.load(ML_MODELS_DIR + '/ert_bound/TruncatedSVD_hmm.joblib')
        model_clf = joblib.load(ML_MODELS_DIR + '/ert_bound/ExtraTreesClassifier_hmm.joblib')
    else:
        model_scaler = joblib.load(ML_MODELS_DIR + '/dnn_bound/MaxAbsScaler_hmm.joblib')
        model_svd = joblib.load(ML_MODELS_DIR + '/dnn_bound/TruncatedSVD_hmm.joblib')
        model_clf = joblib.load(ML_MODELS_DIR + '/dnn_bound/MLPClassifier_hmm.joblib')

    boundaries_list = []
    regions_dataframes_to_keep = []
    idx_to_remove = set()

    for k, (rdf, X_hmm) in enumerate(zip(regions_dataframes, regions_gen_hmm_list)):
        i_signature = np.where(rdf['RefSignature'] == 'yes')[0]
        assert len(i_signature) == 1
        i_signature = i_signature[0]
        signature_id = rdf.index[i_signature]

        run = True
        for l in range(len(regions_dataframes_to_keep)):
            s, e = boundaries_list[l]
            if signature_id in regions_dataframes_to_keep[l].index[s:e+1]:
                idx_to_remove.add(k)
                run = False
                break

        if run:
            i_signature_rep = [i_signature] * X_hmm.shape[0]
            X_sig = X_hmm[i_signature_rep]
            X_hmm = sparse.hstack((X_sig, X_hmm))

            X_hmm = model_scaler.transform(X_hmm)
            X_hmm = model_svd.transform(X_hmm)
            y_pred = model_clf.predict(X_hmm)
            boundaries = find_boundaries(y_pred, args.max_gap, i_signature)
            boundaries_list.append(boundaries)
            regions_dataframes_to_keep.append(rdf)          

    regions_dataframes = regions_dataframes_to_keep
    regions_cas_hmm_list = [x for k, x in enumerate(regions_cas_hmm_list) if not (k in idx_to_remove)]
    regions_prop_list = [x for k, x in enumerate(regions_prop_list) if not (k in idx_to_remove)]

    assert len(regions_dataframes) == len(regions_cas_hmm_list) == len(regions_prop_list)
    
    print('Labeling Cas proteins.')

    if args.model == 'ert':
        model_maxabs_scaler = joblib.load(ML_MODELS_DIR + '/ert_prot/MaxAbsScaler_hmm2019_prop.joblib')
        model_minmax_scaler = joblib.load(ML_MODELS_DIR + '/ert_prot/MinMaxScaler_hmm2019_prop.joblib')
        model_protein_clf = joblib.load(ML_MODELS_DIR + '/ert_prot/ClassifierWrapper_hmm2019_prop.joblib')
    else:
        model_maxabs_scaler = joblib.load(ML_MODELS_DIR + '/dnn_prot/MaxAbsScaler_hmm2019_prop.joblib')
        model_minmax_scaler = joblib.load(ML_MODELS_DIR + '/dnn_prot/MinMaxScaler_hmm2019_prop.joblib')
        model_protein_clf = joblib.load(ML_MODELS_DIR + '/dnn_prot/ClassifierWrapper_hmm2019_prop.joblib')

    final_predictions = []

    for i, (reg_df, X_hmm, X_prop) in enumerate(zip(regions_dataframes, regions_cas_hmm_list, regions_prop_list)):
        b_start, b_end = boundaries_list[i]

        if abs(b_start - b_end) + 1 >= MIN_GENES:
            start = max(0, b_start - 2)
            added_start = b_start - start
            end = min(reg_df.shape[0] - 1, b_end + 2)
            added_end = end - b_end
            
            reg_df = reg_df.iloc[start:end+1]
            hmm_labels = label_best_hmm_hits(reg_df, [cas_output_dir, sig_output_dir])

            X_hmm = regions_cas_hmm_list[i].toarray()
            X_prop = regions_prop_list[i]

            X_hmm = X_hmm[start:end+1]
            X_hmm = model_maxabs_scaler.transform(X_hmm)
            X_prop = X_prop[start:end+1]
            X_prop = model_minmax_scaler.transform(X_prop)
            X = np.hstack((X_hmm, X_prop))
            y_pred = model_protein_clf.predict(X)

            consensus = []
            for k in range(len(reg_df)):
                if k < added_start or k >= len(reg_df) - added_end or (y_pred[k] != 'unknown' and hmm_labels[reg_df.index[k]][0] == 'unknown'):
                    consensus.append(y_pred[k])
                else:
                    consensus.append(hmm_labels[reg_df.index[k]][0])
            
            # refinement steps
            cs = added_start
            s = b_start
            while cs > 0 and consensus[cs - 1] != 'unknown':
                cs -= 1
                s -= 1

            ce = len(consensus) - added_end - 1
            e = b_end
            while ce < len(consensus) - 1 and consensus[ce + 1] != 'unknown':
                ce += 1
                e += 1

            rdf = regions_dataframes[i].iloc[s:e+1]
            consens = consensus[cs:ce+1]
            rdf = rdf.assign(CasType=consensus[cs:ce+1])
            final_predictions.append(rdf)

    print('Merging and saving final predictions.')
    with open(proteins_fasta_file, 'r') as file_:
        all_sequences = SeqIO.to_dict(SeqIO.parse(file_, 'fasta'))
    
    print('Decomposing into modules.')
    decompose_into_modules(final_predictions)
    
    print('Saving predictions.')
    cassettes_dir = args.output_dir + '/predictions'
    write_final_predictions(final_predictions, info_df, all_sequences, cassettes_dir, args.max_gap)

    if args.draw_cassettes:
        print('Drawing cassettes.')
        files = glob.glob(cassettes_dir + '/*.csv')

        for f in files:
            draw_CRISPR_proteins(f)
    
    end_time = time.time()
    print(f'Total runtime in seconds: {end_time - start_time:.2f}')
