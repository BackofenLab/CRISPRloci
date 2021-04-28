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

from Bio.Data import IUPACData
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils import ProtParamData # Dipeptide Instability Weight Values (DIWV)
                                       # and Kyte and Doolittle's hydrophobicity scale.

from collections import OrderedDict
from itertools import product

"""The methods below are adapted from the biopython package.
Here we implement workarounds for ambiguous amino acids
(which are not handled by biopython).
"""

__AMINOACIDS = frozenset('ACDEFGHIKLMNPQRSTVWYJXOZBU')
__UNKNOWN = {'X', 'O'}

def seq_format(seq):
    return ''.join(str(seq).split()).upper()  # Do the minimum formatting

def __calc_ambiguous_values(aa_dict):
    # Rough approximation of weights/values for Ambiguous Amino Acids.
    # Here we define the weight/value as the mean weight/value between
    # the ambiguous meanings as suggested in
    # https://stackoverflow.com/questions/42159712/biopython-amino-acid-sequence-contains-j-and-cant-calculate-the-molecular-we
    # 
    # We are using https://www.samformat.info/IUPAC-ambiguity-codes and
    # https://www.dnabaser.com/articles/IUPAC%20ambiguity%20codes.html
    # to define ambiguities.
    
    return {"B" : (aa_dict['D'] + aa_dict['N']) / 2,
            "Z" : (aa_dict['E'] + aa_dict['Q']) / 2,
            "J" : (aa_dict['I'] + aa_dict['L']) / 2}

def molecular_weight(seq, monoisotopic=False):
    # Adapted from https://github.com/biopython/biopython/blob/master/Bio/SeqUtils/__init__.py
    seq = seq_format(seq)

    if monoisotopic:
        unambiguous_aa_weights = IUPACData.monoisotopic_protein_weights
        water_weight = 18.010565
    else:
        unambiguous_aa_weights = IUPACData.protein_weights
        water_weight = 18.0153
    
    # rough amino acid weight approximation when monoisotopic is True
    # weight['B'] = (weight['D'] (133.037508) + weight['N'] (132.053492)) / 2
    # weight['Z'] = (weight['E'] (147.053158) + weight['Q'] (146.069142)) / 2
    # weight['J'] = (weight['I'] (131.094629) + weight['L'] (131.094629)) / 2
    #
    # rough amino acid weight when monoisotopic is False
    # weight['B'] = (weight['D'] (133.1027) + weight['N'] (132.1179)) / 2
    # weight['Z'] = (weight['E'] (147.1293) + weight['Q'] (146.1445)) / 2
    # weight['J'] = (weight['I'] (131.1729) + weight['L'] (131.1729)) / 2
    ambiguous_aa_weights = __calc_ambiguous_values(unambiguous_aa_weights)
    
    aa_weights = {**unambiguous_aa_weights, **ambiguous_aa_weights}

    for aa in __UNKNOWN:
        seq = seq.replace(aa, '')
    
    return sum(aa_weights[aa] for aa in seq) - (len(seq) - 1) * water_weight

def instability_index(seq):
    # adapted from https://github.com/biopython/biopython/blob/master/Bio/SeqUtils/ProtParam.py
    # original reference: https://academic.oup.com/peds/article/4/2/155/1491271
    seq = seq_format(seq)

    score = 0.0

    for i in range(len(seq) - 1):
        current_aa, next_aa = seq[i:i+2]

        # considering only valid dipeptides and ignoring the invalid ones
        # as suggested in https://www.biostars.org/p/17833/
        if current_aa in ProtParamData.DIWV and next_aa in ProtParamData.DIWV[current_aa]:
            dipeptide_value = ProtParamData.DIWV[current_aa][next_aa]
            score += dipeptide_value

    return (10.0 / len(seq)) * score

def hydrophobicity(seq):
    # adapted from https://github.com/biopython/biopython/blob/master/Bio/SeqUtils/ProtParam.py
    seq = seq_format(seq)

    # rough kd hydrophobicity approximation for ambiguous values
    # kd['B'] = (kd['D'] (-3.5) + kd['N'] (-3.5)) / 2
    # kd['Z'] = (kd['E'] (-3.5) + kd['Q'] (-3.5)) / 2
    # kd['J'] = (kd['I'] ( 4.5) + kd['L'] ( 3.8)) / 2
    ambiguous_aa_kd = __calc_ambiguous_values(ProtParamData.kd)
    kd = {**ProtParamData.kd, **ambiguous_aa_kd}
    return sum(kd[aa] for aa in seq if aa in kd) / len(seq)

def correct(seq):
    if '*' in seq:
        assert seq.index('*') == len(seq) - 1
        return seq.replace('*', '')
    return seq

def get_protein_features(seq):
    seq = correct(seq)
    prot_analysis = ProteinAnalysis(seq)
    prot_weight = molecular_weight(seq)
    pI = prot_analysis.isoelectric_point()
    aa_count = prot_analysis.count_amino_acids()
    neg_charged_residues = aa_count['D'] + aa_count['E']
    pos_charged_residues = aa_count['K'] + aa_count['R']
    extinction_coefficient_1 = aa_count['Y'] * 1490 + aa_count['W'] * 5500
    extinction_coefficient_2 = aa_count['Y'] * 1490 + aa_count['W'] * 5500 + aa_count['C'] * 125
    instability_idx = instability_index(seq)
    gravy = hydrophobicity(seq)
    secondary_structure_fraction = [frac for frac in prot_analysis.secondary_structure_fraction()]

    names = ['length', 'weight', 'pI', 'neg_charged_residues', 'pos_charged_residues',
            'extinction_coeff1', 'extinction_coeff2', 'instability_index', 'gravy',
            'helix', 'turn', 'sheet']    
    
    return names, [len(seq), prot_weight, pI, neg_charged_residues, pos_charged_residues,
            extinction_coefficient_1, extinction_coefficient_2, instability_idx,
            gravy, *secondary_structure_fraction]
