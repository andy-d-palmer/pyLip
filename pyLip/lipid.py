import pylab as pl
import numpy as np
import re
import sys
import os
import fnmatch

from cpyMSpec.legacy_interface import complete_isodist
from pyMSpec.pyisocalc.pyisocalc import parseSumFormula, ParseError, InvalidFormulaError

# Build a lipid Database
class Lipid():
    # This class should take a lipid shorthand (as used by swiss lipids and defined in Liebisch et al) at sub_species level and return
    # species, class, functional information + mass, sum formula, isotope pattern etc.
    def __init__(self, sub_species_name):
        # Controlled Vocabulary/Mass Definitions
        self.lipid_class = {'Glycerolipids': ('GL',  # abbreviation
                                              (('Gly',),),  # compulsory component(s)
                                              (('FA',), ('FA',), ('FA',)),  # side chain options  #Fatty Acyls':'FA'
                                              ('MG', 'DG', 'TG')),  # sub class abbreviaitons
                            'Glycerophospholipids': ('GP',
                                                     (('Gly',),),
                                                     (('FA',), ('FA',), ('PO4', 'X')),
                                                     ('BMP', 'CL', 'PA', 'PC', 'PE', 'PG', 'PGP', 'PI', 'PS')),
                            'Sphingolipids': ('SP',
                                              (('Sys',),),
                                              (('S',), ('FA',), ('PO4', 'X',)),
                                              (
                                              'Cer', 'SM', 'C1P', 'SPH', 'S1P', 'HexCer', 'GlcCer', 'GalCer', 'Hex2Cer',
                                              'LacCer')),
                            }  # 'Sterollipids'        :('ST',) - todo
        self.head_group_defs = {'Gly': 'C3H3O2',
                                'Sys': 'C3NH6',}
        self.head_groups = {}
        self.head_groups['X'] = {'PO4': 'PO4H2',
                                 'BMP': 'C6H8O5P2',
                                 # 'CL':'',
                                 'PA': 'H2',
                                 'PC': 'C5H13N',
                                 'PE': 'C2H7N',
                                 'PG': 'C3H8O2',
                                 'PGP': 'C3H6O5P',
                                 'PI': 'C6O5H12',
                                 'PS': 'C3H7O2N',
                                 'Cer': 'H',
                                 # 'C1P':'',
                                 # 'SPH':'',
                                 # 'S1P':'',
                                 'SM': 'C2H4',
                                 # 'HexCer':'',
                                 # 'GlcCer':'',
                                 # 'GalCer':'',
                                 # 'Hex2Cer':'',
                                 # 'LacCer':''
                                 }
        self.fragments = {}
        self.modifier_lookup = {'': (0, 0, 0), 'O-': (-1, +2, 0), 'dO-': (-2, +4, 0), 'tO-': (-3, +6, 0),
                                "d": (0, 3, 1), "t": (1, 3, 1)}  # (O,H,N)
        self.parse_subspecies_name(sub_species_name)

    # break name into component parts
    def parse_subspecies_name(self, sub_species_name):
        # todo - input checking, error handling
        sub_species_name = re.sub('[()]', ' ', sub_species_name)
        self.sub_class, chains = sub_species_name.split(' ', 1)
        chain_strs = chains.split('/')

        self.lipid_category = self.get_category()
        if self.lipid_category is None:
            raise KeyError(' lipid category not found')

        self.fatty_acids = []
        for chain in chain_strs:
            self.fatty_acids.append(self.parse_fatty_acid(chain))

        self.sum_formula = self.generate_sum_formula()
        self.mw = self.generate_molecular_weights(self.sum_formula)

    # get carbons and double bonds from chain
    def parse_fatty_acid(self, chain_element):
        # check for modifier
        import re
        ix = re.search("(\d)", chain_element).start()
        modifier = chain_element[0:ix]
        chain_element = chain_element[ix:]
        chain_length, n_dbl_bonds = chain_element.split(':')
        return chain_length, n_dbl_bonds, modifier

    # calculate mass from chain description
    def get_FA_formula(self, chain_length, n_dbl_bonds, modifier):
        n_oxygen = 1 + self.modifier_lookup[modifier][0]
        n_carbon = int(chain_length)
        n_hydrogen = n_carbon * 2 - 2 * int(n_dbl_bonds) + self.modifier_lookup[modifier][1] - 1
        n_nitrogen = 0 + self.modifier_lookup[modifier][2]
        sf = 'C{}H{}O{}N{}'.format(n_carbon, n_hydrogen, n_oxygen, n_nitrogen)
        return sf

    # lookup which category a lipid sub class is in
    def get_category(self):
        for element in self.lipid_class:
            if self.lipid_class[element][3].count(self.sub_class) > 0:
                category = element
        return category

    # construct sum formula for lipid
    def generate_sum_formula(self):
        sum_form = ''
        # compulsory parts
        for part in self.lipid_class[self.lipid_category][1]:
            for subpart in part:
                sum_form += self.head_group_defs[subpart]
        # varying parts
        fa_count = 0
        for part in self.lipid_class[self.lipid_category][2]:
            for subpart in part:
                if subpart == 'FA':
                    sum_form += self.get_FA_formula(self.fatty_acids[fa_count][0], self.fatty_acids[fa_count][1],
                                                    self.fatty_acids[fa_count][2])
                    fa_count += 1
                elif subpart == 'S':
                    #    sum_form += self.get_FA_formula(self.fatty_acids[fa_count][0], self.fatty_acids[fa_count][1], 's'+self.fatty_acids[fa_count][2])
                    #    fa_count+=1
                    sum_form += self.get_FA_formula(self.fatty_acids[fa_count][0], self.fatty_acids[fa_count][1],
                                                    self.fatty_acids[fa_count][2])
                    fa_count += 1
                elif subpart == '':  # nothing to add
                    sum_form = sum_form  # do nothing
                elif subpart == 'X':
                    subpart = self.sub_class
                    sum_form += self.head_groups['X'][subpart]
                else:
                    sum_form += self.head_groups['X'][subpart]

        # collect atomic terms
        self.atoms = self.get_atoms(sum_form)
        # make it into a string
        sum_formula_string = ''
        for a in sorted(self.atoms.keys()):
            sum_formula_string += '{}{}'.format(a, self.atoms[a])
        return sum_formula_string

    def get_atoms(self, sum_form):
        # collate atoms within a string - can appear multiple times, will be combined
        segments = re.findall('[A-Z][a-z]*[0-9]*', sum_form)
        atom = {}
        for segment in segments:
            a = re.findall('[A-Z][a-z]*', segment)[0]
            n = re.findall('[0-9]+', segment)
            if n == []:
                n = 1
            else:
                n = int(n[0])
            if a in atom:
                atom[a] += n
            else:
                atom[a] = n
        return atom

    def generate_molecular_weights(self, sf_string, charge=0):
        atoms = self.parse_sf_string(sf_string)
        sf_string = ''
        for a in sorted(atoms.keys()):
            if atoms[a] > 0:
                sf_string += '{}{}'.format(a, atoms[a])
            elif atoms[a] < 0:
                raise ValueError('negative elements for {} with {}'.format(a, atoms[a]))
        ms_output = complete_isodist(parseSumFormula(sf_string), charge=charge, output='', plot=False, sigma=0.01,
                                     resolution=2000000, cutoff=0.0001, do_centroid=False)
        return ms_output.get_spectrum(source='centroids')[0]

    def add_atoms(self, atoms, tmp_atom):
        for a in tmp_atom:
            if a in atoms:
                atoms[a] += tmp_atom[a]
            else:
                atoms[a] = tmp_atom[a]
        return atoms

    def get_mzs(self, adduct, charge):
        sf_string = self.sum_formula + adduct
        mzs = self.generate_molecular_weights(sf_string, charge=charge)
        return mzs

    def parse_sf_string(self, sf_string):
        # string can be of form A1B2C3-D4E5+F
        from pyMSpec.pyisocalc.pyisocalc import parseSumFormula
        if all([sf_string[0] != '+', sf_string[0] != '-']):
            sf_string = '+' + sf_string
        sf = parseSumFormula(sf_string[1:])
        atoms = {}
        for segment in sf.get_segments():
            atoms[segment.element().name()] = segment.amount()
        return atoms

    def generate_fragments(self, adduct, charge, frag_dict):
        # Probably fragments are:
        # neutral loss of head group components
        # netural loss of FAs (each and both)
        # FAs
        if adduct in self.fragments:
            mz_list = self.fragments[adduct]
        else:
            mz_parent = self.get_mzs(adduct, charge)[0]
            mz_list = {}
            for frag_key in frag_dict:
                frag = frag_dict[frag_key]
                sf = ''
                signs = re.findall('[+-]', frag)
                if not any((frag.startswith('+'), frag.startswith('-'))):
                    signs.insert(0, '+')
                frag = re.split('[+-]', frag)
                # Generate sum formula for residual piece, relative to parent
                try:
                    for f, s in zip(frag, signs):
                        if f.startswith('$'):
                            if f.startswith('$M'):
                                f = self.sum_formula
                            elif f.startswith('$A'):
                                f = adduct
                            elif f.startswith('$SN'):
                                if int(f[-1]) > len(self.fatty_acids):  # can't make this adduct
                                    raise ValueError('fatty acid too long')
                                fa = self.fatty_acids[int(f[-1]) - 1]
                                f = self.get_FA_formula(fa[0], fa[1], fa[2])
                            elif f.startswith('$HG'):
                                f = self.head_group_defs[f[3:]]
                            elif f.startswith('$X'):
                                f = self.head_groups['X'][self.sub_class]
                            else:
                                raise KeyError('{} not recognised'.format(f))
                        sf += '{}{}'.format(s, f)
                    sf = sf.replace("--", "+")
                    sf = sf.replace("+-", "-")
                    sf = sf.replace("-+", "-")
                    sf = sf.replace("++", "+")
                    mz_list[frag_key] = self.generate_molecular_weights(sf, charge=charge)[0]
                    self.fragments[adduct] = mz_list
                except ValueError, e:  # doto option to add verbosity to this
                    if e.message.startswith('negative elements'):
                        continue
                    elif e.message.startswith('fatty acid too long'):
                        continue
                    else:
                        raise
                except ParseError, e:
                    pass
                    # print sf, e
                except InvalidFormulaError, e:
                    pass
                    # print sf, e
                except KeyError, e:
                    if e.message == self.sub_class:  # tried to get an 'X' that isn't there.
                        continue
                    else:
                        raise
        return mz_list

    def get_lipid_classes(self):
        return self.lipid_classes.keys()