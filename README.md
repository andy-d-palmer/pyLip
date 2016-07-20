# README [![Documentation Status](https://readthedocs.org/projects/pyMSpec/badge/?version=latest)](http://pyMSpec.readthedocs.org/en/latest/?badge=latest)#

This repository contains a Python library for generating lipid fragmentation spectra

# Installation

install from github with pip
```
sudo pip install git+https://github.com/alexandrovteam/pyLip
```
# Package Contents

This library is currently a skeleton that includes
* a class for generating fragmentaton spectra from a lipid name
* simple search for experimental data against a database of synthetic fragments

# Example
```python
from pyLip import Lipid
# Instantiate a lipid and produce it's sum formula
lipid = Lipid(PC(18:0/16:0)')
print lipid_.generate_sum_formula()

# Generate some fragments
frag_list = [('$M+$A','Precursor ion'),
      ('$M+$A-$A-CH3',' loss of adduct and CH3'),
      ('$M+$A-$X','Loss of headgroup'),
   }
adduct = "+H"
charge = 1
lipid_.fragments = {}
lipid_.generate_fragments(adduct, charge, frag_dict)

```

