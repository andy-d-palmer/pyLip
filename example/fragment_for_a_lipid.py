import argparse

frag_list = [('$M+$A', 'Precursor ion'),
             ('$M+$A-$A-CH3', ' loss of adduct and CH3'),
             ('$M+$A-$X', 'Loss of headgroup'),
             ('$M+$A-$A-$X+H', ' Loss of headgroup and adduct from precursor ion +H on oxygen'),
             ('$M+$A-$SN1+OH', 'Loss of sn1 acyl chain as ketene (RCH=C=O)'),
             ('$M+$A-$SN2+OH', 'Loss of sn2 acyl chain as ketene (RCH=C=O)'),
             ('$M+$A-$SN1+OH-C3H6O2', 'Loss of sn1 acyl chain as ketene (RCH=C=O) and glycerol'),
             ('$M+$A-$SN2+OH-C3H6O2', 'Loss of sn2 acyl chain as ketene (RCH=C=O) and glycerol'),
             ('$M+$A-$SN1+OH+H2-$X', 'Loss of sn1 acyl chain as ketene (RCH=C=O) and headgroup'),
             ('$M+$A-$SN2+OH+H2-$X', 'Loss of sn1 acyl chain as ketene (RCH=C=O) and headgroup'),
             (
             '$M+$A-$A-$SN1+OH-CH3', 'Loss of sn1 acyl chain as ketene (RCH=C=O), CH3 and chloride from precursor ion'),
             (
             '$M+$A-$A-$SN2+OH-CH3', 'Loss of sn2 acyl chain as ketene (RCH=C=O), CH3 and chloride from precursor ion'),
             ('$M+$A-$SN1-H', 'Loss of sn1 acyl chain RCOOH group'),
             ('$M+$A-$SN1-H-$X+H2', 'Loss of sn1 acyl chain RCOOH group and headgroup'),
             ('$M+$A-$SN2-H-$X+H2', 'Loss of sn2 acyl chain as RCOOH group and headgroup'),
             ('$M+$A-$SN1-H-C3H6O2', 'Loss of sn1 acyl chain RCOOH group and glycerol'),
             ('$M+$A-$SN2-H-C3H6O2', 'Loss of sn2 acyl chain as RCOOH group and glycerol'),
             ('$M+$A-$A-$SN1-H-CH3', 'Loss of sn1 acyl chain RCOOH group, CH3 and adduct'),
             ('$M+$A-$A-$SN2-H-CH3', 'Loss of sn2 acyl chain as RCOOH group, CH3 and adduct'),
             ('$M+$A-$A-$SN1-$HGGly-PO4H2-$X', 'SN2 RCOO- ion'),
             ('$M+$A-$A-$SN2-$HGGly-PO4H2-$X', 'SN1 RCOO- ion'),
             ('$M+$A-$A-$SN3-$HGGly-PO4H2-$X', 'SN3 RCOO- ion'),
             ('$M+$A-$A-$SN1-$SN2-$HGGly', 'Head group ion, -adduct'),
             ('$M+$A-$SN1-$SN2-$HGGly', 'Head group ion'),
             ('$M+$A-$SN1-$SN2+H2O2', 'loss of sn1 & sn2 as ketene (Glycerol-3-phosphate ion)'),
             ('$M+$A-$SN1-$SN2-$X+H2O2-H2O', 'loss of sn1 & sn2 as ketene (Glycerol-3-phosphate ion) loss of H2O'),
             ('$M+$A-$A+H-$SN1-$SN2-$X-$HGGly', 'H2PO4- ion (from phosphate)'),
             ('$M+$A-$A+H-$SN1-$SN2-$X-$HGGly-H2O', 'PO3- ion (from phosphate)'),
             ('$M-$SN1-$SN2-C3H6+$A', ''),
             ('$M+$A-$SN1-$SN2+H2O2-H2O', 'loss of sn1 & sn2 as ketene -H2O'),
             ('$M+$A-$SN1-$SN2+H2O2-(H2O)2', 'loss of sn1 & sn2 as ketene -(H2O)2'),
             ('$M+$A-$A-$SN1-$SN2-$HGGly-CH4', 'Phosphocholine with loss of CH3'),
             ('$M+$A-H2O', ''),
             ('$M+$A-H2O-H2O', ''),
             ('$M+$A-CH2O-H2O', ''),
             ('$M+$A-CH3', ''),
             ('$M+$A-(CH3)2', ''),
             ('$M+$A-(CH3)3', ''),
             ('$M+$A-(CH3)4', ''),
             ('$M+$A-$A-CH3', ''),
             ('$M+$A-$A-(CH3)2', ''),
             ('$M+$A-$A-(CH3)3', ''),
             ('$M+$A-$A-(CH3)4', ''),
             ('$X+$A', ''),
             ('$SN1-O-H2O', ''),
             ('$SN2-O-H2O', ''),
             ('$SN3-O-H2O', ''),
             ('$SN1-O', ''),
             ('$SN2-O', ''),
             ('$SN3-O', ''),
             ('$SN1+$A-OH', ''),
             ('$SN2+$A-OH', ''),
             ('$SN3+$A-OH', ''),
             ('$M+$A-$A-$SN3', ''),
             ('$M+$A-$A-C2HO-CH3', ''),
             ('$M+$A-$SN1-$SN2-$HGGly', ''),
             ('$A+$HGGly+HO4+PO4H2', ''),
             ('$M+$A-$X-PO4H', 'NL headgroup, phosphate group'),
             ]
def get_fragments(lipid_str, adduct="+H", charge=1, ):
    from pyLip.lipid import Lipid
    frag_dict = {}
    for ii, tup in enumerate(frag_list):
        frag_dict['LF{:06.0f}'.format(ii + 1)] = tup[0]
    lipid = Lipid(lipid_str)
    lipid.generate_fragments(adduct, charge, frag_dict)
    for frag_id in sorted(lipid.fragments[adduct]):
        print frag_id, frag_dict[frag_id], lipid.fragments[adduct][frag_id]


'PC(16:0/18:0)'

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="convert centroids from Bruker sqlite to imzML")
    parser.add_argument('lipid_str', type=str, help="")
    parser.add_argument('adduct', type=str, help="")
    parser.add_argument('charge', type=str, help="")

    args = parser.parse_args()
    get_fragments(args.lipid_str, args.adduct, int(args.charge))

