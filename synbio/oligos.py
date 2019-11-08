"""Functions for oligos. Tm calc."""

import math
from typing import Dict, List, Tuple

from Bio.Seq import Seq


DNA_NN = {
    "init": (0.2, -5.7),
    "init_G/C": (0.0, 0.0),
    "init_A/T": (2.2, 6.9),
    "sym": (0, -1.4),
    "AA/TT": (-7.6, -21.3),
    "AT/TA": (-7.2, -20.4),
    "TA/AT": (-7.2, -20.4),
    "CA/GT": (-8.5, -22.7),
    "GT/CA": (-8.4, -22.4),
    "CT/GA": (-7.8, -21.0),
    "GA/CT": (-8.2, -22.2),
    "CG/GC": (-10.6, -27.2),
    "GC/CG": (-9.8, -24.4),
    "GG/CC": (-8.0, -19.0),
}
"""
SantaLucia & Hicks (2004), Annu. Rev. Biophys. Biomol. Struct 33: 415-440
"""

DNA_INTERNAL_MM = {
    "AG/TT": (1.0, 0.9),
    "AT/TG": (-2.5, -8.3),
    "CG/GT": (-4.1, -11.7),
    "CT/GG": (-2.8, -8.0),
    "GG/CT": (3.3, 10.4),
    "GG/TT": (5.8, 16.3),
    "GT/CG": (-4.4, -12.3),
    "GT/TG": (4.1, 9.5),
    "TG/AT": (-0.1, -1.7),
    "TG/GT": (-1.4, -6.2),
    "TT/AG": (-1.3, -5.3),
    "AA/TG": (-0.6, -2.3),
    "AG/TA": (-0.7, -2.3),
    "CA/GG": (-0.7, -2.3),
    "CG/GA": (-4.0, -13.2),
    "GA/CG": (-0.6, -1.0),
    "GG/CA": (0.5, 3.2),
    "TA/AG": (0.7, 0.7),
    "TG/AA": (3.0, 7.4),
    "AC/TT": (0.7, 0.2),
    "AT/TC": (-1.2, -6.2),
    "CC/GT": (-0.8, -4.5),
    "CT/GC": (-1.5, -6.1),
    "GC/CT": (2.3, 5.4),
    "GT/CC": (5.2, 13.5),
    "TC/AT": (1.2, 0.7),
    "TT/AC": (1.0, 0.7),
    "AA/TC": (2.3, 4.6),
    "AC/TA": (5.3, 14.6),
    "CA/GC": (1.9, 3.7),
    "CC/GA": (0.6, -0.6),
    "GA/CC": (5.2, 14.2),
    "GC/CA": (-0.7, -3.8),
    "TA/AC": (3.4, 8.0),
    "TC/AA": (7.6, 20.2),
    "AA/TA": (1.2, 1.7),
    "CA/GA": (-0.9, -4.2),
    "GA/CA": (-2.9, -9.8),
    "TA/AA": (4.7, 12.9),
    "AC/TC": (0.0, -4.4),
    "CC/GC": (-1.5, -7.2),
    "GC/CC": (3.6, 8.9),
    "TC/AC": (6.1, 16.4),
    "AG/TG": (-3.1, -9.5),
    "CG/GG": (-4.9, -15.3),
    "GG/CG": (-6.0, -15.8),
    "TG/AG": (1.6, 3.6),
    "AT/TT": (-2.7, -10.8),
    "CT/GT": (-5.0, -15.8),
    "GT/CT": (-2.2, -8.4),
    "TT/AT": (0.2, -1.5),
    "AI/TC": (-8.9, -25.5),
    "TI/AC": (-5.9, -17.4),
    "AC/TI": (-8.8, -25.4),
    "TC/AI": (-4.9, -13.9),
    "CI/GC": (-5.4, -13.7),
    "GI/CC": (-6.8, -19.1),
    "CC/GI": (-8.3, -23.8),
    "GC/CI": (-5.0, -12.6),
    "AI/TA": (-8.3, -25.0),
    "TI/AA": (-3.4, -11.2),
    "AA/TI": (-0.7, -2.6),
    "TA/AI": (-1.3, -4.6),
    "CI/GA": (2.6, 8.9),
    "GI/CA": (-7.8, -21.1),
    "CA/GI": (-7.0, -20.0),
    "GA/CI": (-7.6, -20.2),
    "AI/TT": (0.49, -0.7),
    "TI/AT": (-6.5, -22.0),
    "AT/TI": (-5.6, -18.7),
    "TT/AI": (-0.8, -4.3),
    "CI/GT": (-1.0, -2.4),
    "GI/CT": (-3.5, -10.6),
    "CT/GI": (0.1, -1.0),
    "GT/CI": (-4.3, -12.1),
    "AI/TG": (-4.9, -15.8),
    "TI/AG": (-1.9, -8.5),
    "AG/TI": (0.1, -1.8),
    "TG/AI": (1.0, 1.0),
    "CI/GG": (7.1, 21.3),
    "GI/CG": (-1.1, -3.2),
    "CG/GI": (5.8, 16.9),
    "GG/CI": (-7.6, -22.0),
    "AI/TI": (-3.3, -11.9),
    "TI/AI": (0.1, -2.3),
    "CI/GI": (1.3, 3.0),
    "GI/CI": (-0.5, -1.3),
}
"""
Internal mismatch and inosine table (DNA)
Allawi & SantaLucia (1997), Biochemistry 36: 10581-10594
Allawi & SantaLucia (1998), Biochemistry 37: 9435-9444
Allawi & SantaLucia (1998), Biochemistry 37: 2170-2179
Allawi & SantaLucia (1998), Nucl Acids Res 26: 2694-2701
Peyret et al. (1999), Biochemistry 38: 3468-3477
Watkins & SantaLucia (2005), Nucl Acids Res 33: 6258-6267
"""

DNA_TERMINAL_MM = {
    "AA/TA": (-3.1, -7.8),
    "TA/AA": (-2.5, -6.3),
    "CA/GA": (-4.3, -10.7),
    "GA/CA": (-8.0, -22.5),
    "AC/TC": (-0.1, 0.5),
    "TC/AC": (-0.7, -1.3),
    "CC/GC": (-2.1, -5.1),
    "GC/CC": (-3.9, -10.6),
    "AG/TG": (-1.1, -2.1),
    "TG/AG": (-1.1, -2.7),
    "CG/GG": (-3.8, -9.5),
    "GG/CG": (-0.7, -19.2),
    "AT/TT": (-2.4, -6.5),
    "TT/AT": (-3.2, -8.9),
    "CT/GT": (-6.1, -16.9),
    "GT/CT": (-7.4, -21.2),
    "AA/TC": (-1.6, -4.0),
    "AC/TA": (-1.8, -3.8),
    "CA/GC": (-2.6, -5.9),
    "CC/GA": (-2.7, -6.0),
    "GA/CC": (-5.0, -13.8),
    "GC/CA": (-3.2, -7.1),
    "TA/AC": (-2.3, -5.9),
    "TC/AA": (-2.7, -7.0),
    "AC/TT": (-0.9, -1.7),
    "AT/TC": (-2.3, -6.3),
    "CC/GT": (-3.2, -8.0),
    "CT/GC": (-3.9, -10.6),
    "GC/CT": (-4.9, -13.5),
    "GT/CC": (-3.0, -7.8),
    "TC/AT": (-2.5, -6.3),
    "TT/AC": (-0.7, -1.2),
    "AA/TG": (-1.9, -4.4),
    "AG/TA": (-2.5, -5.9),
    "CA/GG": (-3.9, -9.6),
    "CG/GA": (-6.0, -15.5),
    "GA/CG": (-4.3, -11.1),
    "GG/CA": (-4.6, -11.4),
    "TA/AG": (-2.0, -4.7),
    "TG/AA": (-2.4, -5.8),
    "AG/TT": (-3.2, -8.7),
    "AT/TG": (-3.5, -9.4),
    "CG/GT": (-3.8, -9.0),
    "CT/GG": (-6.6, -18.7),
    "GG/CT": (-5.7, -15.9),
    "GT/CG": (-5.9, -16.1),
    "TG/AT": (-3.9, -10.5),
    "TT/AG": (-3.6, -9.8),
}
"""
Terminal mismatch table (DNA)
SantaLucia & Peyret (2001) Patent Application WO 01/94611
"""

DNA_DE = {
    "AA/.T": (0.2, 2.3),
    "AC/.G": (-6.3, -17.1),
    "AG/.C": (-3.7, -10.0),
    "AT/.A": (-2.9, -7.6),
    "CA/.T": (0.6, 3.3),
    "CC/.G": (-4.4, -12.6),
    "CG/.C": (-4.0, -11.9),
    "CT/.A": (-4.1, -13.0),
    "GA/.T": (-1.1, -1.6),
    "GC/.G": (-5.1, -14.0),
    "GG/.C": (-3.9, -10.9),
    "GT/.A": (-4.2, -15.0),
    "TA/.T": (-6.9, -20.0),
    "TC/.G": (-4.0, -10.9),
    "TG/.C": (-4.9, -13.8),
    "TT/.A": (-0.2, -0.5),
    ".A/AT": (-0.7, -0.8),
    ".C/AG": (-2.1, -3.9),
    ".G/AC": (-5.9, -16.5),
    ".T/AA": (-0.5, -1.1),
    ".A/CT": (4.4, 14.9),
    ".C/CG": (-0.2, -0.1),
    ".G/CC": (-2.6, -7.4),
    ".T/CA": (4.7, 14.2),
    ".A/GT": (-1.6, -3.6),
    ".C/GG": (-3.9, -11.2),
    ".G/GC": (-3.2, -10.4),
    ".T/GA": (-4.1, -13.1),
    ".A/TT": (2.9, 10.4),
    ".C/TG": (-4.4, -13.1),
    ".G/TC": (-5.2, -15.0),
    ".T/TA": (-3.8, -12.6),
}
"""DNA dangling ends

Bommarito et al. (2000), Nucl Acids Res 28: 1929-1934
"""

# the energies are the same for each loop stack in the
# reverse complementary direction
DNA_NN.update({k[::-1]: v for k, v in DNA_NN.items()})
DNA_INTERNAL_MM.update({k[::-1]: v for k, v in DNA_INTERNAL_MM.items()})
DNA_TERMINAL_MM.update({k[::-1]: v for k, v in DNA_TERMINAL_MM.items()})
DNA_DE.update({k[::-1]: v for k, v in DNA_DE.items()})

DNA_TRI_TETRA_LOOPS = {
    "AGAAT": (-1.5, 0.0),
    "AGCAT": (-1.5, 0.0),
    "AGGAT": (-1.5, 0.0),
    "AGTAT": (-1.5, 0.0),
    "CGAAG": (-2.0, 0.0),
    "CGCAG": (-2.0, 0.0),
    "CGGAG": (-2.0, 0.0),
    "CGTAG": (-2.0, 0.0),
    "GGAAC": (-2.0, 0.0),
    "GGCAC": (-2.0, 0.0),
    "GGGAC": (-2.0, 0.0),
    "GGTAC": (-2.0, 0.0),
    "TGAAA": (-1.5, 0.0),
    "TGCAA": (-1.5, 0.0),
    "TGGAA": (-1.5, 0.0),
    "TGTAA": (-1.5, 0.0),
    "AAAAAT": (0.5, 0.6),
    "AAAACT": (0.7, -1.6),
    "AAACAT": (1.0, -1.6),
    "ACTTGT": (0.0, -4.2),
    "AGAAAT": (-1.1, -1.6),
    "AGAGAT": (-1.1, -1.6),
    "AGATAT": (-1.5, -1.6),
    "AGCAAT": (-1.6, -1.6),
    "AGCGAT": (-1.1, -1.6),
    "AGCTTT": (0.2, -1.6),
    "AGGAAT": (-1.1, -1.6),
    "AGGGAT": (-1.1, -1.6),
    "AGGGGT": (0.5, -0.6),
    "AGTAAT": (-1.6, -1.6),
    "AGTGAT": (-1.1, -1.6),
    "AGTTCT": (0.8, -1.6),
    "ATTCGT": (-0.2, -1.6),
    "ATTTGT": (0.0, -1.6),
    "ATTTTT": (-0.5, -1.6),
    "CAAAAG": (0.5, 1.3),
    "CAAACG": (0.7, 0.0),
    "CAACAG": (1.0, 0.0),
    "CAACCG": (0.0, 0.0),
    "CCTTGG": (0.0, -2.6),
    "CGAAAG": (-1.1, 0.0),
    "CGAGAG": (-1.1, 0.0),
    "CGATAG": (-1.5, 0.0),
    "CGCAAG": (-1.6, 0.0),
    "CGCGAG": (-1.1, 0.0),
    "CGCTTG": (0.2, 0.0),
    "CGGAAG": (-1.1, 0.0),
    "CGGGAG": (-1.0, 0.0),
    "CGGGGG": (0.5, 1.0),
    "CGTAAG": (-1.6, 0.0),
    "CGTGAG": (-1.1, 0.0),
    "CGTTCG": (0.8, 0.0),
    "CTTCGG": (-0.2, 0.0),
    "CTTTGG": (0.0, 0.0),
    "CTTTTG": (-0.5, 0.0),
    "GAAAAC": (0.5, 3.2),
    "GAAACC": (0.7, 0.0),
    "GAACAC": (1.0, 0.0),
    "GCTTGC": (0.0, -2.6),
    "GGAAAC": (-1.1, 0.0),
    "GGAGAC": (-1.1, 0.0),
    "GGATAC": (-1.6, 0.0),
    "GGCAAC": (-1.6, 0.0),
    "GGCGAC": (-1.1, 0.0),
    "GGCTTC": (0.2, 0.0),
    "GGGAAC": (-1.1, 0.0),
    "GGGGAC": (-1.1, 0.0),
    "GGGGGC": (0.5, 1.0),
    "GGTAAC": (-1.6, 0.0),
    "GGTGAC": (-1.1, 0.0),
    "GGTTCC": (0.8, 0.0),
    "GTTCGC": (-0.2, 0.0),
    "GTTTGC": (0.0, 0.0),
    "GTTTTC": (-0.5, 0.0),
    "GAAAAT": (0.5, 3.2),
    "GAAACT": (1.0, 0.0),
    "GAACAT": (1.0, 0.0),
    "GCTTGT": (0.0, -1.6),
    "GGAAAT": (-1.1, 0.0),
    "GGAGAT": (-1.1, 0.0),
    "GGATAT": (-1.6, 0.0),
    "GGCAAT": (-1.6, 0.0),
    "GGCGAT": (-1.1, 0.0),
    "GGCTTT": (-0.1, 0.0),
    "GGGAAT": (-1.1, 0.0),
    "GGGGAT": (-1.1, 0.0),
    "GGGGGT": (0.5, 1.0),
    "GGTAAT": (-1.6, 0.0),
    "GGTGAT": (-1.1, 0.0),
    "GTATAT": (-0.5, 0.0),
    "GTTCGT": (-0.4, 0.0),
    "GTTTGT": (-0.4, 0.0),
    "GTTTTT": (-0.5, 0.0),
    "TAAAAA": (0.5, -0.3),
    "TAAACA": (0.7, -1.6),
    "TAACAA": (1.0, -1.6),
    "TCTTGA": (0.0, -4.2),
    "TGAAAA": (-1.1, -1.6),
    "TGAGAA": (-1.1, -1.6),
    "TGATAA": (-1.6, -1.6),
    "TGCAAA": (-1.6, -1.6),
    "TGCGAA": (-1.1, -1.6),
    "TGCTTA": (0.2, -1.6),
    "TGGAAA": (-1.1, -1.6),
    "TGGGAA": (-1.1, -1.6),
    "TGGGGA": (0.5, -0.6),
    "TGTAAA": (-1.6, -1.6),
    "TGTGAA": (-1.1, -1.6),
    "TGTTCA": (0.8, -1.6),
    "TTTCGA": (-0.2, -1.6),
    "TTTTGA": (0.0, -1.6),
    "TTTTTA": (-0.5, -1.6),
    "TAAAAG": (0.5, 1.6),
    "TAAACG": (1.0, -1.6),
    "TAACAG": (1.0, -1.6),
    "TCTTGG": (0.0, -3.2),
    "TGAAAG": (-1.0, -1.6),
    "TGAGAG": (-1.0, -1.6),
    "TGATAG": (-1.5, -1.6),
    "TGCAAG": (-1.5, -1.6),
    "TGCGAG": (-1.0, -1.6),
    "TGCTTG": (-0.1, -1.6),
    "TGGAAG": (-1.0, -1.6),
    "TGGGAG": (-1.0, -1.6),
    "TGGGGG": (0.5, -0.6),
    "TGTAAG": (-1.5, -1.6),
    "TGTGAG": (-1.0, -1.6),
    "TTTCGG": (-0.4, -1.6),
    "TTTTAG": (-1.0, -1.6),
    "TTTTGG": (-0.4, -1.6),
    "TTTTTG": (-0.5, -1.6),
}
"""Experimental delta H and delta S for tri/tetra loops

Supplemental Material: Annu.Rev.Biophs.Biomol.Struct.33:415-40
doi: 10.1146/annurev.biophys.32.110601.141800
The Thermodynamics of DNA Structural Motifs
SantaLucia and Hicks, 2004

delta S was computed using delta G and delta H and is in cal / (K x mol)
(versus delta H in kcal / mol)
"""

DNA_MIN_HAIRPIN_LEN = 3
"""Cannot have extremely sharp angles in DNA.
This limit is from Nussinov, et al. (1980)
"""


DNA_INTERNAL_LOOPS = {
    1: (0, 0),
    2: (0, 0),
    3: (0, -10.3),
    4: (0, -11.6),
    5: (0, -12.9),
    6: (0, -14.2),
    7: (0, -14.8),
    8: (0, -15.5),
    9: (0, -15.8),
    10: (0, -15.8),
    11: (0, -16.1),
    12: (0, -16.8),
    13: (0, -16.4),
    14: (0, -17.4),
    15: (0, -17.7),
    16: (0, -18.1),
    17: (0, -18.4),
    18: (0, -18.7),
    19: (0, -18.7),
    20: (0, -19.0),
    21: (0, -19.0),
    22: (0, -19.3),
    23: (0, -19.7),
    24: (0, -20.0),
    25: (0, -20.3),
    26: (0, -20.3),
    27: (0, -20.6),
    28: (0, -21.0),
    29: (0, -21.0),
    30: (0, -21.3),
}
"""Enthalpy and entropy increments for length dependence of internal loops

Were calculated from delta G Table 4 of SantaLucia, 2004:

Annu.Rev.Biophs.Biomol.Struct.33:415-40
doi: 10.1146/annurev.biophys.32.110601.141800
The Thermodynamics of DNA Structural Motifs
SantaLucia and Hicks, 2004

Additional loop sizes are accounted for with the Jacobson-Stockmayer
entry extrapolation formula in paper:
delta G (loop-n) = delta G (loop-x) + 2.44 x R x 310.15 x ln(n / x)

Additional correction is applied for asymmetric loops in paper:
delta G (asymmetry) = |length A - length B| x 0.3 (kcal / mol)
where A and B are lengths of both sides of loop
"""

DNA_BULGE_LOOPS = {
    1: (0, -12.9),
    2: (0, -9.4),
    3: (0, -10.0),
    4: (0, -10.3),
    5: (0, -10.6),
    6: (0, -11.3),
    7: (0, -11.9),
    8: (0, -12.6),
    9: (0, -13.2),
    10: (0, -13.9),
    11: (0, -14.2),
    12: (0, -14.5),
    13: (0, -14.8),
    14: (0, -15.5),
    15: (0, -15.8),
    16: (0, -16.1),
    17: (0, -16.4),
    18: (0, -16.8),
    19: (0, -16.8),
    20: (0, -17.1),
    21: (0, -17.4),
    22: (0, -17.4),
    23: (0, -17.7),
    24: (0, -17.7),
    25: (0, -18.1),
    26: (0, -18.1),
    27: (0, -18.4),
    28: (0, -18.7),
    29: (0, -18.7),
    30: (0, -19.0),
}
"""Enthalpy and entropy increments for length depedence of bulge loops

Were calculated from delta G Table 4 of SantaLucia, 2004:

Annu.Rev.Biophs.Biomol.Struct.33:415-40
doi: 10.1146/annurev.biophys.32.110601.141800
The Thermodynamics of DNA Structural Motifs
SantaLucia and Hicks, 2004

For bulge loops of size 1, the intervening NN energy is used.
Closing AT penalty is applied on both sides
"""

DNA_HAIRPIN_LOOPS = {
    1: (0, 0.0),
    2: (0, 0.0),
    3: (0, -11.3),
    4: (0, -11.3),
    5: (0, -10.6),
    6: (0, -12.9),
    7: (0, -13.5),
    8: (0, -13.9),
    9: (0, -14.5),
    10: (0, -14.8),
    11: (0, -15.5),
    12: (0, -16.1),
    13: (0, -16.1),
    14: (0, -16.4),
    15: (0, -16.8),
    16: (0, -17.1),
    17: (0, -17.4),
    18: (0, -17.7),
    19: (0, -18.1),
    20: (0, -18.4),
    21: (0, -18.7),
    22: (0, -18.7),
    23: (0, -19.0),
    24: (0, -19.3),
    25: (0, -19.7),
    26: (0, -19.7),
    27: (0, -19.7),
    28: (0, -20.0),
    29: (0, -20.0),
    30: (0, -20.3),
}
"""Enthalpy and entropy increments for length depedence of hairpin loops

Were calculated from delta G Table 4 of SantaLucia, 2004:

Annu.Rev.Biophs.Biomol.Struct.33:415-40
doi: 10.1146/annurev.biophys.32.110601.141800
The Thermodynamics of DNA Structural Motifs
SantaLucia and Hicks, 2004

For hairpins of length 3 and 4, the entropy values are looked up
in the DNA_TRI_TETRA_LOOPS Dict

From formula 8-9 of the paper:
An additional 1.6 delta entropy penalty if the hairpin is closed by AT
"""


def calc_tm(seq1: str, seq2: str = "", pcr: bool = True) -> float:
    """Calculate the annealing temperature between seq1 and seq2.

    If seq2 is not provided, it's exact complement is used.
    In otherwords, it's assumed to be an exact match. This tm
    calculate does not account for pseudoknots or anything other
    than an exact, unpadded alignment between seq1 and seq2.

    This is largley influenced by Bio.SeqUtils.MeltingTemp with
    some different defaults. Here, the reaction mixture is assumed to
    be PCR and concentrations for Mg, Tris, K, and dNTPs are included
    that match a typical PCR reaction according to Thermo and NEB. Additionally,
    the salt correction formula from IDT's Owczarzy et al. (2008) is used.

    NEB: https://www.neb.com/tools-and-resources/usage-guidelines/guidelines-for-pcr-optimization-with-taq-dna-polymerase
    ThermoFisher: https://www.thermofisher.com/order/catalog/product/18067017?SID=srch-srp-18067017

    NOTE: Sequences are assumed not to be symmetrical. Oligo not binding to self.

    Args:
        seq1: The seq whose tm is calculated

    Keyword Args:
        seq2: The seq that seq1 anneals to in 3' -> 5' direction
        pcr: Whether tm is being calculated for the oligo is in a
            PCR reaction mixture. If so, ion and Tris concentrations
            that match a typical NEB/ThermoFisher PCR mixture are used

    Returns:
        The estimated tm as a float
    """

    if not seq2:
        seq2 = str(Seq(seq1).complement())

    if len(seq1) != len(seq2):
        raise ValueError(
            f"length mismatch between seq1 {len(seq1)} and seq2 {len(seq2)}"
        )

    # sum enthalpy and entropy. Enthaply is first value of each tuple and
    # entropy is the second value of each tuple in:
    # SantaLucia & Hicks (2004), Annu. Rev. Biophys. Biomol. Struct 33: 415-440

    # start with initiation enthalpy and entropy
    delta_h, delta_s = DNA_NN["init"]

    # add in initial A/T and initial G/Cs
    init = seq1[0] + seq1[-1]
    init_at = init.count("A") + init.count("T")
    init_gc = init.count("G") + init.count("C")
    init_at_h, init_at_s = DNA_NN["init_A/T"]
    init_gc_h, init_gc_s = DNA_NN["init_G/C"]
    delta_h += init_at * init_at_h + init_gc * init_gc_h
    delta_s += init_at * init_at_s + init_gc * init_gc_s

    # work through each nearest neighbor pair
    for i in range(len(seq1) - 1):
        pair1 = seq1[i] + seq1[i + 1]
        pair2 = seq2[i] + seq2[i + 1]
        pair = pair1 + "/" + pair2

        # assuming internal neighbor pair
        pair_delta_h, pair_delta_s = 0.0, 0.0
        if pair in DNA_NN:
            pair_delta_h, pair_delta_s = DNA_NN[pair]
        elif pair in DNA_INTERNAL_MM:
            pair_delta_h, pair_delta_s = DNA_INTERNAL_MM[pair]

        # overwrite if it's a terminal pair
        if i in (0, len(seq1) - 2):
            if pair in DNA_TERMINAL_MM:
                pair_delta_h, pair_delta_s = DNA_TERMINAL_MM[pair]

        delta_h += pair_delta_h
        delta_s += pair_delta_s

    # adjust salt based on mode
    if pcr:
        seq1_conc = 250.0
        seq2_conc = 0.0
        Na = 0
        K = 50
        Tris = 2  # see Thermo
        Mg = 1.5  # see NEB
        dNTPs = 0.2  # see NEB
    else:
        seq1_conc = 25.0
        seq2_conc = 25.0
        Na = 50
        K = 0
        Tris = 0
        Mg = 0
        dNTPs = 0

    # salt correction for deltaS
    # copied-pasted from Bio.SeqUtils' use of a decision tree by:
    # Owczarzy et al. (2008), Biochemistry 4 7: 5336-5353
    Mon = Na + K + Tris / 2.0  # monovalent ions
    mg = Mg * 1e-3  # Lowercase ions (mg, mon, dntps) are molar
    mon = Mon * 1e-3

    # coefficients to a multi-variate from the paper
    a, b, c, d, e, f, g = 3.92, -0.911, 6.26, 1.42, -48.2, 52.5, 8.31

    if dNTPs > 0:
        dntps = dNTPs * 1e-3
        ka = 3e4  # Dissociation constant for Mg:dNTP
        # Free Mg2+ calculation:
        mg = (
            -(ka * dntps - ka * mg + 1.0)
            + math.sqrt((ka * dntps - ka * mg + 1.0) ** 2 + 4.0 * ka * mg)
        ) / (2.0 * ka)
    if Mon > 0:
        R = math.sqrt(mg) / mon
        if R < 0.22:
            corr = (4.29 * _gc(seq1) / 100 - 3.95) * 1e-5 * math.log(
                mon
            ) + 9.40e-6 * math.log(mon) ** 2
            return corr
        elif R < 6.0:
            a = 3.92 * (0.843 - 0.352 * math.sqrt(mon) * math.log(mon))
            d = 1.42 * (1.279 - 4.03e-3 * math.log(mon) - 8.03e-3 * math.log(mon) ** 2)
            g = 8.31 * (0.486 - 0.258 * math.log(mon) + 5.25e-3 * math.log(mon) ** 3)
    corr = (
        a
        + b * math.log(mg)
        + (_gc(seq1) / 100) * (c + d * math.log(mg))
        + (1 / (2.0 * (len(seq1) - 1))) * (e + f * math.log(mg) + g * math.log(mg) ** 2)
    ) * 1e-5

    # tm with concentration consideration
    k = (seq1_conc - (seq2_conc / 2.0)) * 1e-9
    R = 1.9872
    tm = (delta_h * 1000.0) / (delta_s + R * math.log(k)) - 273.15

    # add in salt correction
    tm = 1 / (1 / (tm + 273.15) + corr) - 273.15

    return tm


Cache = Dict[Tuple[int, int], float]
"""A map from i, j tuple to a value."""


def fold(seq: str, temp: float = 32.0) -> float:
    """Fold the DNA sequence and return lowest free energy score.

    Based on the approach described in:
    Nussinov R, Jacobson AB.
    Fast algorithm for predicting the secondary structure of single-stranded RNA.
    Proc Natl Acad Sci U S A. 1980;77(11):6309–6313. doi:10.1073/pnas.77.11.6309

    Args:
        seq: The sequence to fold

    Keyword Args:
        temp: The temperature the fold takes place in, in Celcius

    Returns:
        The lowest free energy fold possible for the sequences.
    """

    seq = seq.upper()
    temp = temp + 273.15  # kelvin

    # inductive fill step
    e_cache: Cache = {}
    k_cache: Cache = {}
    # for increasing fragment length
    for f_len in range(DNA_MIN_HAIRPIN_LEN + 1, len(seq)):
        # for increasing start index
        for i in range(len(seq) - f_len):
            # fill the energy and backtracking cache
            _e(seq, i, i + f_len, temp, e_cache, k_cache)

    # backtracking for structure
    pairs = []

    min_e = 0.0
    min_key = None
    for key, v in e_cache.items():
        if v < min_e:
            min_e = v
            min_key = key
    print(round(min_e, 2), min_key)

    _backtrack(seq, 0, len(seq) - 1, k_cache, e_cache, pairs)

    # print the folding results
    print("\nfolding:", seq)
    print("min free energy structure:", min_e)
    print("pairs:")
    for pair in sorted(pairs, key=lambda x: x[1]):
        print(pair[0] + " " + str(pair[1] + 1) + "-" + str(pair[2] + 1) +" " +  str(pair[3]))

    return min_e


def _e(seq: str, i: int, j: int, temp: float, e_cache: Cache, k_cache: Cache) -> float:
    """Find, store and return the minimum free energy of the structure between i and j

    Based on formula E(i,j) (2) from Nussinov, et al. (1980)

    Args:
        seq: The sequence being folded
        i: The start index
        j: The end index (inclusive)
        temp: The temperature in Kelvin
        e_cache: A free energy cache with key (start, end)
            to value (min energy structure)
        k_cache: A cache for which value of i, j is the best base
            for separating the "upper" and "lower" sections of the RNA (see paper)

    Returns:
        The minimum energy folding structure possible between i and j on seq
    """

    key = (i, j)
    if key in e_cache:
        return e_cache[key]

    # check whether it's a small hairpin with a pre-computed energy
    if j - i <= DNA_MIN_HAIRPIN_LEN + 1:
        d_g = 1000.0
        e_cache[key] = d_g
        k_cache[key] = i
        return d_g

    # check all possible interactions between k and j
    min_e, min_k = _e(seq, i, j - 1, temp, e_cache, k_cache), -1
    for k in range(i, j - (DNA_MIN_HAIRPIN_LEN + 2)):
        # calculate E_kj, which depends on the energies of
        # the inward basepairs. See Fig 5 of paper.
        # Be it a pair (A), single strand bulge on k side (B),
        # single strand bulge on k side (C), internal loop (D),
        # or branched structure (E)
        k_i = k + 1  # k inwards
        j_i = j - 1  # j inwards
        k_1 = k_cache[(k_i, j_i)]  # index of inward k with lowest divide
        while k_1 < 0:
            j_i -= 1
            k_1 = k_cache[(k_i, j_i)]

        pair = seq[k] + seq[k_1] + "/" + seq[j] + seq[j_i]
        bulge_left = k_1 > k + 1 and pair in DNA_NN
        bulge_right = j_i < j - 1 and pair in DNA_NN

        if (
            k_1 == k + 1 and
            j_i == j - 1 and
            (pair in DNA_NN or pair in DNA_INTERNAL_MM)
        ):
            # Fig 5A: it's a neighboring pair in a helix (simplest case)
            e_kj = _pair(pair, seq, k, j, temp)
        elif bulge_left and bulge_right:
            # Fig 5D: it's an internal loop
            loop_left = seq[k : k_1 + 1]
            loop_right = seq[j_i : j + 1]
            e_kj = _internal_loop(loop_left, loop_right, temp)
        elif bulge_left or bulge_right:
            # Fig 5B: it's a bulge on the k side OR
            # Fig 5C: it's a bulge on the j side
            if bulge_left:
                e_kj = _bulge(pair, seq, k, k_1, temp)
            elif bulge_right:
                e_kj = _bulge(pair, seq, j_i, j, temp)
        else:
            e_kj = _hairpin(seq, k, j, temp)

        # sum the energy and compare to alternatives
        e_left = _e(seq, i, k - 1, temp, e_cache, k_cache)
        e_right = _e(seq, k + 1, j - 1, temp, e_cache, k_cache)
        e_total = e_left + e_kj + e_right

        if i == 0 and j == 30:
            # print(i, j)
            pass

        if e_total < min_e:
            min_e = e_total
            min_k = k

    e_cache[key] = min_e
    k_cache[key] = min_k
    return min_e


def _d_g(d_h: float, d_s: float, temp: float) -> float:
    """Find the free energy given delta h, s and temp

    Args:
        d_h: The enthalpy increment in kcal / mol
        d_s: The entropy increment in cal / mol
        temp: The temperature in Kelvin

    Returns:
        The free energy increment in kcal / (mol x K)
    """

    return d_h - temp * (d_s / 1000.0)


def _j_s(n: int, x: int, d_g_x: float, temp: float) -> float:
    """Estimate the free energy of length n based on one of length x.

    The Jacobson-Stockmayer entry extrapolation formula is used
    for bulges, hairpins, etc that fall outside the 30nt upper limit
    for pre-calculated free-energies. See SantaLucia and Hicks (2014).

    Args:
        n: Length of element without known free energy value
        x: Length of element with known free energy value (d_g_x)
        d_g_x: The free energy of the element x
        temp: Temperature in Kelvin

    Returns:
        The free energy for a structure of length n
    """

    R = 1.9872e-3
    return d_g_x + 2.44 * R * temp * math.log(n / float(x))


def _pair(pair: str, seq: str, k: int, j: int, temp: float) -> float:
    """Get the free energy for a pair.

    Using the indexes k and j, check whether it's at the end of
    the sequence or internal. Then check whether it's a match
    or mismatch, and return.

    Args:
        pair: The pair sequence, ex: (AG/TC)
        seq: The full folding sequence
        k: The start index on left size of helix
        j: The end index on right size of helix
        temp: Temperature in Kelvin

    Returns:
        The free energy of the NN pairing
    """

    if k > 0 and j < len(seq) - 1:
        # it's internal
        d_h, d_s = DNA_NN[pair] if pair in DNA_NN else DNA_INTERNAL_MM[pair]
        return _d_g(d_h, d_s, temp)

    # it's terminal
    d_h, d_s = DNA_NN[pair] if pair in DNA_NN else DNA_TERMINAL_MM[pair]
    return _d_g(d_h, d_s, temp)


def _hairpin(seq: str, i: int, j: int, temp: float) -> float:
    """Calculate the free energy of a hairpin.

    Args:
        seq: The sequence we're folding
        i: The index of start of hairpin
        j: The index of end of hairpin
        temp: Temperature in Kelvin

    Returns:
        The free energy increment from the hairpin structure
    """

    if j - i <= DNA_MIN_HAIRPIN_LEN:
        return 1000.0

    hairpin = seq[i : j + 1]
    hairpin_len = len(hairpin) - 2
    pair = hairpin[0] + hairpin[1] + "/" + hairpin[-1] + hairpin[-2]

    if pair not in DNA_TERMINAL_MM or hairpin_len < DNA_MIN_HAIRPIN_LEN:
        # not known terminal pair, nothing to close "hairpin", treat as bulge
        return _bulge(pair, seq, i, j, temp)

    d_g = 0
    if hairpin in DNA_TRI_TETRA_LOOPS:
        # it's a pre-known hairpin with known value
        d_h, d_s = DNA_TRI_TETRA_LOOPS[hairpin]
        d_g = _d_g(d_h, d_s, temp)

    # add penalty based on size
    if hairpin_len in DNA_HAIRPIN_LOOPS:
        d_h, d_s = DNA_HAIRPIN_LOOPS[hairpin_len]
        d_g += _d_g(d_h, d_s, temp)
    else:
        # it's too large, extrapolate
        hairpin_max_len = max(DNA_HAIRPIN_LOOPS.keys())
        d_h, d_s = DNA_HAIRPIN_LOOPS[hairpin_max_len]
        d_g_inc = _d_g(d_h, d_s, temp)
        d_g += _j_s(hairpin_len, hairpin_max_len, d_g_inc, temp)

    # add penalty for a terminal mismatch
    d_h, d_s = DNA_TERMINAL_MM[pair]
    d_g += _d_g(d_h, d_s, temp)

    # add penalty if length 3 and AT closing, formula 8 from SantaLucia, 2004
    if hairpin_len == 3 and (pair.startswith("A") or pair.startswith("T")):
        d_g += 0.5  # TODO: convert to entropy

    return d_g


def _bulge(pair: str, seq: str, i: int, j: int, temp: float) -> float:
    """Calculate the free energy associated with a bulge.

    Args:
        pair: The NN pair outside the bulge
        seq: The full folding DNA sequence
        i: The start index of the bulge
        j: The end index of the bulge
        temp: Temperature in Kelvin

    Returns:
        The increment in free energy from the bulge
    """

    n = j - i - 1
    if n <= 0:
        return 0.0

    # add penalty based on size
    if n in DNA_BULGE_LOOPS:
        d_h, d_s = DNA_BULGE_LOOPS[n]
        d_g = _d_g(d_h, d_s, temp)
    else:
        # it's too large for pre-calculated list, extrapolate
        bulge_max_len = max(DNA_BULGE_LOOPS.keys())
        d_h, d_s = DNA_BULGE_LOOPS[bulge_max_len]
        d_g = _d_g(d_h, d_s, temp)
        d_g = _j_s(n, bulge_max_len, d_g, temp)

    if n == 1 and (pair in DNA_NN or pair in DNA_INTERNAL_MM):
        # if len 1, include the delta G of intervening NN (SantaLucia 2004)
        d_g += _pair(pair, seq, i, j, temp)

    # penalize AT terminal bonds
    if pair.count("A"):
        d_g += 0.5

    return d_g


def _internal_loop(left: str, right: str, temp: float) -> float:
    """Calculate the free energy of an internal loop.

    The first and last bp of both left and right sequences
    are not themselves parts of the loop, but are the terminal
    bp on either side of it. They are needed for when there's
    a single internal looping bp (where just the mismatching
    free energies are used)

    Note that both left and right sequences are in 5' to 3' direction

    Args:
        left: The sequence on the left side
        right: The sequence on the right side
        temp: Temperature in Kelvin

    Returns:
        The free energy associated with the internal loop
    """

    loop_left = len(left) - 2
    loop_right = len(right) - 2

    # apply a penalty based on loop size
    loop_len = loop_left + loop_right
    if loop_len in DNA_INTERNAL_LOOPS:
        d_h, d_s = DNA_INTERNAL_LOOPS[loop_len]
        d_g = _d_g(d_h, d_s, temp)
    else:
        # it's too large an internal loop, extrapolate
        loop_max_len = max(DNA_INTERNAL_LOOPS.keys())
        d_h, d_s = DNA_INTERNAL_LOOPS[loop_max_len]
        d_g = _d_g(d_h, d_s, temp)
        d_g = _j_s(loop_len, loop_max_len, d_g, temp)

    # apply an asymmetry penalty
    loop_asymmetry = abs(loop_left - loop_right)
    d_g += 0.3 * loop_asymmetry

    # apply penalty based on the mismatching pairs on either side the loop
    pair_left = left[:2] + "/" + right[-2::][::-1]
    d_h, d_s = DNA_NN[pair_left] if pair_left in DNA_NN else DNA_TERMINAL_MM[pair_left]
    d_g += _d_g(d_h, d_s, temp)

    pair_right = left[-2:] + "/" + right[:2][::-1]
    d_h, d_s = (
        DNA_NN[pair_right] if pair_right in DNA_NN else DNA_TERMINAL_MM[pair_right]
    )
    d_g += _d_g(d_h, d_s, temp)

    return d_g


def _backtrack(
    seq: str, i: int, j: int, k_cache: Cache, e_cache: Cache, pairs: List[Tuple[int, int]]
):
    """Backtrack and aquire the list of basepairs and total minimum energy

    Args:
        seq: The DNA sequence being folded
        i: The index of the left-side of span
        j: The index of the right-side of span
        k_cache: The low energy track cache
        e_cache: The min energy in a branch cache
        pairs: The list of basepairs in the min-energy structure
    """

    key = (i, j)

    if key not in k_cache:
        return

    k = k_cache[key]

    while k == i or k < 0:
        j -= 1
        key = (i, j)

        if key not in k_cache:
            return

        k = k_cache[key]

    if len(set([i, j, k])) != 3:
        return

    pair = seq[k] + seq[j]
    pairs.append((pair, k, j, round(e_cache[key], 2)))

    _backtrack(seq, i, k - 1, k_cache, e_cache, pairs)
    _backtrack(seq, k, j, k_cache, e_cache, pairs)


def _gc(seq: str) -> float:
    """Return the GC ratio of a sequence."""

    seq = seq.upper()

    return float(seq.count("G") + seq.count("C")) / float(len(seq))
