from itertools import zip_longest, takewhile
from functools import partial
import numpy as np

## In theory bytes are significantly faster but requires consistent use
_complement_trans_bytes = bytes.maketrans(b'acgtACGT', b'TGCATGCA')
_complement_trans_str = str.maketrans('acgtACGT', 'TGCATGCA')

genetic_code = {
  'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
  'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
  'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
  'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
  'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
  'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
  'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
  'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
  'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
  'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
  'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
  'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
  'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
  'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
  'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
  'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
}

genetic_code_wo_stop = {k: v for k, v in genetic_code.items() if v != '*'}

as_ascii = partial(bytes, encoding='ascii')

genetic_code_byte = {as_ascii(k): as_ascii(v) for k, v in genetic_code.items()}


def dna_rev_comp(sequence: str) -> str:
  """
  Returns the uppercase reverse complement of the input DNA
  String version
  """
  return sequence.translate(_complement_trans_str)[::-1]


def dna_rev_comp_byte(sequence: bytes) -> bytes:
  """
  Returns the uppercase reverse complement of the input DNA
  String version
  """
  return sequence.translate(_complement_trans_bytes)[::-1]


def translate_dna2aa(orf: str) -> str:
  """
  Ignorant DNA to protein translation.
  Stops on anything not from the genetic code and does not check for length.
  """
  protein = ''
  # Interestingly enough timing revealed that bad practices (loops and string addition)
  # resulted in the fastest translation (for a 150 base dummy sequence)
  for i in range(0, (len(orf) // 3) * 3, 3):
    try:
      protein += genetic_code[orf[i:i + 3]]
    except KeyError:
      break
  return protein


def translate_dna2aa_byte(orf: bytes) -> bytes:
  """
  Ignorant DNA to protein translation.
  Stops on anything not from the genetic code and does not check for length.
  """
  protein = b''
  # Interestingly enough timing revealed that bad practices (loops and string addition)
  # resulted in the fastest translation (for a 150 base dummy sequence)
  for i in range(0, (len(orf) // 3) * 3, 3):
    try:
      protein += genetic_code_byte[orf[i:i + 3]]
    except KeyError:
      break
  return protein

def convert_phred_byte(qual_bytes):
  """
  Takes bytes containing ascii encoded PHRED scores and returns the values as np.array
  (not the probabilities)
  """
  return np.frombuffer(qual_bytes, dtype=np.uint8) - 33

