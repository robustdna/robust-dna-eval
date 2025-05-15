# Robust DNA Evaluation Toolkit

This repository provides tools to generate and evaluate adversarial perturbations on DNA sequences. The attack methods are designed to simulate biologically plausible mutations, and evaluation metrics help quantify their effect on model performance and biological validity.

## Installation

```bash
pip install numpy pandas tqdm
```

Note RNAfold is required for MFE calculation and must be installed separately.

---

## Attack Methods

All attacks take a pandas Series of DNA sequences as input and return a mutated version. Each function supports `mutation_rate` and `iteration` parameters.

### 1. Nucleotide Substitution

Randomly replaces individual nucleotides in the sequence.

```python
from attacks.nucleotide import nucleotide_attack
mutated = nucleotide_attack(sequences, mutation_rate=0.1, iteration=1)
```

---

### 2. Codon Substitution

Replaces 3-mer codons with random codons, disregarding amino acid preservation.

```python
from attacks.codon import codon_attack
mutated = codon_attack(sequences, mutation_rate=0.1, iteration=1)
```

---

### 3. Synonymous Codon Substitution

Mutates codons to biologically synonymous alternatives that code for the same amino acid.

```python
from attacks.synonymous_codon import synonymous_codon_attack
mutated = synonymous_codon_attack(sequences, mutation_rate=0.1, iteration=1)
```

---

### 4. GC-Guided Synonymous Codon Substitution

Synonymous mutation that minimizes GC-content deviation from the original sequence.

```python
from attacks.GC_guide_synonymous_codon import synonymous_codon_attack_gc
mutated = synonymous_codon_attack_gc(sequences, mutation_rate=0.1, lambda_gc=1.0, iteration=1)
```

---

### 5. Backtranslation

Translates DNA to amino acids and randomly re-encodes to DNA with synonymous codons. This preserves the full protein sequence.

```python
from attacks.backtranslation import backtranslation_attack
mutated = backtranslation_attack(sequences, iteration=1)
```

Note Mutation rate is not used for backtranslation. Perturbation is controlled by the number of iterations.

---

## Evaluation Metrics

Three evaluation criteria are supported for comparing original and perturbed sequences.

### 1. Attack Success Rate ASR

```python
import numpy as np
true_labels = [...]
original_preds = [...]
adv_preds = [...]

correct = np.array(original_preds) == np.array(true_labels)
changed = np.array(original_preds) != np.array(adv_preds)
asr = (np.sum(correct & changed) / np.sum(correct)) * 100
print(f"Attack Success Rate ASR {asr:.2f} percent")
```

---

### 2. GC Content Deviation

```python
from evaluation.gc_content import calculate_gc_content

gc1 = calculate_gc_content(original_sequence)
gc2 = calculate_gc_content(perturbed_sequence)
delta_gc = abs(gc1 - gc2)
print(f"Delta GC {delta_gc:.2f} percent")
```

---

### 3. Minimum Free Energy Change

```python
from evaluation.mfe import get_mfe

mfe1 = get_mfe(original_sequence)
mfe2 = get_mfe(perturbed_sequence)
delta_mfe = abs(mfe1 - mfe2)
print(f"Delta MFE {delta_mfe:.2f}")
```

---

## Folder Structure

```
attacks/
├── nucleotide.py
├── codon.py
├── synonymous_codon.py
├── GC_guide_synonymous_codon.py
└── backtranslation.py

evaluation/
├── asr.py
├── gc_content.py
└── mfe.py
```

---

