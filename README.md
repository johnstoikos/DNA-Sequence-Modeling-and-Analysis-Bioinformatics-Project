# DNA-Sequence-Modeling-and-Analysis-Bioinformatics-Project


## Overview
This project, developed in **Python**, focuses on the **modeling and analysis of DNA sequences** using core **bioinformatics methods**.  
The program automatically generates datasets, performs multiple sequence alignment, builds a profile Hidden Markov Model (HMM), and evaluates new DNA sequences using the **Viterbi algorithm**.
**The following project was created for a university class**

---

## Authors
- **Anargyrou Lamprou Aikaterini**  
- **Stoikos Ioannis Panagiotis**

---

## Project Stages

### Sequence Synthesis
- The DNA alphabet `['A', 'C', 'G', 'T']` and predefined **patterns** were used as the base.  
- Each sequence is generated with:
  - 1–3 random characters added at the start  
  - 1–2 random mutations (deletions or substitutions) applied to the patterns  
  - 1–2 random characters appended at the end  
- All sequences are collected and split into datasets using `random.sample()`:
  - **Dataset A** → 10 sequences  
  - **Dataset B** → 70 sequences  
  - **Dataset C** → 20 sequences  

---

### Multiple Sequence Alignment (MSA)
- Implemented using the **Needleman–Wunsch algorithm** with **NumPy** for dynamic programming.  
- Penalty values for match, mismatch, and gap were determined based on the last digit of the student ID (α = 1).  
- The **longest sequence from Dataset A** was chosen as the central sequence for the star alignment strategy.  
- The function `update_alignments()` was developed to:
  - Maintain global alignment consistency
  - Insert gaps dynamically when necessary
  - Align all sequences relative to the central one

---

### Profile Hidden Markov Model (HMM)
- The profile HMM represents **statistical behavior per position** in the alignment.  
- For each column:
  - The frequency of each base (A, C, G, T) is calculated.
  - Probabilities are computed by normalizing counts over the total number of sequences.
- The result is a **probability matrix** describing the DNA model per position.  
- The profile HMM is generated using **Dataset A** as input.

---

### Viterbi Algorithm & Sequence Evaluation
- The **Viterbi algorithm** was implemented **from scratch**, without using external libraries.  
- The algorithm computes the most probable state path for a given sequence according to the profile HMM.  
- For every sequence in **Dataset C** and for **20 random sequences**, the algorithm:
  - Calculates the **log-score** of the most probable path
  - Returns both the path and its score  
- This comparison demonstrates the model’s ability to **distinguish related sequences** from random ones.

---

## Results & Conclusions
Through this project, we gained practical understanding of:
- How **bioinformatics algorithms** (Needleman–Wunsch, HMMs, Viterbi) work in DNA analysis  
- The power of **profile HMMs** in identifying sequence similarity  
- The challenges of **alignment gap management** to maintain uniform sequence alignment

Overall, the implementation successfully demonstrates how bioinformatics tools can model, align, and evaluate biological sequences.

---

## Technologies Used
- **Language:** Python 3  
- **Libraries:** `numpy`, `random`  
- **Algorithms:** Needleman–Wunsch, Profile HMM, Viterbi  

---

## How to Run
```bash
python main.py
