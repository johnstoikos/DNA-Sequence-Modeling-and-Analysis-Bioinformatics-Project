import random
import numpy as np
from collections import Counter, defaultdict

# θέμα 1

import random
# εδώ δημιουργούμε τους αρχικούς πίνακες ορίζουμε το αλφάβητο και τα μοτίβα
alphavito = ['A', 'C', 'G', 'T']
patterns = ["ATTAGA","ACGCATTT","AGGACTCAA","ATTTCAGT"]

# η συναρτηση κανει αλλαγη στο μοτιβο, είτε με διαγραφη είτε με αντικατάσταη
def metavoli_pat(pattern, megisto_metavoles=2):
    pattern = list(pattern)
    plithos_metavolon = random.randint(0, megisto_metavoles)
    theseis = random.sample(range(len(pattern)), plithos_metavolon) if plithos_metavolon > 0 else []
    for thesi in theseis:
        energia = random.choice(["Αντικατάσταση", "Διαγραφή"])
        if energia == "Αντικατάσταση":
            trexon = pattern[thesi]
            neo_simvolo = random.choice([a for a in alphavito if a != trexon])
            pattern[thesi] = neo_simvolo
        elif energia == "Διαγραφή":
            pattern[thesi] = "" 
    return "".join([s for s in pattern if s != ""])

def sinthetise_akolouthia():
    prothima = "".join(random.choices(alphavito, k=random.randint(1, 3)))
    akolouthia = prothima

    for p in patterns:
        metavlito = metavoli_pat(p)
        akolouthia += metavlito
    epithima = "".join(random.choices(alphavito, k=random.randint(1, 2)))
    akolouthia += epithima
    return akolouthia

oles_akol = [sinthetise_akolouthia() for _ in range(100)]

# εδω βάζουμε όλες τις ακολουθιες σε μια λιστα και μετα ξεχερωσιτα σε άλλες τρεις.
random.shuffle(oles_akol)
datasetA = oles_akol[:10]
datasetB = oles_akol[10:80]
datasetC = oles_akol[80:]

print("datasetA:")
for akolouthia in datasetA:
    print(akolouthia)
print("\ndatasetB:", len(datasetB))
print("datasetC:", len(datasetC))



# θέμα 2

def zeugos_stoixisi(akolouthia1, akolouthia2, antistoixisi=1, timi_mh_antistoixisis=0.5, timi_kenou=1):
   
    n, m = len(akolouthia1), len(akolouthia2)
   
    bathmologia = np.zeros((n+1, m+1))
    deiktis = np.zeros((n+1, m+1), dtype=int)

    for i in range(1, n+1):
        bathmologia[i][0] = bathmologia[i-1][0] - timi_kenou
        deiktis[i][0] = 1  
    for j in range(1, m+1):
        bathmologia[0][j] = bathmologia[0][j-1] - timi_kenou
        deiktis[0][j] = 2  #

    for i in range(1, n+1):
        for j in range(1, m+1):
            if akolouthia1[i-1] == akolouthia2[j-1]:
                diag = bathmologia[i-1][j-1] + antistoixisi
            else:
                diag = bathmologia[i-1][j-1] - timi_mh_antistoixisis
            panw = bathmologia[i-1][j] - timi_kenou
            aristera = bathmologia[i][j-1] - timi_kenou
            epiloges = [diag, panw, aristera]
            bathmologia[i][j] = max(epiloges)
            deiktis[i][j] = epiloges.index(bathmologia[i][j])


    stoixismeni1, stoixismeni2 = "", ""
    
    i, j = n, m
    
    while i > 0 or j > 0:
        if i > 0 and j > 0 and deiktis[i][j] == 0:
            stoixismeni1 = akolouthia1[i-1] + stoixismeni1
            stoixismeni2 = akolouthia2[j-1] + stoixismeni2
            i -= 1
            j -= 1
        elif i > 0 and (j == 0 or deiktis[i][j] == 1):
            stoixismeni1 = akolouthia1[i-1] + stoixismeni1
            stoixismeni2 = "-" + stoixismeni2
            i -= 1
        else:
            stoixismeni1 = "-" + stoixismeni1
            stoixismeni2 = akolouthia2[j-1] + stoixismeni2
            j -= 1
    return stoixismeni1, stoixismeni2




# πολλαπλή στοιχιση με βάση την πρωτη ακολουθία

def poll_akolouthies(akol, antistoixisi=1, timi_antistoix=0.5, timi_kenou=1):
    stoix_akol = [akol[0]]
    for i in range(1, len(akol)):
        symfonia = stoix_akol[0]
        nextAkolouthia = akol[i]
        stoixismeni1, stoixismeni2 = zeugos_stoixisi(symfonia, nextAkolouthia, antistoixisi, timi_antistoix, timi_kenou)

        nees_stoixismenes = []
        idx1 = 0
        for x in stoixismeni1:
            if x == "-":
                for j in range(len(stoix_akol)):
                    stoix_akol[j] = stoix_akol[j][:idx1] + "-" + stoix_akol[j][idx1:]
                idx1 += 1
            else:
                idx1 += 1
        stoix_akol.append(stoixismeni2)
    return stoix_akol

a = 1
antistoixisi = 1
timi_antistoix = a / 2
timi_kenou = a

stoixismenes = poll_akolouthies(datasetA, antistoixisi, timi_antistoix, timi_kenou)

print("\nΑποτέλεσμα πολλαπλής Στοίχησης:\n")


for i, akolouthia in enumerate(stoixismenes):
    print(f"Ακολουθία{i+1}: {akolouthia}")


# θέμα 3

def ftiakzeHmm(dedomena):
    mikos = max(len(akolouthia) for akolouthia in dedomena)
    profil = [defaultdict(int) for _ in range(mikos)]
    for akolouthia in dedomena:
        for i, s in enumerate(akolouthia):
            profil[i][s] += 1
    hmm = []
    for i, metriths in enumerate(profil):
        synolo = sum(metriths[a] for a in alphavito)
        pithanotites = {a: (metriths[a] / synolo if synolo > 0 else 0.0) for a in alphavito}
        hmm.append(pithanotites)
    return hmm

# Πριν κατασκευάσω το HMM, φροντίζω όλα τα alignments να έχουν ίδιο μήκος
max_len = max(len(seq) for seq in stoixismenes)
for i in range(len(stoixismenes)):
    while len(stoixismenes[i]) < max_len:
        stoixismenes[i] += '-'

hmm_profil = ftiakzeHmm(datasetA)

def hmmPinakas(hmm):
    print("Thesi | " + " | ".join(alphavito))
    for i, thesi in enumerate(hmm):
        grammi = [f"{thesi[a]:.2f}" for a in alphavito]
        print(f"{i+1:5d} | " + " | ".join(grammi))

print("HMM profile:")
hmmPinakas(hmm_profil)


# θέμα 4

def viterbi_akol(akolouthia, hmm):
    n = len(akolouthia)
    m = len(hmm)
    V = np.zeros((n, m))
    monopati = np.zeros((n, m), dtype=int)
    
    for j in range(m):
        V[0][j] = hmm[j].get(akolouthia[0], 1e-6)
        monopati[0][j] = 0
    
    for i in range(1, n):
        for j in range(m):
            lista_pith = []
            for k in range(m):
                pith = V[i-1][k] * hmm[j].get(akolouthia[i], 1e-6)
                lista_pith.append((pith, k))
            V[i][j], monopati[i][j] = max(lista_pith)
    
    best_monopati = []
    teleutaio = np.argmax(V[n-1])
    for i in range(n-1, -1, -1):
        best_monopati.append(teleutaio)
        teleutaio = monopati[i][teleutaio]
    best_monopati.reverse()
    return best_monopati, np.max(V[n-1])

print("\nViterbi αποτελέσματα για όλες τις ακολουθίες του datasetC:")
bathmologies_C = []
for i, akolouthia in enumerate(datasetC):
    monopati, bathmos = viterbi_akol(akolouthia, hmm_profil)
    print(f"ΑκολουθίαC_{i+1}: {akolouthia} | Καλύτερο Μονοπάτι: {monopati} | Βαθμός: {bathmos:.5f}\n")
    bathmologies_C.append(bathmos)

print("\nAlignment score για 20 τυχαίες ακολουθίες από datasetC:")

tuxaioi_deiktes = random.sample(range(len(datasetC)), 20)
vathmoi_C = []
for idx, i in enumerate(tuxaioi_deiktes):
    akolouthia = datasetC[i]
    _, bathmos = viterbi_akol(akolouthia, hmm_profil)
    print(f"ΑκολουθίαC_δείγμα_{idx+1}: {akolouthia} | Alignment score: {bathmos:.5f}")
    vathmoi_C.append(bathmos)

mesosBathmos_C = sum(vathmoi_C) / len(vathmoi_C)
print(f"\nΜεσο alignment score για 20 τυχαίες ακολουθιες C: {mesosBathmos_C:.5f}")
