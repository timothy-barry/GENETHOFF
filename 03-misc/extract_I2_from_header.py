#!/usr/bin/env python3
"""
Script pour extraire l'index d'un read FASTQ depuis le header
et remplacer la séquence par cet index.
Traite des fichiers FASTQ.gz en entrée et sortie.
"""

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import gzip
import sys


def extract_index_from_header(header):
    """
    Extrait l'index après le '+' dans le header FASTQ.
    
    Args:
        header: Description du header FASTQ
    
    Returns:
        L'index extrait ou None si non trouvé
    """
    if '+' in header:
        index = header.split('+', 1)[1].strip()
        return index
    return None


def process_fastq_gz(input_file, output_file):
    """
    Traite un fichier FASTQ.gz en remplaçant les séquences par leurs index.
    Écrit les reads au fur et à mesure (streaming).
    
    Args:
        input_file: Chemin du fichier FASTQ.gz d'entrée
        output_file: Chemin du fichier FASTQ.gz de sortie
    
    Returns:
        Nombre de reads traités
    """
    nb_with_index = 0
    nb_without_index = 0
    total_reads = 0
    
    # Compte d'abord le nombre total de reads
    print(f"Comptage des reads dans {input_file}...")
    with gzip.open(input_file, 'rt') as handle:
        for _ in SeqIO.parse(handle, "fastq"):
            total_reads += 1
    
    print(f"Total: {total_reads} reads à traiter\n")
    print(f"Traitement et écriture en cours...")
    
    # Traite et écrit les reads au fur et à mesure
    with gzip.open(input_file, 'rt') as input_handle, \
         gzip.open(output_file, 'wt') as output_handle:
        
        for i, record in enumerate(SeqIO.parse(input_handle, "fastq"), 1):
            # Affiche la progression tous les 1000 reads
            if i % 1000 == 0:
                percentage = (i / total_reads) * 100
                print(f"  Traités: {i}/{total_reads} reads ({percentage:.1f}%)...", end='\r')
            
            # Extrait l'index du header
            index_seq = extract_index_from_header(record.description)
            
            if index_seq:
                # Crée un nouveau record avec l'index comme séquence
                new_record = SeqRecord(
                    Seq(index_seq),
                    id=record.id,
                    name=record.name,
                    description=record.description,
                    letter_annotations={'phred_quality': [40] * len(index_seq)}
                )
                SeqIO.write(new_record, output_handle, "fastq")
                nb_with_index += 1
            else:
                # Écrit le record original si pas d'index trouvé
                SeqIO.write(record, output_handle, "fastq")
                nb_without_index += 1
    
    print(f"  Traités: {total_reads} reads (100.0%) - Terminé!  ")
    
    print(f"\nTraitement terminé:")
    print(f"  - {nb_with_index} reads avec index (séquence remplacée)")
    print(f"  - {nb_without_index} reads sans index (séquence conservée)")
    print(f"  - Total: {total_reads} reads")
    
    return total_reads


if __name__ == "__main__":
    # Vérification des arguments
    if len(sys.argv) != 3:
        print("Usage: python script.py <input.fastq.gz> <output.fastq.gz>")
        print("\nExemple:")
        print("  python script.py input.fastq.gz output.fastq.gz")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    # Traitement du fichier
    try:
        process_fastq_gz(input_file, output_file)
    except FileNotFoundError:
        print(f"Erreur: Le fichier {input_file} n'existe pas.")
        sys.exit(1)
    except Exception as e:
        print(f"Erreur lors du traitement: {e}")
        sys.exit(1)
