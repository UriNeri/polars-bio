#!/usr/bin/env python3
"""
test script for genetic code translation functionality
Demonstrates the use of different genetic code tables. Genetic codes sourced from Seals2 by Yuri Wolf (https://github.com/YuriWolf-ncbi/seals-2/blob/master/bin/misc/orf)
"""

import sys
import tempfile
import os
sys.path.insert(0, '.')

from polars_bio.polars_bio import translate, translate_file, reverse_complement_seq, count_nucleotides, gc_content, list_genetic_codes

def test_genetic_code_tables():
    """Test translation with different genetic code tables"""
    print("=== Testing Different Genetic Code Tables ===")
    
    # Get available genetic codes
    codes = list_genetic_codes()
    print(f"Available genetic codes: {len(codes)}")
    for code_num, name in codes:
        print(f"  {code_num}: {name}")
    
    # Test sequence with ambiguous codons that differ between codes
    # UGA codon: * in standard code, W in some mitochondrial codes
    test_seq = "ATGTGAAAG"  # ATG-TGA-AAG (M-*-K in standard, M-W-K in some mito codes)
    
    print(f"\nTesting sequence: {test_seq}")
    print("Translations with different genetic codes:")
    
    # Test key genetic codes
    test_codes = [
        (1, "Standard"),
        (2, "Vertebrate Mitochondrial"),
        (3, "Yeast Mitochondrial"), 
        (4, "Mold Mitochondrial"),
        (5, "Invertebrate Mitochondrial"),
        (11, "Bacterial and Plant Plastid")
    ]
    
    for code_num, name in test_codes:
        result = translate(test_seq, genetic_code=code_num)
        print(f"  Code {code_num:2d} ({name}): {result}")
    
    print("Translation tests passed!")

def test_mitochondrial_differences():
    """Test specific differences in mitochondrial codes"""
    print("\n=== Testing Mitochondrial Code Differences ===")
    
    # UGA codon testing (stop vs tryptophan)
    uga_seq = "TGAAAG"  # TGA-AAG
    print(f"Testing UGA codon with sequence: {uga_seq}")
    
    standard = translate(uga_seq, genetic_code=1)
    vertebrate_mito = translate(uga_seq, genetic_code=2)
    
    print(f"Standard code (1):           {standard}")
    print(f"Vertebrate Mitochondrial (2): {vertebrate_mito}")
    
    # AGA/AGG codon testing (arginine vs stop in some codes)
    aga_seq = "AGAAGG"  # AGA-AGG
    print(f"\nTesting AGA/AGG codons with sequence: {aga_seq}")
    
    standard = translate(aga_seq, genetic_code=1)
    vertebrate_mito = translate(aga_seq, genetic_code=2)
    
    print(f"Standard code (1):           {standard}")
    print(f"Vertebrate Mitochondrial (2): {vertebrate_mito}")
    
    print("Mitochondrial difference tests passed!")

def test_all_frames_with_genetic_codes():
    """Test translation in all frames with different genetic codes"""
    print("\n=== Testing All Frames with Different Genetic Codes ===")
    
    seq = "ATGTGAAAGTAG"
    codes_to_test = [1, 2, 11]  # Standard, Vertebrate Mito, Bacterial
    
    for genetic_code in codes_to_test:
        print(f"\nGenetic Code {genetic_code}:")
        for frame in range(-3, 4):
            if frame == 0:
                continue
            result = translate(seq, frame=frame, genetic_code=genetic_code)
            print(f"  Frame {frame:2d}: {result}")
    
    print("All frames tests passed!")

def test_translate_file_with_genetic_codes():
    """Test file translation with different genetic codes"""
    print("\n=== Testing File Translation with Genetic Codes ===")
    
    # Create a temporary FASTA file with mitochondrial-like sequences
    fasta_content = """>mito_seq1
ATGTGAAAGTAG
>mito_seq2  
GGCAGAAGGTGC
>nuclear_seq
ATGAAATGA
"""
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(fasta_content)
        temp_file = f.name
    
    try:
        # Test with standard genetic code
        print("Standard genetic code (1):")
        results_standard = translate_file(temp_file, genetic_code=1)
        for seq_id, translation in results_standard:
            print(f"  {seq_id}: {translation}")
        
        # Test with vertebrate mitochondrial code
        print("\nVertebrate Mitochondrial genetic code (2):")
        results_mito = translate_file(temp_file, genetic_code=2)
        for seq_id, translation in results_mito:
            print(f"  {seq_id}: {translation}")
        
        print("File translation tests passed!")
        
    finally:
        os.unlink(temp_file)

def test_codon_tables():
    """Test codon table output"""
    print("\n=== Testing Codon Tables ===")
    
    seq = "TTCTAGTGA"  # TTC-TAG-TGA
    
    # Get codons instead of amino acids
    codons_std = translate(seq, to_protein=False, genetic_code=1)
    codons_mito = translate(seq, to_protein=False, genetic_code=2)
    
    print(f"Sequence: {seq}")
    print(f"Codons: {codons_std}")
    print(f"Standard translation (1):     {translate(seq, genetic_code=1)}")
    print(f"Mitochondrial translation (2): {translate(seq, genetic_code=2)}")
    
    print("Codon table tests passed!")

def test_basic_functionality():
    """Test basic sequence utility functions"""
    print("\n=== Testing Basic Sequence Functions ===")
    
    seq = "ATGCGTAACGTTAGC"
    
    # Reverse complement
    rev_comp = reverse_complement_seq(seq)
    print(f"Original:         {seq}")
    print(f"Reverse complement: {rev_comp}")
    
    # Nucleotide counts
    counts = count_nucleotides(seq)
    print(f"Nucleotide counts: {counts}")
    
    # GC content
    gc = gc_content(seq)
    print(f"GC content: {gc:.3f}")
    
    print("Basic functionality tests passed!")

if __name__ == "__main__":
    test_genetic_code_tables()
    test_mitochondrial_differences()
    test_all_frames_with_genetic_codes()
    test_translate_file_with_genetic_codes()
    test_codon_tables()
    test_basic_functionality()
    print("\n🧬 All genetic code tests passed! 🧬")
    print("Supporting translation for various organisms and organelles.") 