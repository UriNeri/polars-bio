#!/usr/bin/env python3
"""
Simple test script for the translate function
"""

import sys
sys.path.insert(0, '.')

from polars_bio import translate

def test_translate():
    # Test basic translation
    seq = "ATGAAATAG"  # Start codon (ATG), Lysine (AAA), Stop codon (TAG)
    result = translate(seq)
    print(f"Basic translation: {seq} -> {result}")
    assert result == "MK*", f"Expected 'MK*', got '{result}'"
    
    # Test different frames
    seq = "ATGAAATAG"
    for frame in [1, 2, 3]:
        result = translate(seq, frame=frame)
        print(f"Frame {frame}: {seq} -> {result}")
    
    # Test negative frames (reverse complement)
    for frame in [-1, -2, -3]:
        result = translate(seq, frame=frame)
        print(f"Frame {frame}: {seq} -> {result}")
    
    # Test with start and stop positions
    seq = "AAATGAAATAGTTT"
    result = translate(seq, start=3, stop=12)  # Extract "ATGAAATAG"
    print(f"With start/stop: {seq}[3:12] -> {result}")
    assert result == "MK*", f"Expected 'MK*', got '{result}'"
    
    # Test RNA sequence (with U)
    rna_seq = "AUGAAAUAG"
    result = translate(rna_seq)
    print(f"RNA translation: {rna_seq} -> {result}")
    assert result == "MK*", f"Expected 'MK*', got '{result}'"
    
    # Test getting codons instead of amino acids
    result = translate(seq, to_protein=False)
    print(f"Codons: {seq} -> {result}")
    
    print("\nAll tests passed!")

if __name__ == "__main__":
    test_translate() 