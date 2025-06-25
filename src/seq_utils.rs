use std::collections::HashMap;
use pyo3::prelude::*;
use needletail::parse_fastx_file;

// Genetic codes / variables sourced from Seals2 by Yuri Wolf (https://github.com/YuriWolf-ncbi/seals-2/blob/master/bin/misc/orf)
// Original Perl implementation for comprehensive genetic code support

/// Genetic code amino acid sequences (64 codons in order: TTT, TTC, TTA, TTG, TCT, TCC, TCA, TCG, TAT, TAC, TAA, TAG, TGT, TGC, TGA, TGG, CTT, CTC, CTA, CTG, CCT, CCC, CCA, CCG, CAT, CAC, CAA, CAG, CGT, CGC, CGA, CGG, ATT, ATC, ATA, ATG, ACT, ACC, ACA, ACG, AAT, AAC, AAA, AAG, AGT, AGC, AGA, AGG, GTT, GTC, GTA, GTG, GCT, GCC, GCA, GCG, GAT, GAC, GAA, GAG, GGT, GGC, GGA, GGG)
fn get_genetic_code_table() -> HashMap<u32, &'static str> {
    let mut code_table = HashMap::new();
    
    code_table.insert(1, "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG");
    code_table.insert(2, "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG");
    code_table.insert(3, "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG");
    code_table.insert(4, "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG");
    code_table.insert(5, "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG");
    code_table.insert(6, "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG");
    code_table.insert(9, "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG");
    code_table.insert(10, "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG");
    code_table.insert(11, "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG");
    code_table.insert(12, "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG");
    code_table.insert(13, "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG");
    code_table.insert(14, "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG");
    code_table.insert(15, "FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG");
    code_table.insert(16, "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG");
    code_table.insert(21, "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG");
    code_table.insert(22, "FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG");
    code_table.insert(23, "FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG");
    code_table.insert(24, "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG");
    code_table.insert(25, "FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG");
    code_table.insert(26, "FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG");
    code_table.insert(27, "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG");
    code_table.insert(28, "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG");
    code_table.insert(29, "FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG");
    code_table.insert(30, "FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG");
    code_table.insert(31, "FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG");
    code_table.insert(33, "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG");
    
    code_table
}

// /// Genetic code start codon patterns (M for start codons, - for others)
// fn get_start_codon_table() -> HashMap<u32, &'static str> {
//     let mut start_table = HashMap::new();
//     
//     start_table.insert(1, "---M---------------M---------------M----------------------------");
//     start_table.insert(2, "--------------------------------MMMM---------------M------------");
//     start_table.insert(3, "----------------------------------MM----------------------------");
//     start_table.insert(4, "--MM---------------M------------MMMM---------------M------------");
//     start_table.insert(5, "---M----------------------------MMMM---------------M------------");
//     start_table.insert(6, "-----------------------------------M----------------------------");
//     start_table.insert(9, "-----------------------------------M---------------M------------");
//     start_table.insert(10, "-----------------------------------M----------------------------");
//     start_table.insert(11, "---M---------------M------------MMMM---------------M------------");
//     start_table.insert(12, "-------------------M---------------M----------------------------");
//     start_table.insert(13, "-----------------------------------M----------------------------");
//     start_table.insert(14, "-----------------------------------M----------------------------");
//     start_table.insert(15, "-----------------------------------M----------------------------");
//     start_table.insert(16, "-----------------------------------M----------------------------");
//     start_table.insert(21, "-----------------------------------M---------------M------------");
//     start_table.insert(22, "-----------------------------------M----------------------------");
//     start_table.insert(23, "--------------------------------M--M---------------M------------");
//     start_table.insert(24, "---M---------------M---------------M---------------M------------");
//     start_table.insert(25, "---M-------------------------------M---------------M------------");
//     start_table.insert(26, "-------------------M---------------M----------------------------");
//     start_table.insert(27, "-----------------------------------M----------------------------");
//     start_table.insert(28, "-----------------------------------M----------------------------");
//     start_table.insert(29, "-----------------------------------M----------------------------");
//     start_table.insert(30, "-----------------------------------M----------------------------");
//     start_table.insert(31, "-----------------------------------M----------------------------");
//     start_table.insert(33, "---M---------------M---------------M---------------M------------");
//     
//     start_table
// }

/// Genetic code names
fn get_genetic_code_names() -> HashMap<u32, &'static str> {
    let mut names = HashMap::new();
    
    names.insert(1, "Standard");
    names.insert(2, "Vertebrate Mitochondrial");
    names.insert(3, "Yeast Mitochondrial");
    names.insert(4, "Mold Mitochondrial; Protozoan Mitochondrial; Coelenterate Mitochondrial; Mycoplasma; Spiroplasma");
    names.insert(5, "Invertebrate Mitochondrial");
    names.insert(6, "Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear");
    names.insert(9, "Echinoderm Mitochondrial; Flatworm Mitochondrial");
    names.insert(10, "Euplotid Nuclear");
    names.insert(11, "Bacterial and Plant Plastid");
    names.insert(12, "Alternative Yeast Nuclear");
    names.insert(13, "Ascidian Mitochondrial");
    names.insert(14, "Alternative Flatworm Mitochondrial");
    names.insert(15, "Blepharisma Macronuclear");
    names.insert(16, "Chlorophycean Mitochondrial");
    names.insert(21, "Trematode Mitochondrial");
    names.insert(22, "Scenedesmus obliquus mitochondrial");
    names.insert(23, "Thraustochytrium mitochondrial");
    names.insert(24, "Pterobranchia mitochondrial");
    names.insert(25, "Candidate Division SR1 and Gracilibacteria");
    names.insert(26, "Pachysolen tannophilus Nuclear Code");
    names.insert(27, "Karyorelict Nuclear Code");
    names.insert(28, "Condylostoma Nuclear Code");
    names.insert(29, "Mesodinium Nuclear Code");
    names.insert(30, "Peritrich Nuclear Code");
    names.insert(31, "Blastocrithidia Nuclear Code");
    names.insert(33, "Cephalodiscidae Mitochondrial UAA-Tyr Code");
    
    names
}

/// Convert codon index to amino acid using genetic code table
fn codon_to_amino_acid(codon: &str, genetic_code_table: &str) -> char {
    // Codon order: TTT, TTC, TTA, TTG, TCT, TCC, TCA, TCG, TAT, TAC, TAA, TAG, TGT, TGC, TGA, TGG, 
    //              CTT, CTC, CTA, CTG, CCT, CCC, CCA, CCG, CAT, CAC, CAA, CAG, CGT, CGC, CGA, CGG,
    //              ATT, ATC, ATA, ATG, ACT, ACC, ACA, ACG, AAT, AAC, AAA, AAG, AGT, AGC, AGA, AGG,
    //              GTT, GTC, GTA, GTG, GCT, GCC, GCA, GCG, GAT, GAC, GAA, GAG, GGT, GGC, GGA, GGG
    
    let bases = ['T', 'C', 'A', 'G'];
    let codon_chars: Vec<char> = codon.chars().collect();
    
    if codon_chars.len() != 3 {
        return 'X';
    }
    
    let pos1 = bases.iter().position(|&b| b == codon_chars[0]).unwrap_or(0);
    let pos2 = bases.iter().position(|&b| b == codon_chars[1]).unwrap_or(0);
    let pos3 = bases.iter().position(|&b| b == codon_chars[2]).unwrap_or(0);
    
    let index = pos1 * 16 + pos2 * 4 + pos3;
    
    genetic_code_table.chars().nth(index).unwrap_or('X')
}

/// Generate reverse complement of a DNA sequence
fn reverse_complement(sequence: &str) -> String {
    sequence
        .chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A', 
            'C' => 'G',
            'G' => 'C',
            _ => c,
        })
        .collect()
}

/// Translate a DNA or RNA sequence into amino acid sequence
/// 
/// # Arguments
/// * `sequence` - DNA or RNA sequence string
/// * `frame` - Reading frame (-3, -2, -1, 1, 2, 3). Negative frames use reverse complement
/// * `start` - Start position in the sequence (0-based)
/// * `stop` - Stop position in the sequence (0-based, exclusive). If None, uses end of sequence
/// * `to_protein` - If true, return amino acid sequence. If false, return codons
/// * `genetic_code` - Genetic code table number (default: 1 for standard code)
#[pyfunction]
#[pyo3(signature = (sequence, frame=1, start=0, stop=None, to_protein=true, genetic_code=1))]
pub fn translate(
    sequence: String,
    frame: i8,
    start: usize,
    stop: Option<usize>,
    to_protein: bool,
    genetic_code: u32,
) -> PyResult<String> {
    if sequence.is_empty() {
        return Ok(String::new());
    }
    
    // Convert to uppercase and replace U with T for RNA
    let mut seq = sequence.to_uppercase().replace('U', "T");
    
    // Validate sequence contains only valid nucleotides
    for c in seq.chars() {
        if !matches!(c, 'A' | 'T' | 'C' | 'G' | 'U') {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                "Sequence contains invalid nucleotides. Only A, T, C, G, U are allowed."
            ));
        }
    }
    
    // Apply start and stop positions  
    let end_pos = stop.unwrap_or(seq.len());
    if start < seq.len() && start < end_pos {
        seq = seq[start..end_pos.min(seq.len())].to_string();
    } else {
        return Ok(String::new());
    }
    
    // Handle negative frames (reverse complement)
    let frame_abs = if frame < 0 {
        seq = reverse_complement(&seq);
        frame.abs()
    } else {
        frame
    };
    
    // Validate frame is in range 1-3
    if !(1..=3).contains(&frame_abs) {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            "Frame must be between -3 and 3 (excluding 0)."
        ));
    }
    
    // Adjust for reading frame (1-based to 0-based)
    let frame_offset = (frame_abs - 1) as usize;
    if frame_offset >= seq.len() {
        return Ok(String::new());
    }
    
    seq = seq[frame_offset..].to_string();
    
    // Get genetic code table
    let code_tables = get_genetic_code_table();
    let genetic_code_str = code_tables.get(&genetic_code).unwrap_or(&code_tables[&1]);
    
    // Translate sequence
    let mut result = Vec::new();
    let mut codons = Vec::new();
    
    for i in (0..seq.len()).step_by(3) {
        if i + 2 < seq.len() {
            let codon = &seq[i..i+3];
            codons.push(codon.to_string());
            let amino_acid = codon_to_amino_acid(codon, genetic_code_str);
            result.push(amino_acid);
        }
    }
    
    if to_protein {
        Ok(result.into_iter().collect())
    } else {
        Ok(codons.join(" "))
    }
}

/// Translate sequences from a FASTA/FASTQ file
/// 
/// # Arguments
/// * `file_path` - Path to the FASTA/FASTQ file
/// * `frame` - Reading frame (-3, -2, -1, 1, 2, 3)
/// * `to_protein` - If true, return amino acid sequences. If false, return codons
/// * `genetic_code` - Genetic code table number (default: 1 for standard code)
#[pyfunction]
#[pyo3(signature = (file_path, frame=1, to_protein=true, genetic_code=1))]
pub fn translate_file(
    file_path: String,
    frame: i8,
    to_protein: bool,
    genetic_code: u32,
) -> PyResult<Vec<(String, String)>> {
    let mut results = Vec::new();
    
    let mut reader = parse_fastx_file(&file_path)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Error reading file: {}", e)))?;
    
    while let Some(record) = reader.next() {
        let record = record
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Error parsing record: {}", e)))?;
        
        let id = String::from_utf8_lossy(record.id()).to_string();
        let sequence = String::from_utf8_lossy(&record.seq()).to_string();
        
        let translated = translate(sequence, frame, 0, None, to_protein, genetic_code)?;
        results.push((id, translated));
    }
    
    Ok(results)
}

/// Get available genetic code tables
#[pyfunction]
pub fn list_genetic_codes() -> PyResult<Vec<(u32, String)>> {
    let names = get_genetic_code_names();
    let mut codes: Vec<(u32, String)> = names.iter()
        .map(|(&code, &name)| (code, name.to_string()))
        .collect();
    codes.sort_by_key(|&(code, _)| code);
    Ok(codes)
}

/// Get reverse complement of a DNA sequence
#[pyfunction]
pub fn reverse_complement_seq(sequence: String) -> PyResult<String> {
    let seq = sequence.to_uppercase().replace('U', "T");
    
    // Validate sequence contains only valid nucleotides
    for c in seq.chars() {
        if !matches!(c, 'A' | 'T' | 'C' | 'G') {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                "Sequence contains invalid nucleotides. Only A, T, C, G, U are allowed."
            ));
        }
    }
    
    Ok(reverse_complement(&seq))
}

/// Count nucleotides in a sequence
#[pyfunction]
pub fn count_nucleotides(sequence: String) -> PyResult<HashMap<char, usize>> {
    let mut counts = HashMap::new();
    
    for c in sequence.to_uppercase().chars() {
        if matches!(c, 'A' | 'T' | 'C' | 'G' | 'U' | 'N') {
            *counts.entry(c).or_insert(0) += 1;
        }
    }
    
    Ok(counts)
}

/// Calculate GC content of a sequence
#[pyfunction]
pub fn gc_content(sequence: String) -> PyResult<f64> {
    let seq = sequence.to_uppercase();
    let total = seq.len() as f64;
    let gc_count = seq.chars().filter(|&c| c == 'G' || c == 'C').count() as f64;
    
    if total == 0.0 {
        Ok(0.0)
    } else {
        Ok(gc_count / total)
    }
} 