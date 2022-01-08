use std::collections::HashMap;

#[derive(Hash, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum Nucleotide {
    A,
    C,
    G,
    T,
}

const NUCLEOTIDES_IN_CODON: usize = 3;

type Codon = (Nucleotide, Nucleotide, Nucleotide);
type Gene = Vec<Codon>;

fn str_to_nucleotides(s: &str) -> Vec<Nucleotide> {
    s.chars()
        .map(|c| match c {
            'A' => Nucleotide::A,
            'C' => Nucleotide::C,
            'G' => Nucleotide::G,
            'T' => Nucleotide::T,
            _ => panic!("Not a nucleotide"),
        })
        .collect()
}

fn nucleotide_frequency(dna: &[Nucleotide]) -> HashMap<Nucleotide, u32> {
    let mut frequencies: HashMap<Nucleotide, u32> = HashMap::from([
        (Nucleotide::A, 0),
        (Nucleotide::T, 0),
        (Nucleotide::G, 0),
        (Nucleotide::C, 0),
    ]);

    for nucleotide in dna {
        increment_nucleotide_count(&mut frequencies, nucleotide);
    }

    return frequencies;
}

fn binary_search_for_codon(gene: &mut Gene, target_codon: &Codon) -> bool {
    gene.sort();

    let mut low = 0;
    let mut high = gene.len() - 1;

    while low <= high {
        let middle_codon_index = (low + high) / 2;
        let middle_codon = &gene[middle_codon_index];

        if middle_codon < target_codon {
            low = middle_codon_index + 1;
        } else if middle_codon > target_codon {
            high = middle_codon_index - 1;
        } else {
            return true;
        }
    }

    return false;
}

fn str_to_gene(s: &str) -> Gene {
    let nucleotides = str_to_nucleotides(s);
    let num_nucleotides_in_codons = nucleotides.len() - (nucleotides.len() % NUCLEOTIDES_IN_CODON);

    return (0..num_nucleotides_in_codons)
        .step_by(NUCLEOTIDES_IN_CODON)
        .map(|i| (nucleotides[i], nucleotides[i + 1], nucleotides[i + 2]))
        .collect();
}

fn increment_nucleotide_count(frequencies: &mut HashMap<Nucleotide, u32>, nucleotide: &Nucleotide) {
    *frequencies.get_mut(nucleotide).unwrap() += 1;
}

fn naive_match(dna: &[Nucleotide], target_sequence: &[Nucleotide]) -> Vec<usize> {
    assert!(dna.len() >= target_sequence.len());

    let mut target_occurance_at = Vec::new();

    let possible_starting_positions = 0..(dna.len() - target_sequence.len() + 1);

    for i in possible_starting_positions {
        for j in 0..target_sequence.len() {
            if dna[i + j] != target_sequence[j] {
                break;
            }

            if j == target_sequence.len() - 1 {
                target_occurance_at.push(i);
            }
        }
    }

    return target_occurance_at;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn binary_search() {
        let gene_str = "ATATCTTAGAGGGAGGGCTGAGGGTTTGAAGTCC";

        // Try to make non immutatble
        let mut gene = str_to_gene(gene_str);

        let ata: Codon = (Nucleotide::A, Nucleotide::T, Nucleotide::A);
        let atc: Codon = (Nucleotide::A, Nucleotide::T, Nucleotide::C);
        let agg: Codon = (Nucleotide::A, Nucleotide::G, Nucleotide::G);
        let tcc: Codon = (Nucleotide::T, Nucleotide::C, Nucleotide::C);

        assert!(binary_search_for_codon(&mut gene, &ata));
        assert!(!binary_search_for_codon(&mut gene, &atc));
        assert!(binary_search_for_codon(&mut gene, &agg));
        assert!(!binary_search_for_codon(&mut gene, &tcc));
    }

    #[test]
    fn nucleotide_frequency_test() {
        let dna_str = "ATATCTTAGAGGGAG";
        let freq = nucleotide_frequency(&str_to_nucleotides(&dna_str));

        assert_eq!(freq.get(&Nucleotide::A), Some(&5));
        assert_eq!(freq.get(&Nucleotide::T), Some(&4));
        assert_eq!(freq.get(&Nucleotide::C), Some(&1));
        assert_eq!(freq.get(&Nucleotide::G), Some(&5));
    }

    #[test]
    fn naive_match_test() {
        let dna = str_to_nucleotides("ATATCTTAGAGGGAGGGAGG");
        let target_sequence = str_to_nucleotides("AGG");

        let match_indices = naive_match(&dna, &target_sequence);

        assert_eq!(match_indices.len(), 3);
        assert_eq!(match_indices[0], 9);
        assert_eq!(match_indices[1], 13);
        assert_eq!(match_indices[2], 17);
    }
}
