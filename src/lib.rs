#[repr(u8)]
#[derive(Copy, Clone, PartialEq, Debug, Eq, PartialOrd, Ord)]
enum Nucleotide {
    A = b'A',
    C = b'C',
    G = b'G',
    T = b'T',
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
            _ => panic!("Yeet"),
        })
        .collect()
}

fn only_nucleotides_in_codon(nucleotides: &[Nucleotide]) -> &[Nucleotide] {
    &nucleotides[0..(nucleotides.len() - (nucleotides.len() % NUCLEOTIDES_IN_CODON))]
}

fn str_to_gene(s: &str) -> Gene {
    let nucleotides = str_to_nucleotides(s);
    let num_nucleotides = only_nucleotides_in_codon(&nucleotides).len();

    return (0..num_nucleotides)
        .step_by(NUCLEOTIDES_IN_CODON)
        .map(|i| (nucleotides[i], nucleotides[i + 1], nucleotides[i + 2]))
        .collect();
}

fn linear_search_for_codon(gene: &Gene, key_codon: &Codon) -> bool {
    for codon in gene.iter() {
        if codon == key_codon {
            return true;
        }
    }
    return false;
}

fn binary_search_for_codon(gene: &mut Gene, key_codon: &Codon) -> bool {
    gene.sort();

    let mut low = 0;
    let mut high = gene.len() - 1;

    while low <= high {
        let mid = (low + high) / 2;

        if &gene[mid] < key_codon {
            low = mid + 1;
        } else if &gene[mid] > key_codon {
            high = mid - 1;
        } else {
            return true;
        }
    }

    return false;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn linear_search() {
        //TODO: make an actual gene maybe?
        let gene_str = "AACACAGTTGTGACACATTGTGTTCACCACGTGTC";
        let gene = str_to_gene(gene_str);

        let cat = (Nucleotide::C, Nucleotide::A, Nucleotide::T);
        let tgc = (Nucleotide::T, Nucleotide::G, Nucleotide::C);

        assert!(linear_search_for_codon(&gene, &cat));
        assert!(!linear_search_for_codon(&gene, &tgc));
    }

    #[test]
    fn binary_search() {
        // Snippet of the human genome: ATATCTTAGAGGGAGGGCTGAGGGTTTGAAGTCC
        let gene_str = "AAATTTAACACAGTTGTGACACATTGTGTTCACCACGTGTC";

        // Try to make non immutatble
        let mut gene = str_to_gene(gene_str);

        let aaa: Codon = (Nucleotide::A, Nucleotide::A, Nucleotide::A);
        let ttt: Codon = (Nucleotide::T, Nucleotide::T, Nucleotide::T);
        let aca: Codon = (Nucleotide::A, Nucleotide::C, Nucleotide::A);
        let cag: Codon = (Nucleotide::C, Nucleotide::A, Nucleotide::G);

        assert!(binary_search_for_codon(&mut gene, &aaa));
        assert!(binary_search_for_codon(&mut gene, &ttt));
        assert!(binary_search_for_codon(&mut gene, &aca));
        assert!(!binary_search_for_codon(&mut gene, &cag));
    }
}
