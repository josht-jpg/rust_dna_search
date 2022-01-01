#[repr(u8)]
#[derive(Copy, Clone, PartialEq, Debug, Eq, PartialOrd, Ord)]
enum Nucleotide {
    A = b'A',
    C = b'C',
    G = b'G',
    T = b'T',
}

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

fn str_to_gene(s: &str) -> Gene {
    let nucleotides = str_to_nucleotides(s);

    //Maybe take a functional approach?
    let mut gene: Gene = Vec::new();
    for i in (0..nucleotides.len()).step_by(3) {
        if i >= nucleotides.len() - 2 {
            return gene;
        }

        let codon: Codon = (nucleotides[i], nucleotides[i + 1], nucleotides[i + 2]);
        gene.push(codon);
    }

    return gene;
}

// Pass by reference?
fn linear_search_for_codon(gene: &Gene, key_codon: &Codon) -> bool {
    for codon in gene.iter() {
        print!("{:?}, {:?}\n", codon, key_codon);
        if codon == key_codon {
            return true;
        }
    }
    return false;
}

/*fn sort_codons_in_gene(gene: &Gene) -> Gene {
    gene.sort()
}*/

fn binary_search_for_codon(gene: &mut Gene, key_codon: &Codon) -> bool {
    gene.sort();

    println!("{:?}\n", gene);

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

/*fn binary_search_for_codon(gene: &mut Gene, key_codon: &Codon) -> bool {
    gene.sort();

    println!("{:?}\n", gene);

    let mut half_index = gene.len() / 2;
    let mut lower_slice: &[Codon] = &gene[..half_index];
    let mut upper_slice: &[Codon] = &gene[half_index..];

    while half_index > 0 {
        println!("{}\n {:?}\n {:?}\n\n", half_index, lower_slice, upper_slice);

        if key_codon == &gene[half_index] {
            return true;
        }

        if key_codon < &gene[half_index] {
            half_index /= 2;
            upper_slice = &lower_slice[half_index..];
            lower_slice = &lower_slice[..half_index];
        } else {
            //  half_index += (gene.len() - half_index) / 2;
            half_index /= 2;
            lower_slice = &upper_slice[..half_index];
            upper_slice = &upper_slice[half_index..];
        }
    }

    return false;
}*/

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
        //TODO: make an actual gene maybe?
        let gene_str = "AAATTTAACACAGTTGTGACACATTGTGTTCACCACGTGTC";

        // Try to make non immutatble
        let mut gene = str_to_gene(gene_str);

        let aaa: Codon = (Nucleotide::A, Nucleotide::A, Nucleotide::A);
        let ttt: Codon = (Nucleotide::T, Nucleotide::T, Nucleotide::T);
        let aca: Codon = (Nucleotide::A, Nucleotide::C, Nucleotide::A);
        let cag: Codon = (Nucleotide::C, Nucleotide::A, Nucleotide::G);

        //assert!(binary_search_for_codon(&mut gene, &aca));
        assert!(binary_search_for_codon(&mut gene, &aaa));
        assert!(binary_search_for_codon(&mut gene, &ttt));
        assert!(binary_search_for_codon(&mut gene, &aca));
        assert!(!binary_search_for_codon(&mut gene, &cag));
    }
}

/*fn string_to_gene(s: &str) -> Gene {
    let nucleotides: Vec<Nucleotide> = s
        .chars()
        .map(|c| match c {
            'A' => Nucleotide::A,
            'C' => Nucleotide::C,
            'G' => Nucleotide::G,
            'T' => Nucleotide::T,
            _ => panic!("Yeet"),
        })
        .collect();

    let mut gene: Gene = Vec::new();
    for i in (0..nucleotides.len()).step_by(3) {
        if i >= nucleotides.len() - 2 {
            return gene;
        }

        let codon: Codon = (nucleotides[i], nucleotides[i + 1], nucleotides[i + 2]);
        gene.push(codon);
    }

    return gene;
}*/
