package org.stepic.bioinformatics;

import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static java.util.stream.Collectors.toList;
import static org.stepic.bioinformatics.Probability.computeEntropy;

public class  GenomeList extends ArrayList<FuzzyGenome> {
    public Sequence medianString(int k){
        return
            Stream.iterate((long)0,i->i+1)
                .limit(Probability.pow(4,k))
                .map(i -> Genome.numberToPattern(i,k))
                .min(Comparator.comparingLong(this::distanceBetweenPatternAndThis))
                .orElse(new Sequence(""));
    }


    public int distanceBetweenPatternAndThis(Sequence pattern) {
        final var k = pattern.length();
        Function <FuzzyGenome,Integer> findMinDistanceWithThePattern =
                (text) -> text.kmerSequenceStream(k)
                        .map(pattern::hammingDistance)
                        .min(Integer::compareTo).orElse(Integer.MAX_VALUE);
        return this.stream()
                .map(findMinDistanceWithThePattern)
                .reduce(0,Integer::sum);
    }

    public GenomeList(final String... texts){
        this.addAll(Arrays.asList(texts).stream()
                .map(str -> new FuzzyGenome(str))
                .collect(toList()));
    }

    public Stream<Map<Character,Long>> countMotifsToStream(){
        var length = this.size() > 0 ? this.get(0).length() : 0;
        return Stream.iterate(0,i->i+1).limit(length)
                .map( i -> this.stream()
                        .filter(g -> g.length() == length)
                        .map( g -> g.charAt(i))
                        .collect(Collectors.groupingBy(
                                Function.identity(),
                                Collectors.counting()
                        ))
                );
    }

    public ProfileList computeProfile() {
        return countMotifsToStream()
                .map( (Map<Character,Long> column) ->{
                    var profile = new HashMap<Character,Double>(column.size());
                    for (char key :
                            new char[]{'A', 'C', 'G', 'T'}) {
                        double value = column.containsKey(key) ? column.get(key) : 0.0;
                        value = value/(double)this.size();
                        profile.put( key, value );
                    }
                    return profile;
                })
                .collect(
                        ProfileList::new,
                        ProfileList::add,
                        ProfileList::addAll
                );
    }

    public Set<Sequence> motifEnumeration (int k, int d){
        var Dna = this;
        if (Dna.isEmpty()) return new HashSet<>();
        var first_dna = Dna.get(0);
        if (Dna.size() == 1) return first_dna.kmerSequenceStream(k).collect(Collectors.toSet());
        Dna.remove(0);
        return first_dna
                .kmerSequenceStream(k)
                .parallel()
                .flatMap(seq -> seq.getNeighbors(d).stream())
                .filter(
                        seq -> Dna.stream().allMatch(
                                dna -> dna.kmerSequenceStream(k)
                                        .anyMatch(
                                                s -> s.hammingDistance(seq) <= d
                                        )
                        )
                ).collect(Collectors.toSet());
    }

//    public int score(){
//        Function<Map,Integer>
//
//        return this.countMotifsToStream()
//                .map(hash -> hash.keySet().)
//                .peek(i -> System.out.println(i))
//                .reduce(Integer::sum).get();
//    }

    public GenomeList greedyMotifSearch(int k) {

        var bestMotifs = new GenomeList();
        var kmers = get(0).kmerSequenceStream(k).collect(Collectors.toList());
        for (var kmer : kmers) {
            var motif0 = new FuzzyGenome(kmer,0,kmer.length());
            var motifs = new GenomeList();
            motifs.add(motif0);

            for (int i = 1; i < this.size(); i++) {
                var profile = motifs.computeProfile();
                var mostPKMer = get(i).profileMostProrableKmer(k,profile);
                var motifI = new FuzzyGenome(mostPKMer,0,mostPKMer.length());
                motifs.add(motifI);
            }
            if (computeEntropy(motifs.computeProfile().getAllValuesStream()) >
                    computeEntropy(bestMotifs.computeProfile().getAllValuesStream()))
                bestMotifs = motifs;
        }
        return bestMotifs;
    }

    public static void main(String[] args) {
        String[] dna = (
                "GGGTAGAAGTGTTGGATCATTAGGGTGATGCGGGTGGGTATATAAATCCGCTAGGGGCTCGTACATGCATAAGGAGCGACGCCCCTGTTGCTTAAGGATAAATCCTCTCTACCCCGTCAGTGCTTCGTAGACGCGACTGACAGCATCGCCGCCCGA\n" +
                        "ATTGCGAAATCCCTAACGATGATACGTATAATCTTCCAGCTATATCACTATGTAAACAAACGATAGCAAGGCCTTGTTCTTGACTTACCGTGAGCGGATCAATAAATCATAAACAGAGTCTCTGCATTGGTACCGATCGCTCTGATATAACCCATG\n" +
                        "ATCACGCTAGATCTTTCGTTGTACGGTGAACTGTTGACAACGTATTTAGATACAAAAATCGTTTACCCGTAGGATTCTGGACCACCGGTAGGCTCCCGAACCCTCCATCGGGCACTAGACTTTTGACTTGTAGGTCAATCCTGCACCATCCATGTT\n" +
                        "CACAAATTCACGACGGACTACACACGATCATCGTTCGACCGCATGACGATCTTCGGCACCTGGTAGAAGACATATAACCCTGACGGTTGACAAATCAACCACAATCGACGAACCACTACCTAGCACAACTGATTACGAAAGCCCACATGGAGTCTC\n" +
                        "CTGGAAGGCACGGTTAGGCCCCCTATCTTACAATTAAAAAACAGTCACACAGGTTGTAGGTGCCGTCGCGCATCCAATGAATTTGACGCAGTAAGGTCCACCAAACTAGGTTGGTTTGTAGATGAACAAATCGGCTAATGTTAATAAAGATGTTGG\n" +
                        "CACTGGCTCCTTACCTGTCCTTTTCTTGGCAGCTTGTCCGGGTCCAAACCTTTTTTTGGAGCTATAAAAATCATGATGTTCCGGATCATCTCCCCTCTGGACTAGCGCGTGTGATGTACAGCATTACGGGCCCCTTACCTAAGCTGCGGCTGCTCT\n" +
                        "TAAAACCTGATTCGCAAAGCCATTGTGTTGGATAGGAGGCTACACCCTTGGGCGTAACAGACTTTGTCCGAACAGACAAGCGCTGTTGAAGAAATCGCAATGTACACGTGTCATAATAAATCCATCGCGTGAACTCTCGACACTACATATAGCCCC\n" +
                        "AGAGGGGTATAATTACTTGCGCCGTTTATCACTTTCTATAGACATTGCACTTCCCCCCGAATGGTACTTCTGGTAATGGTTTCTCTGTACCGGAGACGTACAACTGCTCTTAAGGCGTGCCTGATTTTGCAAGATATAAAAATCCGGGGAGAGGCC\n" +
                        "TATTCAGTACGCTATCGAGGGTCCCGAAGCTAACCGTCTATGTAAAGAGCTCGTCCATGTAAGCTGATAGATGCCTAATCTCCAGGAGAACTCGTCGCCCCATAAATATGAAAAACGTTGGCTTTACAAATCGTCTGATACATACTAAGTAGTATA\n" +
                        "CTATATTATTGTAGACTTTGCCTCAAATCTCTCTACGTTTGAAAAATCAGGTTGGAGCATGGTGGTATGATATCGGTAACAGCTTGCGATGGTCGTGCGGCGGGACGGATATGACGCTATAGTAAGTCCTATAGACGCTATAACTTAGTGTACGGT\n" +
                        "ACTTATTTAAAGCCACCACCGCACGGTATAGAAATCCTTAATTGGATAAGACTTTGATGCACCTCACATCTACGCGCGTGTTGATTGTCGTCCCGGAAGGGAACCATCGGTCGGTACAGATTAGGATGAATTCTCTTCACGCATCGTTTCTCACAG\n" +
                        "TCGTCGACCCATTGTTCTTATAGTGCTGGACAAATCATAATGACGAAACGTGTGAAAGAGCATGGCTGTTACTAAGGGGGCATGTAAATGATTGTTGCTCCGTCTACCAGAACTTTATTAACTCTAGGGCTCCCACTATGAAAAAGTGGGACCCCA\n" +
                        "GTTCGAGAAATCTGGATTTTCATTGTATTTTTGTATCAAACTCATCCGGTACCCTGTCGTCGAAGCGGAAAAAAACGATCCAGAAGTTGGTGAGAAGCCATTGTCATATATGCTGCATTAGGGTGAAGATATGGGTTCAGGGGGGGCGCAGCTGCC\n" +
                        "GCTAGAGAAATCGTTCCAAGCACGGGGGTATTACTCTAATGCAATAAACGGCGAAACACAGCACAGACATCCGCTTCAAGCACACTCGTAACTGGGCCTGGTGAGGCATGTCGAGTAAACTTTCACGCCTCACTTGTGACTGGGCGCACAATACTT\n" +
                        "GAGATAGGCGGGGTTTCACAAATCCCGTCTTATCCATAGCGAGACTTCCACGTCTGATTCACGTCAATATACGAATCTTAGCTACCCTTGTCCCCCGAGGAAAAGTGTATCTCGCGTCGCATTTGGCCTTCTAACATTCAAACTCCCGGGTTACTG\n" +
                        "CTCATGCCCTGCTGTCTAATGCCTGCGGAGGTTATTGATCCACAAATCGTGAACGTGGGTACGTTAAAGAGTCCCGACATTCCTCTTTGGTTCATGCAGTTTCCTTCGCGCCACTATACGGAAGCTACTTCTCTTTGGACGACTTCATCACAGCGG\n" +
                        "GACAAGAACCGACAACGCACCACACACATGTTGAGAAAGGGCATGCTTTGCGTTGTCAGCTGCCGCGGGGGTAGAGGCAAATTCGATCAAGAAATCTAGTTAATAACCCTGCTCGGCAGCTACACTGCATCAGCCGGACCATAGCCTTTAAGCCCT\n" +
                        "GCGGATAGTTATAATGCGCCGACAGAGCTAGCGCTTGCTCGCAATTTTACAACATGTCGACTCTATATAACGGAAAGGTTCTTGTTATGATGCTGACTTTCGCCACCATCGTAACCAACCAAGCGGCCGAGTGTAAAGCCGTCGGTTGAATAAATC\n" +
                        "CTAGGGGCCCCGTTGAAAAAGGGCCTGGATATAAACAATTGACTACTCTCTTTAGGAATTTTCTGGAAAAAAGTTGGACGGCCATGGCGTCAGTCAAACACCTTAACCGTTTTATAAATCTGGCCACGAAGGACGCTTAACGGTATATTACGATAG\n" +
                        "GAGGCAGGTAGTTGCATGGCGTGGATCGGTTCCTGAGAACACTATGGCAAAATTCCACTCGGATTCTACTGGTCTTGAACGACGGAGGGGCTCGGGGAGCGTGAAAGCTGTGAAGATTCCTACTCGGGCCAGGGTTGACAAATCGACAACAATAAC\n" +
                        "CTTTTCATTCCCGCTGCCAGGAGTCGGATCAGAAAAGCTATGATCTCGCTCACAGGCAAAAGGACATGTTACATGCGGGGCGCTGCTCGACAAATCGTCTGATGTCGAGACCCCCACGACTCCGTACTAGAGAACCGGAGCATTGCGGCCCGGACT\n" +
                        "AACAGCACGCAACATTATACCACCGGCGTGAGCCTTGAAGACAAACGTAACCGACTTGAGATTTCGTAGATATGAGCGATCTTAAACATGGAGCTTGATTTAGAAATCCGCAACGGACTTCACCTGTGCCAAGGAACACTTACGCTCTCGTTCGCA\n" +
                        "GTGGCTATTCCCAATGGGGCAAGGCCCGGAAGTGCTCTGTTTCCGACATAGCATACTCCCTGGCCGAGTTAAAAGTGTTTGGTCGATACAGAAATCGGATGCTCCTATCTGCTTGGCTGAGCAATGCCCGCTTTAAAGCAAACTGTTCGTCTCGAT\n" +
                        "AGCTCGGACTCCAGAGCTCCCCGTTCAGAGTGGGCAGAGGAGAACGGCATTAGGTATCGTTGGGCGTGGAAATAATGATACGTGACGCCAGTCGCGATCCTAGCCACCTTGATCGTGGGACAGTGCACCCCAATGGTTCGCTTGGATCCAGAAATC\n" +
                        "TCGCGAACGTTCCAAGTGGGAGCCTTGGGATACACGGTGTCCTACGATGCACTCACGTAAGTATTGGGCTTAGCTCGACAAATCTGCGTTAAACTTATCGAGTGCGATATCGCATCCTTCTCAAAGGAGAAGACTCTGCGCGATGCGTACATGCGG\n").split("\n");
        var result = new GenomeList(dna).greedyMotifSearch(12);
        result.forEach(
                m -> System.out.print(m + " ")
        );
        System.out.println();
    }


}
