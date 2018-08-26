package org.stepic.bioinformatics;

import java.util.Comparator;
import java.util.Map;
import java.util.Random;
import java.util.stream.Stream;

public class GenomeListWithPseudocounts extends GenomeList {

    public GenomeListWithPseudocounts newGenomeList(){
        return new GenomeListWithPseudocounts();
    }
    public GenomeListWithPseudocounts newGenomeList(final String... texts){
        return new GenomeListWithPseudocounts(texts);
    }
    public GenomeListWithPseudocounts newGenomeList(FuzzyGenome genome){
        return new GenomeListWithPseudocounts(genome);
    }
    public GenomeListWithPseudocounts(final String... texts){
        super(texts);
    }
    public GenomeListWithPseudocounts(FuzzyGenome genome){
        super(genome);
    }
    public GenomeListWithPseudocounts(){
        super();
    }
    @Override
    Stream<Map<Character, Long>> countMotifsToStream() {
        return super.countMotifsToStream()
                .map((Map<Character,Long> map) ->{
                    for (char key :
                            new char[]{'A', 'C', 'G', 'T'}) {
                        map.put(key, map.getOrDefault(key, (long)0) + 1);
                    }
                    return map;
                });
    }

    @VisibleForTesting
    GenomeListWithPseudocounts getRandomMotifs(int k){
        var random = new Random();
        return this.stream()
                .map( dna ->
                        new FuzzyGenome(dna,
                            random.nextInt(dna.length() - k + 1 ),
                            k)
                )
                .collect(GenomeListWithPseudocounts::new,
                        GenomeListWithPseudocounts::add,
                        GenomeListWithPseudocounts::addAll);
    }


    @VisibleForTesting
    GenomeListWithPseudocounts getMotifsFromProfile(ProfileList profile) {
        var k = profile.size();
        return this.stream()
                .map( str -> str.profileMostProrableKmer(k,profile) )
                .map(seq -> new FuzzyGenome(seq,0,seq.length()))
                .collect(GenomeListWithPseudocounts::new,
                        GenomeListWithPseudocounts::add,
                        GenomeListWithPseudocounts::addAll);
    }
    public GenomeListWithPseudocounts randomizedMotifSearch(int k, int howManyTimes) {
        return Stream.iterate(
                this.randomizedMotifSearch(k),
                motifs -> this.randomizedMotifSearch(k))
            .limit(howManyTimes)
            .min(Comparator.comparingInt(GenomeListWithPseudocounts::score)).get();
    }
    public GenomeListWithPseudocounts randomizedMotifSearch(int k) {
        var motifs = this.getRandomMotifs(k);
        var bestMotifs = motifs;
        while (true) {
            motifs = this.getMotifsFromProfile(motifs.computeProfile());
            if (motifs.score() < bestMotifs.score()) {
                bestMotifs = motifs;
            } else {
                return bestMotifs;
            }
        }

    }

    public static void main(String[] args) {

        var list = new GenomeListWithPseudocounts(
                ("TCAACGCGCACATTGTAGGAGCTGTAACCCGGTACAGGGATTCGTACGCCCTGGCCTCAAGTATTTCCACGGGGAGGGCGGACGTTCCTTATTGAGGCGCTAGGCGCGGCGGGGCAACATGCGCATGGGTTAGCACGGGGATTCTAAGTGGAATGCAACCGCGGGGGATCGCACAACACGGGGAAATCAACGCGCACATTG\n" +
                        "TAGGAGCTGTAACCCGGTACAGGGATTCGTACGCCCTGGCCTCAAGTATTAACTAGAACTACAGCTCCACGGGGAGGGCGGACGTTCCTTATTGAGGCGCTAGGCGCGGCGGGGCAACATGCGCATGGGTTAGCACGGGGATTCTAAGTGGAATGCAACCGCGGGGGATCGCACAACACGGGGAAATCAACGCGCACATTG\n" +
                        "ATATCGGATGGTATGGTGGCTTAGCACGTAGACGGGAACATGCATGGAAGCACGTGACCAAACGAATGACACATCAACTAAATAAATTGCAGAAGTATGTATATGCATACCTTAACTACACAGTATGAGCCACGTAGTTGAGAAGGTCTGCTCTTCAGGAGTCAGCCGACACGTTCTATCAAGAGCTTACGCATGCCAGTC\n" +
                        "TCGTCTACACTGGGCAACTCCGGCATAACAACTGACGAGTGCTGCGGGGTGTCACCGAAGCCATCAGTTTGCGCTTCCCGATTCCCATAACCTGGCCTACAGCAAATACAAAAAGTGTCGACGAGTCCACAAAGTTATGATTCCAGTCTGGTCGAAATCGATCTCGCTCGCGTTATAAGTACAGCGGTTGCACCCGCAAAC\n" +
                        "ACACTAAAGGGCTGGTATCTCTAGGAAAAGGCGTGTCCATCATTGTCTCGGACGGATGTGCCTGTTGCTATCCCCAGCAACCTTAACTAGTACAGCCAAACCGGCCTTGAACAGGCGCCGTTCGATCACAGTATTTTCAGTACGGGGACCCCTTCCACGAGTGGTCGAATGTTGTCATTGCCCTGGTAGGAGCGTTGCTCT\n" +
                        "TGTCTGTTCCTACCGCGCGGCTATGCGAGTGCCTCGGATCGGCCCCAGCGTCAAAGTCATTAACGACGGGAACCATTCGATGTCGTCTCCAGTAGCGCACTTGTGTGCCTAATCTCTAGACATGAAACCGATCTGGATAACACAACCTTAACTCGCGCACCGACCCGAATACTGTCGTTCAAATAGTGTTTCAAAGATGAC\n" +
                        "GAGCTTCCTTTATCCGGGGGTTTAAGACATATTCATCTAGTCACCTCTCGCTCAGTGGGTACTTCAACCCTCGTCGACACTTCGCGATGACTCGCCCAAACTTCGGGAAGCCAGTGCATGAACCGATGCTTCTAGCCTCAATTAGAGGATTATTCCCTTAACTACAGGGTTTACCCTAACAATCCAGTAATTTTTTCTTAC\n" +
                        "ATTTGAATGAGGGTTTTTCTCCGACCTTGAAACATTGCAAAGCGTCAAAGCAGACCCAGTCACTCACATGGTTCGTTCCAGGAGCACTCCCTTGACCGATGCGCTGTACATATTAAAGTCTAAACCGCTACTACAGCTGAGCAAGCGCAGCTTAACAGCGTTCCAGCAGAATTTGACGCCATATCGTAACAAGCCTAAACG\n" +
                        "CATGAGCGATAACACAACACCACACCCGGTAACCCCGATTCTGCTTTGGTCTATCACCAAGGGTGCCTCACGGGACCCATATTACCTTCCTGGCTACGCGACTAACTAGACTTGCAGCCAACCTTAGTAACAGCGCGAGAGACCGGGAGGTAAGCTGGCCTAACTCGGTCCTTCAGATAGTGGGCTGGATACCCGAGAATG\n" +
                        "ACCAAAAGCGAGAAGCCCGCTCTGGTCGAAAGGTGGAACCCCTTTTTGTCGTTTGTCAGATAAGTGTAGGACCCTAACTAAGGGCGAGCCACGTCCCCATTAGCGGAGACATCGGACCCACTTAACTACAGCTTCAGAGCTCAGTCGTAGGATATGTCGTTGGGGCAGCACTATGTTTGAAGAAGTAAAAAACCTACGCGC\n" +
                        "TTATATCCAATGGGTGCGGAGAATTGGTTATGTGGGTCAGAATGGCTCGGCAGTCTATTATTATGAGGGGCCTCGGGTTTTGGGTATGCTCGCATAAGTCGTTTCTGGAGAGTAATCCCCTACCGGGAATTAGGCCAACCTGAACCCTCCGTCCGAAGTAGTCGTAAGTCCCAACCTTAAACCCAGCACCTTGTCCCGGCT\n" +
                        "CATGGCGATCCACCACTTTACGGAATGATTCTGGGCACACCCTAGATCGAACCTTGTATACAGCGCGGGTCCAGGATACTGTCGTTCCTGGATAAGTTAGCTTCTCTCGAATGGTAGTACGCCCCTGTAGTGACCAGGAATTTTAACTCGCATTGGTAAGGCTTCTCGAGGGCGCGGCACCGAACGGAGTACTATAAATGG\n" +
                        "ACCTTCAGCGATGTGCGCGTTGGAAACTTCTCTAGCAGTGGGGCGACCCCGTGCCAATCATCTTTTGGCAGCTGACGAACCAACCAAGGGCGGGGGGGGCGACTGGTTATCCGGGGGAGAGTTTACCGCACAAATACTCGAAACCTTTCTTACAGCGTCTCCAAGAAGTGACCGCACTTACGGCGAAGTACATATGCCGCC\n" +
                        "CTTCAGGTTCGTCCCCATCCGGGGCGTACTAAAATCCTGACTATTTAGTGTCCTGCCTGGGACTTAGCCTAAGCCGGCCGGCACCCACGCCGTCCATAACCTGGCCGACACTTGAGGACTGTAACGCTTCCGGGTGGGGTCACTGTGCGGCGCCAAAAATTACTCCGCCTGAAACCAGCACTACAGCAAGCTTGACAGTTT\n" +
                        "ATTCTCACAGATGCATACGTCAATCTGCGGTTTGCAGGGTGTTTCTCGCGTTTCGGGGTCTCGAAGCGTTGCCAAGGTCCGTCTACTGTGCGGGTTTCAATAACCTTAACGCGAGCCGTTTTAGCGGGTGCAAGACCATGAGGGTAAGACATCTATGACCGGCTACCCGCAAATCACGGTTAACATCATAGCTCCAAGGAG\n" +
                        "CAATAACAGAGATTAATTGTCGAAAAACGTCCAGTTGAATCGCCCTATTCCCGAAATCAAACACAGGCACCGCGCCAATCATAATGATTTACGGGCTTAGGTTCTAGGTCGATGCGAGTCGATCTTTTTATCAATGAAAGATGAGTCAACCAGCCCAACCGAAAACTACAACTACAGCGCTTTCTTCAAAACATTTTTTCC\n" +
                        "TGCAAGGTCCGATACATCCCAGGGACCATAAAGCTCAAAATGGACACAACCTATGCTACAGCCGAGGGGCTGCTTCACTGGAAATTGGCTTTACTTAGCAAGAATCGACGTCAGTGGATACAGAGCATAATATGCACCGCACCGCAGCCTCACTGTTCCTCTCGTCTGCAGGATAGCTGTCCTATATACGGGTAGTTAATC\n" +
                        "AGATAGTCAAGAGTCGTACCAGGAGTTGTTATAGCTTCGAACGACTTTCTCAGCCTCTATGCACATAACCTTAACTACTTACTAACAAATCCCTGGAGGCGTACGCTTGGGTTAAAAGAGGTCGATCCGGGACGCCAACTCGCGGACATCCACCTACCGCTTCGCCTCGAAGCTGACGAGGGCGCGCCGTTGGCGTGGGTC\n" +
                        "GACGGAGCGGGGCGGATCCCCCGCTATGTTGGCAGAGGTCCTTATGTGTTCCGAACTCAGGTTTGCCCCTGAATTATCTGCAAAGCTCTGTTCACGTCGTCTGGCTAGGCTCAAGAATGCTGTCGGTGAGAAGTCTAACTACAGCCCATATTCCCAAAGACCGGGGAGAAGGATTTGGCGGACAATCGATACTTACTACCG\n" +
                        "TGAAAAAGCCAGCGTGCTGATAATCTACGGGTTACGGATTTATCGCCGTTTTCCACCCTTTTAGTAATGCTGTTACTGGGATATGGACGGTTATGAATGATTCGGATTTTTAACTACAGCGTGCCGTTGACTGACAGCTTTTCTAGTGTAAGTGCTACTGGGCGAGTAATGTCAGATTGATTATCGGTGACCTACATTATT").split("\n")
        );
        var result = list.greedyMotifSearch(15).stream().map(FuzzyGenome::getText);//.sorted().collect(Collectors.toList());
        result.forEach(
                g -> System.out.print(g + " ")
        );
        System.out.println();
    }
}
