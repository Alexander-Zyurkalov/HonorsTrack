package org.stepic.bioinformatics;

import java.time.Duration;
import java.time.LocalDateTime;
import java.util.*;
import java.util.function.BinaryOperator;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class FuzzyGenome extends Genome {
    public FuzzyGenome(final String text) {
        super(text);
    }

    public FuzzyGenome(Sequence seq, int beginning, int length){
        super(seq,beginning,length);
    }

    public List<Integer> findAllPositionsOfTheApproximatePattern(final Sequence pattern, int searchingDifference) {
        int k = pattern.length();
        var reversedPattern = pattern.reverseComplement();
        return Stream
                .iterate(0,i->i+1).limit(this.length() - k + 1)
                .filter( i -> {
                    var seq = new Sequence(this, i, pattern.length());
                    return seq.hammingDistance(pattern) <= searchingDifference ||
                            seq.hammingDistance(reversedPattern) <= searchingDifference;
                })
                .collect(Collectors.toList());
    }

    public int approximatePatternCount(Sequence pattern, int difference){
        return (int)kmerSequenceStream(pattern.length())
                .filter(seq -> seq.hammingDistance(pattern) <= difference)
                .count();
    }

    private Stream<Sequence> computingFrequenciesWithMismatchesStream(int k, int d) {
        return Stream.iterate(0, i-> i + 1).limit(length() - k + 1)
                .parallel()
                .map( i -> new Sequence(this,i,k))
                .flatMap( pattern -> Stream.iterate(pattern,p->p.reverseComplement()).limit(2))
                .flatMap( pattern -> pattern.getNeighbors(d).stream() )                ;
    }

    public Map<Sequence,Long> computingFrequenciesWithMismatches(int k, int d) {
        return computingFrequenciesWithMismatchesStream(k,d)
            .collect(Collectors.groupingBy(
                    pattern -> pattern,
                    Collectors.counting()
            ));
    }


    private static Set<Sequence> frequentWords(Map<Sequence,Long> frequentArray) {
        Sequence max = frequentArray.keySet().stream()
                .max(Comparator.comparingLong(frequentArray::get))
                .orElse(new Sequence(""));
        return frequentArray.keySet().stream()
                .filter(pattern ->
                        frequentArray.get(pattern) == frequentArray.get(max))
                .collect(Collectors.toSet());
    }
    public Set<Sequence> frequentWordsWithMismatch(int k, int d) {
        Map<Sequence,Long> frequentArray = computingFrequenciesWithMismatches(k,d);
        return  frequentWords(frequentArray);
    }

    public static Set<Sequence> motifEnumeration (List<FuzzyGenome> Dna, int k, int d){
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

    public Sequence findMostProbableKmer(final int k,
                                         final List<Double> profileA,
                                         final List<Double> profileC,
                                         final List<Double> profileG,
                                         final List<Double> profileT){
        final var profile = new HashMap<Character,List<Double>>();
        profile.put('A',profileA);
        profile.put('C',profileC);
        profile.put('G',profileG);
        profile.put('T',profileT);

        return this
                .kmerSequenceStream(k)
                .max(Comparator.comparingDouble(
                        kmer -> Probability.probabilityByProfile(kmer.getText(), profile))
                ).get();
    }

    public static void main(String[] args) {
        var start = LocalDateTime.now();
        System.out.println(start);

        var genome = new FuzzyGenome(
                "GAGCACCACGGAGCTAGAGGGCGAAACTGGGCAGTAAGTTCAAATGTGTTACCAGGACAGATGATTACGACAATGCGACGTTATAGTGTGTCCCTGAACACCGATCGGTAATTCGCGTCTGAACAGTGAATGTAATTCTCCTCATTTTCACGGCGTTATAATACTTCCCTTACCAACGCTCTCTTGCTAGCGATCGTTGTGAGAAACCGGGCATGATGAACGTCTAAGTGTCCTTGCGCAAAGACCCCAGCCGCTTCCGGCATCTAAAGGTAGGGTGGAAATCTCGCTGAAACTAAAAAGCCAAACCTTAGGTGATGCGCATGGCCTGGGTTTCCCCGTCTTTTCATATCTGTTTAGCAGGTTGTAGGATATGGACACATTACTGTCAATTTCGCTCTTTTTTGGCGACTTGTTGGCGTCATGGACTAATCTGCTCCTTTCTACAATGGCCGTAACGTCTCGAATGGCCTGTATTCCATGATTCTAGTGCCGGACGGGTGTGAACACTGAGTACGCCCGATTTACCGCTGGTTGTCTACTCACCGGCCGGGTCTTCATCTTGGGAACATAATCGATTATAGAGCGAGAGTCAACGAGTATCTACACAACAGATATTAGAGAAACACTTCAATGCAGACTGACCAAATCCTACTTTGTTGTCGAGGCCTGGTATAGCTGTTAATGGAAATCCGTTTGACCGCTTTTCACCGATACTGGTCCACGTGGGAATCTTTTGAAGCCCACCTTTATTGACCGACACTCTGAAGACTCGAATAGAGGAAGTGTATGTTAGCATCATTTACCGGAAAGGCAGTCGCTGCCCTGACCGTTCGCCCGTAATGCGATAATGCGCAAAAAGCCCATCTACTCCTGGAATCGCTCCGCCTTGATATTCGAAGTCTTCTTTTCGAAATCCATTAATGGAAAGAGGAGGCACAATTCTGCTATCCGGCTCGCGGCCAATCGGGTAAGGATACCGCCTCACGTCTAGGGGAAACGC");
        var result = genome.findMostProbableKmer(14,
                Arrays.asList(0.296, 0.31, 0.169, 0.127, 0.282, 0.31, 0.268, 0.225, 0.282, 0.211, 0.239, 0.268, 0.352, 0.225),
                Arrays.asList(0.127, 0.183, 0.296, 0.296, 0.268, 0.183, 0.155, 0.352, 0.183, 0.31, 0.296, 0.239, 0.239, 0.324),
                Arrays.asList(0.254, 0.282, 0.268, 0.254, 0.127, 0.211, 0.296, 0.225, 0.282, 0.239, 0.324, 0.225, 0.211, 0.225),
                Arrays.asList(0.324, 0.225, 0.268, 0.324, 0.324, 0.296, 0.282, 0.197, 0.254, 0.239, 0.141, 0.268, 0.197, 0.225)
        ).toString();
        System.out.println(result);

        var stop = LocalDateTime.now();
        System.out.println(stop);
        var duration = Duration.between(start,stop);
        System.out.println(duration.getSeconds());

    }
}
