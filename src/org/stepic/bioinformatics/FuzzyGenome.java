package org.stepic.bioinformatics;

import java.util.*;
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
            .flatMap(seq -> seq.getNeighbors(d).stream())
//            .parallel()
            .filter(
                seq -> Dna.stream().allMatch(
                    dna -> dna.kmerSequenceStream(k)
                        .anyMatch(
                            s -> s.hammingDistance(seq) <= d
                        )
                )
            ).collect(Collectors.toSet());
    }

    public static void main(String[] args) {
        var list = new ArrayList<FuzzyGenome>();
        list.add(new FuzzyGenome("ACAGCGAATAAAGGAATGGGCATAC"));
        list.add(new FuzzyGenome("GTGTCAGAACGAACACCGCGCCTAA"));
        list.add(new FuzzyGenome("GTGTTGGACGCATGAGTAGATGAAG"));
        list.add(new FuzzyGenome("GCACAATTCTCATCTTTAGTTTAAT"));
        list.add(new FuzzyGenome("GTACAGTTCATTGCACGGGCTGAGT"));
        list.add(new FuzzyGenome("TCGTATGGAGCCTTTTTTCGGGAGA"));
        motifEnumeration(list,5,2).forEach(
            seq -> System.out.print(seq + " ")
        );
        System.out.println();

    }
}
