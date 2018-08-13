package org.stepic.bioinformatics;

import java.util.*;
import java.util.function.Function;
import java.util.stream.Collector;
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
        return (int)Stream.iterate(0,i->i+1).limit(this.length() - pattern.length() + 1)
                .map(i -> new Sequence(this,i,pattern.length()))
                .filter(seq -> seq.hammingDistance(pattern) <= difference)
                .count();
    }


    public Map<Sequence,Long> computingFrequenciesWithMismatches(int k, int d) {
        return Stream.iterate(0, i-> i + 1).limit(length() - k + 1)
            .parallel()
            .map( i -> new Sequence(this,i,k))
            .flatMap( pattern -> Stream.iterate(pattern,p->p.reverseComplement()).limit(2))
            .flatMap( pattern -> pattern.getNeighbors(d).stream() )
            .collect(Collectors.groupingBy(
                    pattern -> pattern,
                    Collectors.counting()
            ));
    }

    public Set<Sequence> frequentWordsWithMismatch(int k, int d) {
        Map<Sequence,Long> frequentArray = computingFrequenciesWithMismatches(k,d);
        Sequence max = frequentArray.keySet().stream()
                .max(Comparator.comparingLong(frequentArray::get))
                .orElse(new Sequence(""));
        return frequentArray.keySet().stream()
                .filter(pattern ->
                        frequentArray.get(pattern) == frequentArray.get(max))
                .collect(Collectors.toSet());
    }

    public static void main(String[] args) {
        var genome = new FuzzyGenome("CATGCCATTCGCATTGTCCCAGTGA");
        System.out.println(
                genome.approximatePatternCount(new Sequence("CCC"),2));

    }
}
