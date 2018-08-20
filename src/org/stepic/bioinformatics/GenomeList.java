package org.stepic.bioinformatics;

import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static java.util.stream.Collectors.toList;

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

    public List<Map<Character,Double>> computeProfile() {
        var length = this.size() > 0 ? this.get(0).length() : 0;
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
                .collect(toList());
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

//    public

//    public List<Sequence> greedyMotifSearch(int k, int t) {
//        var Dna = this;
//        if (Dna.isEmpty()) return new ArrayList<>();
//        var first_dna = Dna.get(0);
//        if (Dna.size() == 1) return first_dna.kmerSequenceStream(k).collect(Collectors.toList());
//        return first_dna
//                .kmerSequenceStream(k)
//                .parallel()
//
//
//
//    }

    public static void main(String[] args) {
        String[] motifs = (
                "TAGCTATGTTAGTTACCACCTGATATATTAGTTCGACCCGAG\n" +
                "GGTACAGCGTGCCCTGAGTAAGGCGTCGTGTAAACCGCTTGC\n" +
                "TTTGCGCGATTCGTAGTCCCGGAGTGGCCGCGCAATTGGTAA\n" +
                "CCTGAGCCTGTGCGCTGTAAGGCAAGACGTTTCGATTGTCAG\n" +
                "AGAGCCTGGAGGCGCTTGCCAGAGTGTCCGTCTCTAGCGGTC\n" +
                "CACCAACCGGAGCTAGCGAGTCAGCCTCACGTTATGTCGCCC\n" +
                "GATTGTCCAGAGTTAAGAGTTGCCCACGCCTTGATTTGCTTT\n" +
                "AGCATTCCGGAGAAGGGCATCGATAGTGACTCTTTGGGTCCC\n" +
                "GTGGGAAACCTGCAGTCGTCCGTGGGCCCTAACGTGCCAGAG\n" +
                "CCTGAGGGTCGTTGCAGGGCCGACCTGCATACGCAGTTGTAC").split("\n");
        System.out.println(new GenomeList(motifs).medianString(6));
    }


}
