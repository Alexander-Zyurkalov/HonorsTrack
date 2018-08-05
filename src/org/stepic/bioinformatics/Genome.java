package org.stepic.bioinformatics;

import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

public class Genome extends Sequence {


    public Genome(final String text) {
        super(text);
    }
    public Genome(Sequence seq, int beginning, int length){
        super(seq,beginning,length);
    }

    private final ConcurrentHashMap <Sequence,Integer> patternCountCache = new ConcurrentHashMap<>(length());

    public int patternCount(final String pattern){
        return patternCount(new Sequence(pattern));
    }
    public int patternCount(final Sequence pattern) {
        Integer result = patternCountCache.get(pattern);
        if (result == null) {
            result = (int) Stream.iterate(0, i->i+1).limit(this.length() - pattern.length() + 1)
                .filter(
                    i -> this.equalsAt(pattern, i)
                )
                .count();
            patternCountCache.put(pattern,result);
        }
        return result;
    }

    public HashMap<Sequence,Integer> getKmers(final int k){
        return Stream
                .iterate(0,i->i+1).limit(this.length() - k + 1)
                .map( i -> new Sequence(this,i,i+k) )
                .collect(Collectors.toMap(
                        (Sequence seq) -> seq,
                        (Sequence seq) -> patternCount(seq),
                        (Integer old_value, Integer
                                new_value) -> old_value,
                        HashMap::new
                ));
    }
//
    public List<Sequence> frequentWords(final int k) {
        Set<Sequence> set = Stream
            .iterate(0,i->i+1).limit(this.length() - k + 1)
            .map( i -> new Sequence(this,i,k) )
            .collect(
                    HashSet::new,
                    (HashSet<Sequence> hashSet, Sequence seq) -> {
                        if (!hashSet.contains(seq)) {
                            int max = !hashSet.isEmpty() ? patternCount(hashSet.iterator().next()) : 0;
                            int kmer = patternCount(seq);
                            if ( kmer > max) {
                                hashSet.clear();
                                hashSet.add(seq);
                            }
                            else if ( kmer == max ) {
                                hashSet.add(seq);
                            }
                        }
                    },
                    HashSet::addAll
            );
        return new ArrayList<>(set);
    }
    public HashMap<Sequence,Integer> frequentWordsWithFreq(final int k) {
        HashMap<Sequence,Integer> map = Stream
            .iterate(0,i->i+1).limit(this.length() - k + 1)
            .map( i -> new Sequence(this,i,k) )
            .collect(
                    HashMap::new,
                    (HashMap<Sequence,Integer> hashMap, Sequence seq) -> {
                        if (!hashMap.containsKey(seq)) {
                            int max = !hashMap.isEmpty() ? patternCount(hashMap.entrySet().iterator().next().getKey()) : 0;
                            int kmer = patternCount(seq);
                            if ( kmer > max) {
                                hashMap.clear();
                                hashMap.put(seq,kmer);
                            }
                            else if ( kmer == max ) {
                                hashMap.put(seq,kmer);
                            }
                        }
                    },
                    HashMap::putAll
            );
        return map;
    }

    public HashMap<Sequence,Integer> fastFrequentWordsWithFreqFunctionally(final int k) {
        HashMap<Sequence,Integer> map = Stream
                .iterate(0,i->i+1).limit(this.length() - k + 1)
                .map( i -> new Sequence(this,i,k) )
                .collect(
                        HashMap::new,
                        (HashMap<Sequence,Integer> hashMap, Sequence seq) -> {
                            if (!hashMap.containsKey(seq)) {
                                int max = !hashMap.isEmpty() ? patternCount(hashMap.entrySet().iterator().next().getKey()) : 0;
                                int kmer = patternCount(seq);
                                if ( kmer > max) {
                                    hashMap.clear();
                                    hashMap.put(seq,kmer);
                                }
                                else if ( kmer == max ) {
                                    hashMap.put(seq,kmer);
                                }
                            }
                        },
                        HashMap::putAll
                );
        return map;
    }
//
    public List<Integer> findAllPositionOfThe(final Sequence pattern) {
        int k = pattern.length();
        return Stream
            .iterate(0,i->i+1).limit(this.length() - k + 1)
            .filter( i -> this.equalsAt(pattern, i))
            .collect(Collectors.toList());
    }

    public Set<Sequence> clumpFinding(int k, int l, int t) {
        return
                Stream.iterate(0,i->i+1).limit(this.length() - l + 1)
//                .parallel()
                .map(i-> new Genome(this,i,l))
                .map((Genome genome) -> genome.frequentWordsWithFreq(k))
                .filter((HashMap<Sequence,Integer> mers) ->
                        mers.size() >= 0 && mers.entrySet().iterator().next().getValue() >= t)
                .map((HashMap<Sequence,Integer> mers) -> mers.keySet())
                .collect(
                        HashSet::new,
                        HashSet::addAll,
                        HashSet::addAll
                );
    }


    public int[] computingFrequenciesImperatively(int k) {
        var result = new int[pow.apply(k)];
        for (int i = 0; i <= length() - k; i++) {
            var pattern = new Sequence(this,i,k);
            result[pattern.patternToNumber()]++;
        }
        return result;
    }
    public Map<Sequence,Integer> computingFrequenciesFunctionally(int k){
        return Stream.iterate(0,i->i+1).limit(length() - k + 1)
                .map(i-> new Sequence(this,i,k))
                .collect(Collectors.toMap(
                        (Sequence seq) -> seq,
                        (Sequence seq ) -> 1,
                        (value,new_value) -> value + 1,
                        () -> new HashMap<>(pow.apply(k))
                ));
    }

    public List<Integer> findPattern(Sequence pattern) {
        var result = new ArrayList<Integer>();
        for (int i = 0; i < length(); i++) {
            if (this.equalsAt(pattern,i)){
                result.add(i);
            }
        }
        return result;
    }

    public static void main(String[] args) {
        var text = "GGATCCTCTCGCCTGCCACCTGTCACTATCACTCGGAGGGTTAGGCCTTTAACCGAGGATAAGAAAGGTCCATACTTGTGACCATTATCGATCGTAGGACCACTGCACCAATTCTCACCATAACATAAGGCCAGCGGTGTCCCCGGCAGCTGAACTTGACCGGGCGCAGGCTTCCTAGTCCCGGGAACCGAGAGTCTACTTGCTGGCCTGGAGATGCGCTACCGCGCGACCTCACAGTCTTTCACCCCTATTAACAAAAACGCCCGTGATTGGCTAGAGAGCACGGTGGGTCTCAGATCCTGGGGGACCGGACCCGCGCTGGCTAATGGCGGTGGCCCCAGGGAGTCGGAAGAAAGGAAGCCGCGTGCTCGTTATTGACAAGTGATAGCCCGGTGGCATCCTAGGTCTATCCATTCAGTGTAGCCAACTTAATCGCATCAAGACCCTTTCACAAGCATCAGTGTCTAGCTGAAGAGATCGATGCAGGCGTATTGTGAGTTATAGCGCTTCGACAGCACGATTATGTCTTTTTATGCCTGTAAACTATATAGAGTGCGTGAGCGCCGAGTAGCCCCGGCCTGTATCGGAGGGTAATGCATGTAGTTCCCCTTAGCCGGCTGCCCGTGTTCAGCCGAATACCGTCGATGGCATAGAATATATGTGC";
//        int[] array = new Genome(text).computingFrequencies(5);
//        for (int i = 0; i < array.length; i++) {
//            System.out.print(array[i]);
//            System.out.print(' ');
//        }
//        System.out.println();
//        System.out.println(new Genome(text).clumpFinding());
//        var genome = new Genome("GACGATATACGACGATA");
//        System.out.println(genome.findPattern(new Sequence("ATA")));
//        System.out.println(genome.);

    }

}
