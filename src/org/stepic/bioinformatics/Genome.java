package org.stepic.bioinformatics;

import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class Genome extends Sequence {


    public Genome(final String text) {
        super(text);
    }
    public Genome(Sequence seq, int beginning, int length){
        super(seq,beginning,length);
    }

    private final ConcurrentHashMap <Sequence,Integer> patternCountCache = new ConcurrentHashMap<>(length());

    public Stream<Sequence> kmerSequenceStream(int k){
        return Stream
                .iterate(0,i->i+1).limit(this.length() - k + 1)
                .map( i -> new Sequence(this, i, k) );
    }

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
        return kmerSequenceStream(k)
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
        Set<Sequence> set = kmerSequenceStream(k)
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
        HashMap<Sequence,Integer> map = kmerSequenceStream(k)
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
        HashMap<Sequence,Integer> map = kmerSequenceStream(k)
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
    public List<Integer> findAllPositionsOfThe(final Sequence pattern) {
        int k = pattern.length();
        return Stream
            .iterate(0,i->i+1).limit(this.length() - k + 1)
            .filter( i -> this.equalsAt(pattern, i))
            .collect(Collectors.toList());
    }


    public Set<Sequence> clumpFinding(int k, int l, int t) {
        var k_pow = pow.apply(k);
        var frequentPatterns = new HashSet<Sequence>();
        var clump = new int[k_pow];

        var text = new Genome(this, 0, l);
        int[] frequencyArray = text.computingFrequenciesImperatively(k);
        for (int i = 0; i < k_pow; i++) {
            if (frequencyArray[i] >= t){
                clump[i] = 1;
            }
        }
        for (int i = 1; i < length() - l + 1; i++) {
            var firstPattern = new Genome(this, i - 1, k);
            int index = firstPattern.patternToNumber();
            frequencyArray[index]--;

            var lastPattern = new Genome(this, i + l - k, k);
            index = lastPattern.patternToNumber();
            frequencyArray[index]++;
            if (frequencyArray[index] >= t) {
                clump[index] = 1;
            }

        }
        for (int i = 0; i < k_pow; i++) {
            if (clump[i] == 1) {
                frequentPatterns.add(numberToPattern(i,k));
            }
        }
        return frequentPatterns;
    }


    public Set<Sequence> clumpFindingHashes(int k, int l, int t) {
        var pow_k = k <= 14 ? pow.apply(k) : Integer.MAX_VALUE;
        var frequentPatterns = new HashSet<Sequence>();

        var text = new Genome(this, 0, l);

        var frequencyArray = text.computingFrequenciesHashes(k);
        frequencyArray.forEach(
                (Sequence key,Integer value) -> {
                    if (value >= t ) frequentPatterns.add(key);
                }
        );

        for (int i = 1; i < length() - l + 1; i++) {
            var firstPattern = new Genome(this, i - 1, k);
            frequencyArray.put(firstPattern,frequencyArray.get(firstPattern) - 1);
            var lastPattern = new Genome(this, i + l - k, k);
            if (frequencyArray.containsKey(lastPattern))
                frequencyArray.put(lastPattern,frequencyArray.get(lastPattern) + 1);
            else
                frequencyArray.put(lastPattern, 1);
            if (frequencyArray.get(lastPattern) >= t) frequentPatterns.add(lastPattern);
        }

        return frequentPatterns;
    }

    public int[] computingFrequenciesImperatively(int k) {
        var result = new int[pow.apply(k)];
        for (int i = 0; i <= length() - k; i++) {
            var pattern = new Sequence(this,i,k);
            result[pattern.patternToNumber()]++;
        }
        return result;
    }

    public Map<Sequence,Integer>  computingFrequenciesHashes(int k) {
        var result = new HashMap<Sequence,Integer>();
        for (int i = 0; i <= length() - k; i++) {
            var pattern = new Sequence(this,i,k);
            if (result.containsKey(pattern))
                result.put(pattern,result.get(pattern) + 1);
            else
                result.put(pattern,1);
        }
        return result;
    }

    public Map<Sequence,Integer> computingFrequenciesFunctionally(int k){
        var pow_k = k <= 14 ? pow.apply(k) : Integer.MAX_VALUE;
        return kmerSequenceStream(k)
                .collect(Collectors.toMap(
                        (Sequence seq) -> seq,
                        (Sequence seq ) -> 1,
                        (value,new_value) -> value + 1,
                        () -> new HashMap<>(pow_k)
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

    public List<Integer> minimumSkewImp(){
        var result = new ArrayList<Integer>(length()+1);
        int prev_value = 0;
        int new_value = 0;
        int min = length();
        for (int i = 0; i < length(); i++) {
            if (charAt(i) == 'G')
                new_value = prev_value + 1;
            else if (charAt(i) == 'C')
                new_value = prev_value - 1;
            else
                new_value = prev_value;
            if (min > new_value) {
                min = new_value;
                result.clear();
                result.add(i + 1);
            }
            else if (min == new_value)
                result.add(i + 1 );
            prev_value = new_value;
        }
        return result;
    }


    public List<Integer> skew() {
        return Stream.iterate(0, i -> i + 1).limit(length())
                .map(i -> charAt(i))
                .map(ch -> {
                    switch (ch) {
                        case 'G':
                            return 1;
                        case 'C':
                            return -1;
                        default:
                            return 0;
                    }
                })
                .collect(
                        () -> {
                            var list = new ArrayList<Integer>(length() + 1);
                            list.add(0);
                            return list;
                        },
                        (ArrayList<Integer> list, Integer value) -> list.add(list.get(list.size() - 1) + value),
                        ArrayList::addAll
                );
    }


    public static void main(String[] args) {
        System.out.println(
                new Genome("CATTCCAGTACTTCATGATGGCGTGAAGA")
                .skew()
        );
    }

}
