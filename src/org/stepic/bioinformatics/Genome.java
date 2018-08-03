package org.stepic.bioinformatics;

import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class Genome {
    private final String text;

    public String getText() {
        return text;
    }

    public Genome(String text) {
        this.text = text;
    }

    private final ConcurrentHashMap <String,Long> patternCountCache =
            new ConcurrentHashMap<>();

    public final long patternCount(final String pattern) {
        Long result = patternCountCache.get(pattern);
        if (result == null) {
            result = Stream.iterate(0, i->i+1)
                .limit(text.length() - pattern.length() + 1)
                .filter(
                    i -> text.substring(i, i + pattern.length()).equals(pattern)
                )
                .count();
            patternCountCache.put(pattern,result);
        }
        return result;
    }

    public final HashMap<String,Long> getKmers(final int k){
        return Stream
                .iterate(0,i->i+1).limit(text.length() - k + 1)
                .map( i -> text.substring(i,i+k) )
                .collect(Collectors.toMap(
                        (String str) -> str,
                        (String str) -> patternCount(str),
                        (Long old_value, Long new_value) -> old_value,
                        HashMap::new
                ));
    }

    public final List<String> frequentWords(final int k) {
        Set<String> set = Stream
            .iterate(0,i->i+1).limit(text.length() - k + 1)
            .map( i -> text.substring(i,i+k) )
            .collect(
                    HashSet::new,
                    (HashSet<String> hashSet, String str) -> {
                        if (!hashSet.contains(str)) {
                            long max = !hashSet.isEmpty() ? patternCount(hashSet.iterator().next()) : 0;
                            long kmer = patternCount(str);
                            if ( kmer > max) {
                                hashSet.clear();
                                hashSet.add(str);
                            }
                            else if ( kmer == max ) {
                                hashSet.add(str);
                            }
                        }
                    },
                    HashSet::addAll
            );
        return new ArrayList<>(set);
    }

    public final List<String> frequentWordsImp(final int k) {
        HashMap<String, Long> kmers = getKmers(k);
        long max = 0;
        List<String> listOfMax = new ArrayList<>();
        for (String key : kmers.keySet() ) {
            if (kmers.get(key) > max ) {
                listOfMax.clear();
                listOfMax.add(key);
                max = kmers.get(key);
            }
            else if ( kmers.get(key) == max ) {
                listOfMax.add(key);
            }
        }
        return listOfMax;
    }

    public static final String reverseComplement(final String pattern) {
        Function<Character, Character> reverseLetter =
                letter -> {
                    char result = letter;
                    switch (letter) {
                        case 'A':
                            result = 'T';
                            break;
                        case 'a':
                            result = 't';
                            break;
                        case 'T':
                            result = 'A';
                            break;
                        case 't':
                            result = 'a';
                            break;
                        case 'G':
                            result = 'C';
                            break;
                        case 'g':
                            result = 'c';
                            break;
                        case 'C':
                            result = 'G';
                            break;
                        case 'c':
                            result = 'g';
                            break;
                    }
                    return result;
                };

        return Stream.iterate(pattern.length()-1,i->i-1).limit(pattern.length())
                .map(i -> reverseLetter.apply(pattern.charAt(i)))
                .collect(
                        StringBuilder::new,
                        StringBuilder::append,
                        StringBuilder::append
                ).toString();
    }

    public final List<Integer> findAllPossionOfThe(final String pattern) {
        int k = pattern.length();
        String reverse = reverseComplement(pattern);
        return Stream
            .iterate(0,i->i+1).limit(text.length() - k + 1)
            .filter(i-> text.substring(i,i+k).equals(pattern))
            .collect(Collectors.toList());
    }


}
