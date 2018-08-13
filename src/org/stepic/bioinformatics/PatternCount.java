package org.stepic.bioinformatics;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;


public class PatternCount {
    public static final long patternCount(String text, String pattern) {
        return Stream.iterate(0,i->i+1)
                .limit(text.length() - pattern.length() + 1)
                .filter(
                        i -> text.substring(i, i + pattern.length()).equals(pattern)
                )
                .count();
    }


    public static final List<String> frequentWordsAll(String text, int k) {
        HashMap<String, Long> kmers = new HashMap<>();
        List<String> listOfMax = new ArrayList<>();
        long max = 0;
        for (int i = 0; i < text.length() - k; i++) {
            String str = text.substring(i, i + k);
            if (!kmers.containsKey(str)) {
                long kmer = patternCount(text, str);
                kmers.put(str, kmer);
                if (kmer > max) {
                    listOfMax = new ArrayList<>();
                    listOfMax.add(str);
                    max = kmer;
                }
                else if ( kmer == max ) {
                    listOfMax.add(str);
                }
            }
        }
        return listOfMax;
    }
    public static final HashMap<String,Long> getKmers(String text, int k){
    return Stream
            .iterate(0,i->i+1).limit(text.length() - k + 1)
            .map( i -> text.substring(i,i+k) )
            .collect(Collectors.toMap(
                    (String str) -> str,
                    (String str) -> patternCount(text,str),
                    (Long old_value, Long new_value) -> old_value,
                    HashMap::new
            ));
    }

    public static final List<String> frequentWords(String text, int k) {
        HashMap<String, Long> kmers = getKmers(text, k);
        long max = 0;
        List<String> listOfMax = new ArrayList<>();
        for (String key : kmers.keySet() ) {
            if (kmers.get(key) > max ) {
                listOfMax = new ArrayList<>();
                listOfMax.add(key);
                max = kmers.get(key);
            }
            else if ( kmers.get(key) == max ) {
                listOfMax.add(key);
            }
        }
        return listOfMax;
    }

    public static void main(String[] args) {
        System.out.println(getKmers("GACGATATACGACGATA",3));
    }
}
