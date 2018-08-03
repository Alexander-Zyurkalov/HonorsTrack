package org.stepic.bioinformatics;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.stream.Collector;

import static java.util.stream.Collectors.joining;
import static org.stepic.bioinformatics.PatternCount.patternCount;

public class Main {

    public static void main(String[] args) {

//        Genome genome = new Genome("GGACTTACTGACGTACG");
//        System.out.println(genome.frequentWords(3));
//        System.out.println(Genome.reverseComplement("tgATGATCAAG"));
        String data = new String();
        try (FileReader fileReader = new FileReader(
                "Vibrio_cholerae.txt");
             BufferedReader bufferedReader = new BufferedReader(fileReader)
        ) {
            String line;
            do {
                line = bufferedReader.readLine();
                data = data + line;
            } while (line != null);

        } catch (IOException e) {
            e.printStackTrace();
        }
        Genome genome = new Genome(data);

        System.out.println(
                genome.findAllPossionOfThe("CTTGATCAT")
                .stream()
                .map(i->String.valueOf(i))
                .collect(joining(" "))
        );

    }

}
