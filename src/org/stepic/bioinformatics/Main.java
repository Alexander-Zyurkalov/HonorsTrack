package org.stepic.bioinformatics;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.time.Duration;
import java.time.LocalDateTime;

import static java.util.stream.Collectors.joining;

public class Main {

    public static void main(String[] args) {
        var sb = new StringBuilder();
        try (var fileReader = new FileReader("Salmonella_enterica.txt");
             final var bufferedReader = new BufferedReader(fileReader);
        ) {
            while (true) {
                var line = bufferedReader.readLine();
                if ( line == null) break;
                if (line.matches("^\\>"))
                    continue;
                sb.append(line);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        var genome = new FuzzyGenome(sb.toString());
        var start = LocalDateTime.now();
        System.out.println(start);

        var skew = genome.minimumSkewImp();
        skew.forEach(
                (v) -> System.out.print("skew = " + v + " ")
        );
        System.out.println();

        var window = new FuzzyGenome(genome,skew.get(0)-500,1000);
        System.out.println(window.frequentWordsWithMismatch(9,1));
        window.findAllPositionsOfTheApproximatePattern(new Sequence("CCGGATCCG"),1);


        var stop = LocalDateTime.now();
        System.out.println(stop);
        var duration = Duration.between(start,stop);
        System.out.println(duration.getSeconds());

//        System.out.println(clump.size());
//        System.out.println(clump.stream().map(s->s.toString()).collect(Collectors.joining(" ")));
    }

}
