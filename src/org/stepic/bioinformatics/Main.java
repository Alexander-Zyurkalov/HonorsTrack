package org.stepic.bioinformatics;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.time.Duration;
import java.time.LocalDateTime;
import java.time.Period;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collector;
import java.util.stream.Collectors;

import static java.util.stream.Collectors.joining;

public class Main {

    public static void main(String[] args) {

        String text = "";
        try (var fileReader = new FileReader("E_coli.txt ");
             var bufferedReader = new BufferedReader(fileReader)
        ) {
            text = bufferedReader.readLine();
        } catch (IOException e) {
            e.printStackTrace();
        }
        var genome = new Genome(text);
        var start = LocalDateTime.now();
        System.out.println(start);

        var clump = genome.clumpFinding(9,500,3);

        var stop = LocalDateTime.now();
        System.out.println(stop);
        var duration = Duration.between(start,stop);
        System.out.println(duration.getSeconds());
        System.out.println(clump.size());
//        System.out.println(clump.stream().map(s->s.toString()).collect(Collectors.joining(" ")));
    }

}
