package org.stepic.bioinformatics;

import java.time.Duration;
import java.time.LocalDateTime;
import java.util.function.Consumer;

public class TimeMeasure {
    public static void measure(Runnable runnable){
        var start = LocalDateTime.now();
        System.out.println(start);

        runnable.run();

        var stop = LocalDateTime.now();
        System.out.println(stop);
        var duration = Duration.between(start,stop);
        System.out.println(duration.getNano());
    }
}
