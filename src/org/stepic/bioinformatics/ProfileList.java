package org.stepic.bioinformatics;

import java.util.ArrayList;
import java.util.Map;
import java.util.stream.Stream;

public class ProfileList extends ArrayList<Map<Character,Double>> {
    public ProfileList(){
        super();
    }
    public Stream<Double> getAllValuesStream() {
        return this.stream()
                .flatMap(hash -> hash.values().stream());
    }
}
