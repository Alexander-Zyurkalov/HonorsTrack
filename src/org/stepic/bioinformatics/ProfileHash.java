package org.stepic.bioinformatics;

import java.util.HashMap;
import java.util.List;

public class ProfileHash extends HashMap<Character,List<Double>> {
    public ProfileList toProfileList(){
        return Probability.hashToList(this);
    }

}