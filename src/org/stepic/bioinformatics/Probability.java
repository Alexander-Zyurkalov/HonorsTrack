package org.stepic.bioinformatics;

import java.math.BigDecimal;
import java.util.Arrays;
import java.util.List;
import java.util.stream.DoubleStream;
import java.util.stream.Stream;

public class Probability {


    private static BigDecimal sumNonMutuallyExclusive(BigDecimal p1, BigDecimal p2){
        return p1.add(p2).add( p1.multiply(p2).negate());
//        return p1.multiply(p2).negate();
    }
    public static BigDecimal  sumNonMutuallyExclusive(List<Double> probabilities) {
        return sumNonMutuallyExclusive(probabilities.stream().mapToDouble(n->n));
    }
    public static BigDecimal sumNonMutuallyExclusive(DoubleStream stream) {
        return stream.mapToObj(BigDecimal::valueOf)
                .reduce(Probability::sumNonMutuallyExclusive)
                .orElse(BigDecimal.ZERO);
    }

    public static void main(String[] args) {
        final var one_mer = BigDecimal.valueOf(0.25);
        final var nine_mer = Stream
                .iterate(one_mer,i->one_mer).limit(9)
                .reduce(BigDecimal::multiply).orElse(one_mer);
//        System.out.println(nine_mer);
        System.out.println(
            sumNonMutuallyExclusive(
                DoubleStream.iterate(nine_mer.doubleValue(),i->nine_mer.doubleValue())
                        .limit((1000-9+1) * 500)
            ).multiply(BigDecimal.valueOf((1000-9+1) * 500))
        );
    }
}
