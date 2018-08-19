package org.stepic.bioinformatics;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Function;
import java.util.stream.IntStream;
import java.util.stream.Stream;

public class Sequence {
    private final String text;
    private final int beginning;
    private final int end;
    private Sequence parent;
    int hash = -1;
    private int P = 4;

    private static char[] TO_LETTER = {'A','C','G','T'};

    public static final Function<Character,Byte> TO_NUM = ch -> {
        byte b;
        switch (ch) {
            case 'A':
            case 'a':
                b = 0;
                break;
            case 'C':
            case 'c':
                b = 1;
                break;
            case 'G':
            case 'g':
                b = 2;
                break;
            case 'T':
            case 't':
                b = 3;
                break;
            default:
                b = 0;
        }
        return b;
    };
    public final Function<Integer,Integer> pow = n -> IntStream.range(0,n).reduce(1,(m, i) -> m*P);

    public long patternToNumberLong(){
        long s = 0;
        for (int i = 0; i < length(); i++) {
            s += TO_NUM.apply(charAt(i)) * pow.apply(length() - i - 1);
        }
        return s;
    }

    public int patternToNumber(){
        int result =
                IntStream.range(0,length())
                .map(i-> TO_NUM.apply(charAt(i)) * pow.apply(length() - i - 1)  )
                .reduce(0,(s,i) -> s + i );
        return result;
    }
    public long patternToNumberRecursively(){
        if (length() == 0)
            return 0;
        else {
            var prefix = new Sequence(this,0,length()-1);
            var result = 4 * prefix.patternToNumberRecursively() + TO_NUM.apply(charAt(length()-1));
            System.out.println(result);
            return result;
        }

    }

    public static Sequence numberToPattern(long num, int k){
        return new Sequence(
            Stream.iterate(num,i-> i / 4 ).limit(k)
                .map(i-> TO_LETTER[(int)(i % 4)])
                .collect(
                        StringBuilder::new,
                        StringBuilder::append,
                        StringBuilder::append
                ).reverse().toString()
        );
    }

    @Override
    public int hashCode() {
//        return hash == -1 ? hash = patternToNumber() : hash;
        return getText().hashCode();
    }

    @Override
    public boolean equals(Object anObject) {
        if (this == anObject) {
            return true;
        }
        if (anObject instanceof Sequence) {
            var seq = (Sequence) anObject;
            if (seq.length() == this.length())
                return this.equalsAt(seq,0);
        }
        return false;
    }

    public boolean equalsAt(String pattern, int beginning){
        return equalsAt(new Sequence(pattern),beginning);
    }

    public boolean equalsAt(Sequence seq, int beginning){
        return this.length() - beginning >=  seq.length()  &&
                    IntStream.range(0,seq.length())
                        .allMatch(i->this.charAt(i+beginning) == seq.charAt(i));
    }

    public char charAt(int i){
        return text.charAt(beginning+i);
    }

    public Sequence(final Sequence sequence, final int beginning, final int length) {
        this.parent = sequence;
        this.text = sequence.text;
        var new_beginning = sequence.beginning + beginning;
        this.beginning = new_beginning <= sequence.end ? new_beginning : sequence.end;
        var new_end = this.beginning + length - 1;
        this.end = new_end <= sequence.end ? new_end : sequence.end;
    }

    public Sequence(final String text) {
        beginning = 0;
        end = text.length() - 1;
        this.text = text;
    }

    public String getText() {
//        if (beginning == 0 && end + 1 == text.length())
//            return text;  // !!!!!!! the block is not needed because substring does the same
//        else
            return text.substring(beginning,end+1);
    }

    public int length(){
        return end - beginning + 1;
    }

    @Override
    public String toString() {
        return getText();
    }

    public final Sequence reverseComplement() {
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

        return
            new Sequence(
                Stream.iterate(this.length()-1,i->i-1).limit(this.length())
                        .map(i -> reverseLetter.apply(this.charAt(i)))
                        .collect(
                                StringBuilder::new,
                                StringBuilder::append,
                                StringBuilder::append
                        ).toString()
            );
    }
    public boolean theyAreTheSame(Sequence seq){
        return this.getText() == seq.text;
    }

    public int hammingDistance(Sequence another){
        if (length()!= another.length())
            return length() ;
        return (int)Stream.iterate(0,i->i+1).limit(length())
                .filter(i->charAt(i) != another.charAt(i))
                .count();
    }

    public static int hammingDistance(String s1, String s2){
        var seq = new Sequence(s1);
        return seq.hammingDistance(new Sequence(s2));
    }

    public List<Sequence> getNeighbors(int distance){
        var neighborhood = new ArrayList<Sequence>();
        if (distance == 0) {
            neighborhood.add(this);
            return neighborhood;
        }
        if ((length()) == 1) {
            neighborhood.add(new Sequence("A"));
            neighborhood.add(new Sequence("C"));
            neighborhood.add(new Sequence("G"));
            neighborhood.add(new Sequence("T"));
            return neighborhood;
        }


        var suffixPattern = new Sequence(
                this,
                1,
                length() - 1);
        suffixPattern.getNeighbors(distance).forEach(
            text -> {

                if (suffixPattern.hammingDistance(text) < distance){
                    neighborhood.add(new Sequence('A' + text.getText()));
                    neighborhood.add(new Sequence('C' + text.getText()));
                    neighborhood.add(new Sequence('G' + text.getText()));
                    neighborhood.add(new Sequence('T' + text.getText()));
                }
                else
                    neighborhood.add(new Sequence(charAt(0) + text.getText()));


            }
        );
        return neighborhood;
    }

    public static void main(String[] args) {
        var seq = new Sequence("CCAGTCAATG");
        System.out.println(seq.getNeighbors(1).size());
        seq.getNeighbors(1).forEach(
                sequence -> System.out.println(sequence)
        );
    }



}
