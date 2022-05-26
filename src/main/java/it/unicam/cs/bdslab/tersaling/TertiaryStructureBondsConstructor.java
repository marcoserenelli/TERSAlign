package it.unicam.cs.bdslab.tersaling;

import org.biojava.nbio.structure.contact.Pair;

import java.util.ArrayList;
import java.util.Objects;

public class TertiaryStructureBondsConstructor extends TertiaryStructureBondsBaseListener {

    private final ArrayList<Pair<Integer>> bondsList;

    public TertiaryStructureBondsConstructor(){
        this.bondsList = new ArrayList<>();
    }


    public ArrayList<Pair<Integer>> getBondsList() {
        return bondsList;
    }

    @Override
    public void enterBondsContinue(TertiaryStructureBondsParser.BondsContinueContext ctx){
        //add bonds to the list
        int left = Integer.parseInt(ctx.bond().INDEX(0).getText());
        int right = Integer.parseInt(ctx.bond().INDEX(1).getText());
        this.bondsList.add(new Pair<>(left, right));
    }

    @Override
    public void exitBond(TertiaryStructureBondsParser.BondContext ctx) {
        // everything has been added to the structure, now sort every pair in ascending order
        this.bondsList.sort((o1, o2) -> {
            //if the first number is equals
            if (Objects.equals(o1.getFirst(), o2.getFirst())) {
                //we check the second
                if (o1.getSecond() < o2.getSecond()) {
                    //o1 is more in the left
                    return -1;
                } else if (o1.getSecond() > o2.getSecond()) {
                    //o1 is more in the right
                    return 1;
                }
                //error first and second shouldn't both be equals
                throw new IllegalArgumentException("error first and second shouldn't both be equals");
            } else if (o1.getFirst() < o2.getFirst()) {
                //o1 is more in the left
                return -1;
            }
            //o1 is more in the right
            return 1;
        });
    }

}
