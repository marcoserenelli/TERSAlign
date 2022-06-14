package it.unicam.cs.bdslab.tersaling;
import org.antlr.v4.runtime.ParserRuleContext;
import org.antlr.v4.runtime.tree.ErrorNode;
import org.antlr.v4.runtime.tree.TerminalNode;
import org.biojava.nbio.structure.contact.Pair;

import java.util.ArrayList;
import java.util.Objects;

public class TertiaryStructureBondsOptionalSequenceConstructor extends TertiaryStructureBondsOptionalSequenceBaseListener {
    private String sequence;
    private final StringBuffer sequenceBuffer;
    private final ArrayList<Pair<Integer>> bondsList;

    public TertiaryStructureBondsOptionalSequenceConstructor() {
        this.sequence = "";
        this.sequenceBuffer = new StringBuffer();
        this.bondsList = new ArrayList<>();
    }

    public String getSequence() {
        return sequence;
    }

    public ArrayList<Pair<Integer>> getBondsList() {
        return bondsList;
    }

    @Override
    public void enterSequenceContinue(TertiaryStructureBondsOptionalSequenceParser.SequenceContinueContext ctx) {
        System.out.println("ENTER SEQ CONTINUE");
        this.sequenceBuffer.append(ctx.LETTERS().getText());
    }

    @Override
    public void exitSequenceEnd(TertiaryStructureBondsOptionalSequenceParser.SequenceEndContext ctx) {
        System.out.println("EXIT SEQ END");
        this.sequenceBuffer.append(ctx.LETTERS().getText());
        this.sequence = this.sequenceBuffer.toString();
    }

    @Override
    public void enterBondsContinue(TertiaryStructureBondsOptionalSequenceParser.BondsContinueContext ctx){
        //add bonds to the list
        System.out.println("ENTERCONTINUE");
        int left = Integer.parseInt(ctx.bond().INDEX(0).getText());
        int right = Integer.parseInt(ctx.bond().INDEX(1).getText());
        this.bondsList.add(new Pair<>(left, right));
        System.out.println(bondsList);
    }

    @Override
    public void exitBondsEnd(TertiaryStructureBondsOptionalSequenceParser.BondsEndContext ctx) {
        System.out.println("EXITEND");
    }

    @Override
    public void enterBondsEnd(TertiaryStructureBondsOptionalSequenceParser.BondsEndContext ctx) {
        System.out.println("ENTEREND");
    }

    @Override
    public void enterBond(TertiaryStructureBondsOptionalSequenceParser.BondContext ctx) {
        System.out.println("ENTERBOND");
    }

    @Override
    public void exitBond(TertiaryStructureBondsOptionalSequenceParser.BondContext ctx) {
        System.out.println("EXITBOND");
    }

    @Override
    public void exitBondsContinue(TertiaryStructureBondsOptionalSequenceParser.BondsContinueContext ctx) {
        System.out.println("EXITCONTINUE");
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
