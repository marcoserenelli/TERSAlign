package it.unicam.cs.bdslab.tersaling;

import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.contact.Pair;

import java.util.ArrayList;
import java.util.Objects;

public class TertiaryStructuresComparator {

    public static void main(String[] args) {
        try{
            Structure structure = StructureIO.getStructure("3mge");
            TertiaryStructure tertiaryStructure = new TertiaryStructure(structure);
            System.out.println(tertiaryStructure.getSequence().length());
            System.out.println(tertiaryStructure.getBondList());
            StructuralTree treeBuilder = new StructuralTree(tertiaryStructure);
            System.out.println(TreeOutputter.treeToString(treeBuilder.getStructuralRNATree()));
            Structure structure2 = StructureIO.getStructure("5r7y");
            TertiaryStructure tertiaryStructure2 = new TertiaryStructure(structure2);
            System.out.println(tertiaryStructure2.getSequence().length());
            System.out.println(tertiaryStructure2.getBondList());
            StructuralTree treeBuilder2 = new StructuralTree(tertiaryStructure2);
            System.out.println(TreeOutputter.treeToString(treeBuilder2.getStructuralRNATree()));
        } catch (Exception e){
            e.printStackTrace();
        }
    }
}

