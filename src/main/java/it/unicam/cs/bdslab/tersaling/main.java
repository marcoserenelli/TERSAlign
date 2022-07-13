package it.unicam.cs.bdslab.tersaling;

import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Locale;

/**
 * @author Filippo Lampa
 */
public class main {
    public static void main(String[] args) throws StructureException, IOException {
        Structure struc = StructureIO.getStructure("5a9q");
        TertiaryStructure ts = new TertiaryStructure(struc);
        ArrayList<String> idsArray = new ArrayList<String>();
        idsArray.add("0");
        idsArray.add("1");
        ts.setSpecifiedChains(idsArray);
        ts.setDistanceMatrixCalculationMethod("centerofmass");
        ts.getDistanceMatrix();
        ts.printDistanceMatrix();
        /*struc.getEntityInfos().forEach(i -> {
            Chain selectedChain;
            if(i.getDescription().equalsIgnoreCase("NUCLEOPORIN NUP43")) {
                selectedChain = struc.getChain(i.getChainIds().get(0));
                selectedChain.getAtomGroups().forEach(group -> System.out.println(group.getPDBName()));
            }
        });*/
    }
}
