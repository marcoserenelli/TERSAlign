import org.biojava.nbio.phylo.DistanceMatrixCalculator;
import org.biojava.nbio.structure.*;
import org.forester.evoinference.matrix.distance.DistanceMatrix;

import java.sql.SQLOutput;
import java.util.Arrays;

public class RMSDistanceMatrixGenerator {

    public static DistanceMatrix calculateDistanceMatrix(Structure struc) {
        double[][] rmsdMat = RMSDistanceMatrixGenerator.calculateRMSDMat(struc);
        return DistanceMatrixCalculator.structuralDistance(rmsdMat,0,0,0);
    }

    public static double[][] calculateRMSDMat(Structure struc){
        double[][] rmsdMat = new double[StructureTools.getNrAtoms(struc)][StructureTools.getNrAtoms(struc)];
        for(Chain currentChain : struc.getChains()){
            int groupCounter = 0;
            for(Group currentGroup : currentChain.getAtomGroups()){
                for(Chain comparisonChain : struc.getChains()) {
                    int comparisonGroupCounter = 0;
                    for (Group comparisonGroup : comparisonChain.getAtomGroups()) {
                        rmsdMat[groupCounter][comparisonGroupCounter] = Calc.rmsd(currentGroup.getAtoms().toArray(new Atom[0]), comparisonGroup.getAtoms().toArray(new Atom[0]));
                        comparisonGroupCounter++;
                    }
                }
                groupCounter++;
            }
        }
        return rmsdMat;
    }
}