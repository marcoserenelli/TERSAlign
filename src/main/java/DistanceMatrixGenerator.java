import org.biojava.nbio.structure.*;

public class DistanceMatrixGenerator {

    public static Double[][] calcuateDistanceMatrix(Structure struc){

        Double[][] distanceMatrix = generateMatrix(struc);
        int moleculeCount = 0;
        for(Chain currentChain: struc.getChains()) {
            Atom mainMoleculeAtom;
            for (Group currentMolecule : currentChain.getAtomGroups()) {
                if(!(currentMolecule.getType() == GroupType.HETATM)) {
                    if (currentChain.isNucleicAcid() || currentChain.isProtein()) {
                        mainMoleculeAtom = currentChain.isProtein() ? currentMolecule.getAtom("CA") : currentMolecule.getAtom("P");
                        int comparedMoleculeCount = 0;
                        //Compare the main molecule's atom's (CA or P) position with every main molecule's atom's position of every chain inside the structure
                        for (Chain comparisonChain : struc.getChains()) {
                            for (Group comparisonMolecule : comparisonChain.getAtomGroups()) {
                                if(!(comparisonMolecule.getType() == GroupType.HETATM)) {
                                    Atom comparedMainMoleculeAtom;
                                    if (currentChain.isNucleicAcid() || currentChain.isProtein()) {
                                        comparedMainMoleculeAtom = comparisonChain.isProtein() ? comparisonMolecule.getAtom("CA") : comparisonMolecule.getAtom("P");
                                        distanceMatrix[moleculeCount][comparedMoleculeCount] = Calc.getDistance(mainMoleculeAtom, comparedMainMoleculeAtom);
                                        comparedMoleculeCount++;
                                    }
                                }
                            }
                        }
                    }
                    moleculeCount++;
                }
            }
        }
        return distanceMatrix;
    }

    private static Double[][] generateMatrix(Structure struc) {
        int totalSequenceLength = 0;
        for(Chain currentChain : struc.getChains()) {
            for(Group currentGroup : currentChain.getAtomGroups()){
                if(!(currentGroup.getType() == GroupType.HETATM))
                    totalSequenceLength++;
            }
        }
        return new Double[totalSequenceLength][totalSequenceLength];
    }

    public static void printAmminoacidDistanceMatrix(Double[][] distanceMatrix){
        for (int i=0; i<distanceMatrix.length; i++) {
            for (int j=0; j<distanceMatrix.length; j++) {
                System.out.print("pos " + i + " " + j + ": " + distanceMatrix[i][j] + " ");
            }
            System.out.println();
        }
    }

}
