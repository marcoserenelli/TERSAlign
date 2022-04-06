import org.biojava.nbio.structure.*;

public class RepAtmsDistanceMatrixGenerator {

    /**
     * Create distance matrix based on a structure, for nucleotides takes distance between P as comparison method,
     * for amino acid takes distance between CA as comparison method
     * @param struc structure
     * @return distance matrix
     */
    public static double[][] calculateDistanceMatrix(Structure struc){
        Atom[] representativeAtomsArray = StructureTools.getRepresentativeAtomArray(struc);
        double[][] distanceMatrix = new double[representativeAtomsArray.length][representativeAtomsArray.length];
        for(int i=0; i<representativeAtomsArray.length; i++)
            for(int j=0; j<representativeAtomsArray.length; j++)
                distanceMatrix[i][j] = Calc.getDistance(representativeAtomsArray[i],representativeAtomsArray[j]);
        return distanceMatrix;
    }

    /**
     * Create distance matrix based on a structure, considering aminos / nucleotide's center of mass
     * @param struc structure
     * @return distance matrix
     */
    public static double[][] calculateDistanceMatrixCenterOfMass(Structure struc){
        int groupsNumber = StructureTools.getNrGroups(struc);
        double[][] distanceMatrix = new double[groupsNumber][groupsNumber];
        int moleculeCount = 0;
        for(Chain currentChain: struc.getChains()) {
            for (Group currentMolecule : currentChain.getAtomGroups()) {
                int comparedMoleculeCount = 0;
                for (Chain comparisonChain : struc.getChains()) {
                    for (Group comparisonMolecule : comparisonChain.getAtomGroups()) {
                        distanceMatrix[moleculeCount][comparedMoleculeCount] = Calc.getDistance(Calc.centerOfMass(currentMolecule.getAtoms().toArray(new Atom[0])), Calc.centerOfMass(comparisonMolecule.getAtoms().toArray(new Atom[0])));
                        comparedMoleculeCount++;
                    }
                }
                moleculeCount++;
            }
        }
        return distanceMatrix;
    }

    /**
     * Print the matrix
     * @param distanceMatrix matrix to print
     */
    public static void printAmminoacidDistanceMatrix(double[][] distanceMatrix){
        for (int i=0; i<distanceMatrix.length; i++) {
            for (int j=0; j<distanceMatrix.length; j++) {
                System.out.print("pos " + i + " " + j + ": " + distanceMatrix[i][j] + " ");
            }
            System.out.println();
        }
    }

}
