import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.GroupType;
import org.biojava.nbio.structure.Structure;

public class TertiaryStructure {

    private final Structure structure;
    private double threshold;

    public TertiaryStructure(Structure structure, double threshold) {
        this.structure = structure;
        this.threshold = threshold;
    }

    public double[][] getDistanceMatrix(){
        return calculateDistanceMatrix(this.structure);
    }

    /**
     * Create distance matrix based on a structure, for nucleotides takes distance between P as comparison method,
     * for amino acid takes distance between CA as comparison method
     * @param struc structure
     * @return distance matrix
     */
    public double[][] calculateDistanceMatrix(Structure struc){
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
    public double[][] calculateDistanceMatrixCenterOfMass(Structure struc){
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
    public static void printDistanceMatrix(double[][] distanceMatrix){
        for (int i=0; i<distanceMatrix.length; i++) {
            for (int j=0; j<distanceMatrix.length; j++) {
                System.out.print("pos " + i + " " + j + ": " + distanceMatrix[i][j] + " ");
            }
            System.out.println();
        }
    }

    /**
     * Return a boolean matrix, values are true if their distance (taken from distanceMatrix)
     * is less than threshold value.
     * @return boolean contact matrix
     */
    public boolean[][] getContactMatrix(){
        double[][] distanceMatrix = getDistanceMatrix();
        boolean[][] distanceBool = new boolean[distanceMatrix.length][distanceMatrix.length];
        for (int i=0; i<distanceMatrix.length; i++) {
            for (int j = 0; j < distanceMatrix.length; j++) {
                distanceBool[i][j] = distanceMatrix[i][j] <= this.threshold;
                System.out.print(distanceBool[i][j] + " distance is " + distanceMatrix[i][j] + "  ");
            }
            System.out.println();
        }
        return distanceBool;
    }

    /**
     * Returns the predominantGroupType of a structure chain
     * @return group type
     */
    public GroupType getType(){
        for(Chain c : this.structure.getChains()){
            if(c.getPredominantGroupType() == GroupType.AMINOACID){
                return GroupType.AMINOACID;
            }
            else if(c.getPredominantGroupType() == GroupType.NUCLEOTIDE){
                return GroupType.NUCLEOTIDE;
            }
        }
        return null;
    }

    public double getThreshold() {
        return threshold;
    }

    public void setThreshold(double threshold) {
        this.threshold = threshold;
    }

    public Structure getStructure() {
        return structure;
    }
}
