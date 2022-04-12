import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.GroupType;
import org.biojava.nbio.structure.Structure;
import org.jgrapht.alg.util.Pair;

import java.util.ArrayList;

public class TertiaryStructure {

    private final Structure structure;
    private double threshold;

    public TertiaryStructure(Structure structure) {
        this.structure = structure;
        this.threshold = 4.0;
    }

    /**
     * Return the distance matrix using centerOfMass
     * @return distance matrix using center of mass method
     */
    public double[][] getDistanceMatrixCenterOfMassIgnoreHeta(){
        return calculateDistanceMatrixCenterOfMassIgnoreHeta(this.structure);
    }

    public double[][] getDistanceMatrixCenterOfMass(){
        return calculateDistanceMatrixCenterOfMass(this.structure);
    }

    /**
     * Return the distance matrix using default calculation method
     * @return distance matrix using default calculation method
     */
    public double[][] getDistanceMatrixDefault(){
        return calculateDistanceMatrix(this.structure);
    }

    /**
     * Generates a list containing the indexes of bonded nucleotides/aminos, that is all nucleotides/aminos closer than the specified threshold.
     * @return bond list
     */
    public ArrayList<Pair<Integer, Integer>> getBondList(){
        boolean[][]contactMap = this.getContactMatrixCenterOfMassIgnoreHeta();
        ArrayList<Pair<Integer, Integer>>bondList = new ArrayList<>();
        int colCount = 0;
        for(int i=0; i<contactMap.length; i++) {
            for (int j = colCount; j < contactMap.length; j++)
                if (contactMap[i][j]) {
                    bondList.add(new Pair<>(i, j));
                }
            colCount++;
        }
        return bondList;
    }

    /**
     * Return a boolean matrix, values are true if their distance (taken from distanceMatrix)
     * is less than threshold value.
     * @return boolean contact matrix
     */
    public boolean[][] getContactMatrixCenterOfMassIgnoreHeta(){
        double[][] distanceMatrix = getDistanceMatrixCenterOfMassIgnoreHeta();
        boolean[][] contactMatrix = new boolean[distanceMatrix.length][distanceMatrix.length];
        for (int i=0; i<distanceMatrix.length; i++) {
            for (int j = 0; j < distanceMatrix.length; j++) {
                contactMatrix[i][j] = distanceMatrix[i][j] <= this.threshold;
            }
        }
        return contactMatrix;
    }

    public boolean[][] getContactMatrixCenterOfMass(){
        double[][] distanceMatrix = getDistanceMatrixCenterOfMass();
        boolean[][] contactMatrix = new boolean[distanceMatrix.length][distanceMatrix.length];
        for (int i=0; i<distanceMatrix.length; i++) {
            for (int j = 0; j < distanceMatrix.length; j++) {
                contactMatrix[i][j] = distanceMatrix[i][j] <= this.threshold;
            }
        }
        return contactMatrix;
    }

    /**
     * Return a boolean matrix, values are true if their distance (taken from default calculation)
     * is less than threshold value.
     * @return boolean contact matrix
     */
    public boolean[][] getContactMatrixDefault(){
        double[][] distanceMatrix = getDistanceMatrixDefault();
        boolean[][] contactMatrix = new boolean[distanceMatrix.length][distanceMatrix.length];
        for (int i=0; i<distanceMatrix.length; i++) {
            for (int j = 0; j < distanceMatrix.length; j++) {
                contactMatrix[i][j] = distanceMatrix[i][j] <= this.threshold;
            }
        }
        return contactMatrix;
    }

    /**
     * Create distance matrix based on a structure, for nucleotides takes distance between P as comparison method,
     * for amino acid takes distance between CA as comparison method
     * @param struc structure
     * @return distance matrix
     */
    public double[][] calculateDistanceMatrix(Structure struc){
        Atom[] representativeAtomsArray = StructureTools.getRepresentativeAtomArray(struc);    //questa operazione ignora tutti gli hetatms
        double[][] distanceMatrix = new double[representativeAtomsArray.length][representativeAtomsArray.length];
        for(int i=0; i<representativeAtomsArray.length; i++)
            for(int j=0; j<representativeAtomsArray.length; j++)
                distanceMatrix[i][j] = Calc.getDistance(representativeAtomsArray[i], representativeAtomsArray[j]);
        printDistanceMatrix(distanceMatrix);
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
        for(Chain currentChain: struc.getChains())
            for (Group currentMolecule : currentChain.getAtomGroups()) {
                int comparedMoleculeCount = 0;
                for (Chain comparisonChain : struc.getChains())
                    for (Group comparisonMolecule : comparisonChain.getAtomGroups()) {
                        if(currentMolecule.getType() != GroupType.HETATM && comparisonMolecule.getType() != GroupType.HETATM)
                            distanceMatrix[moleculeCount][comparedMoleculeCount] = Calc.getDistance(Calc.centerOfMass(currentMolecule.getAtoms().toArray(new Atom[0])), Calc.centerOfMass(comparisonMolecule.getAtoms().toArray(new Atom[0])));
                        else distanceMatrix[moleculeCount][comparedMoleculeCount] = Double.MAX_VALUE;
                        comparedMoleculeCount++;
                    }
                moleculeCount++;
            }
        printDistanceMatrix(distanceMatrix);
        return distanceMatrix;
    }


    /**
     * Create distance matrix based on a structure, considering aminos / nucleotide's center of mass
     * @param struc structure
     * @return distance matrix
     */
    public double[][] calculateDistanceMatrixCenterOfMassIgnoreHeta(Structure struc){
        int groupsNumber = getNonHetatmGroupsCounter(struc);
        double[][] distanceMatrix = new double[groupsNumber][groupsNumber];
        int moleculeCount = 0;
        for(Chain currentChain: struc.getChains())
            for (Group currentMolecule : currentChain.getAtomGroups()) {
                if (currentMolecule.getType() != GroupType.HETATM) {
                    int comparedMoleculeCount = 0;
                    for (Chain comparisonChain : struc.getChains()) {
                        for (Group comparisonMolecule : comparisonChain.getAtomGroups()) {
                            if (comparisonMolecule.getType() != GroupType.HETATM) {
                                distanceMatrix[moleculeCount][comparedMoleculeCount] = Calc.getDistance(Calc.centerOfMass(currentMolecule.getAtoms().toArray(new Atom[0])), Calc.centerOfMass(comparisonMolecule.getAtoms().toArray(new Atom[0])));
                                comparedMoleculeCount++;
                            }
                        }
                    }
                    moleculeCount++;
                }
            }
        printDistanceMatrix(distanceMatrix);
        return distanceMatrix;
    }

    private int getNonHetatmGroupsCounter(Structure struc){
        int nonHetatmGroupsCounter = 0;
        for(Chain currentChain : struc.getChains())
            for(Group currentGroup : currentChain.getAtomGroups())
                if(currentGroup.getType() != GroupType.HETATM)
                    nonHetatmGroupsCounter++;
        return nonHetatmGroupsCounter;
    }
    /**
     * Print the distance matrix
     * @param distanceMatrix matrix to print
     */
    public void printDistanceMatrix(double[][] distanceMatrix){
        for (int i=0; i<distanceMatrix.length; i++) {
            for (int j=0; j<distanceMatrix.length; j++) {
                System.out.print("pos " + i + " " + j + ": " + distanceMatrix[i][j] + " ");
            }
            System.out.println();
        }
    }

    /**
     * Print the contact matrix
     * @param contactMatrix contact matrix to print
     */
    public void printContactMatrix(boolean[][] contactMatrix){
        for (int i=0; i< contactMatrix.length; i++) {
            for (int j=0; j<contactMatrix.length; j++) {
                System.out.print("pos " + i + " " + j + ": " + contactMatrix[i][j] + " ");
            }
            System.out.println();
        }
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
