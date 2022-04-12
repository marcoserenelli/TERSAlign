import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.GroupType;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.secstruc.SecStrucCalc;
import org.biojava.nbio.structure.secstruc.SecStrucTools;
import org.jgrapht.alg.util.Pair;

import java.util.ArrayList;
import java.util.Locale;

public class TertiaryStructure {

    private final Structure structure;
    private double threshold; //Value between 4.5 and 12 ångström
    private SecStrucCalc secondaryStructure;
    private ArrayList<Pair<Integer, Integer>> bondList;
    private boolean[][] contactMatrix;
    private double[][] distanceMatrix;
    private String distanceMatrixCalculationMethod;

    public TertiaryStructure(Structure structure) {
        this.structure = structure;
        this.threshold = 8;
        this.secondaryStructure = null;
        this.bondList = null;
        this.contactMatrix = null;
        this.distanceMatrix = null;
        this.distanceMatrixCalculationMethod = "default";
    }

    /**
<<<<<<< HEAD
     * Returns a list containing the indexes of bonded nucleotides/aminos, that is all nucleotides/aminos closer than the specified threshold.
     * @return bond list
     */
    public ArrayList<Pair<Integer, Integer>> getBondList(){
        if(this.bondList == null)
            calculateBondList();
        return this.bondList;
    }

    public void calculateBondList() {
        boolean[][]contactMap = this.getContactMatrix();
=======
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
>>>>>>> master
        ArrayList<Pair<Integer, Integer>>bondList = new ArrayList<>();
        int colCount = 0;
        for(int i=0; i<contactMap.length; i++) {
            for (int j = colCount; j < contactMap.length; j++)
                if (contactMap[i][j]) {
                    bondList.add(new Pair<>(i, j));
                }
            colCount++;
        }
        this.bondList = bondList;
    }

    /**
     * Returns a boolean matrix, values are true if their distance (taken from default calculation)
     * is less than threshold value.
     * @return boolean contact matrix
     */
<<<<<<< HEAD
    public boolean[][] getContactMatrix(){
        if(this.contactMatrix == null)
            calculateContactMatrix();
        return this.contactMatrix;
=======
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
>>>>>>> master
    }

    private void calculateContactMatrix(){
        double[][] distanceMatrix = getDistanceMatrix();
        boolean[][] contactMatrix = new boolean[distanceMatrix.length][distanceMatrix.length];
        for (int i=0; i<distanceMatrix.length; i++) {
            for (int j = 0; j < distanceMatrix.length; j++) {
                contactMatrix[i][j] = distanceMatrix[i][j] <= this.threshold;
            }
        }
        this.contactMatrix = contactMatrix;
    }


    /**
     * Calculates the structure's distance matrix, considering aminos / nucleotide's center of mass or taking distance between P/CA atoms as comparison method,
     * depending on distance matrix calculation method
     * @return distance matrix
     */
    public double[][] getDistanceMatrix(){
        if(this.distanceMatrix == null) {
            if(this.distanceMatrixCalculationMethod.equals("default"))
                this.calculateDistanceMatrixDefault();
            else if (this.distanceMatrixCalculationMethod.equals("centerofmass"))
                this.calculateDistanceMatrixCenterOfMass();
        }
        return this.distanceMatrix;
    }

    private void calculateDistanceMatrixCenterOfMass(){
        int groupsNumber = getNonHetatmGroupsCounter(this.structure);
        double[][] distanceMatrix = new double[groupsNumber][groupsNumber];
        int moleculeCount = 0;
        for(Chain currentChain: this.structure.getChains())
            for (Group currentMolecule : currentChain.getAtomGroups()) {
                if (currentMolecule.getType() != GroupType.HETATM) {
                    int comparedMoleculeCount = 0;
                    for (Chain comparisonChain : this.structure.getChains()) {
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
        this.distanceMatrix = distanceMatrix;
    }

    private void calculateDistanceMatrixDefault(){
        Atom[] representativeAtomsArray = StructureTools.getRepresentativeAtomArray(this.structure);
        double[][] distanceMatrix = new double[representativeAtomsArray.length][representativeAtomsArray.length];
        for(int i=0; i<representativeAtomsArray.length; i++)
            for(int j=0; j<representativeAtomsArray.length; j++)
                distanceMatrix[i][j] = Calc.getDistance(representativeAtomsArray[i], representativeAtomsArray[j]);
        this.distanceMatrix = distanceMatrix;
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
     */
    public void printDistanceMatrix(){
        for (int i=0; i<this.distanceMatrix.length; i++) {
            for (int j=0; j<this.distanceMatrix.length; j++) {
                System.out.print("pos " + i + " " + j + ": " + this.distanceMatrix[i][j] + " ");
            }
            System.out.println();
        }
    }

    /**
     * Print the contact matrix
     */
    public void printContactMatrix(){
        for (int i=0; i<this.contactMatrix.length; i++) {
            for (int j=0; j<this.contactMatrix.length; j++) {
                System.out.print("pos " + i + " " + j + ": " + this.contactMatrix[i][j] + " ");
            }
            System.out.println();
        }
    }

    /**
     * Returns the structure's type as string
     * @return structure type (RNA/Protein)
     */
    public String getType(){
        for(Chain c : this.structure.getChains()){
            if(c.getPredominantGroupType() == GroupType.AMINOACID){
                return "PROTEIN";
            }
            else if(c.getPredominantGroupType() == GroupType.NUCLEOTIDE){
                return "RNA";
            }
        }
        return null;
    }

    /**
     * Returns the secondary structure predicted from the structure associated to this object
     * @return secondary structure
     */
    public SecStrucCalc getSecondaryStructure(){
        if(this.secondaryStructure == null)
            calculateSecondaryStructure();
        System.out.println(this.secondaryStructure.printDSSP());
        SecStrucTools.getSecStrucInfo(structure);
        return this.secondaryStructure;
    }

    private void calculateSecondaryStructure(){
        try {
            SecStrucCalc secStrucCalc = new SecStrucCalc();
            secStrucCalc.calculate(this.structure, true);
            this.secondaryStructure = secStrucCalc;
        }catch (Exception e){
            e.printStackTrace();
        }
    }

    public double getThreshold() {
        return threshold;
    }

    public void setThreshold(double threshold) {
        this.threshold = threshold;
        this.calculateContactMatrix();
    }

    public Structure getStructure() {
        return structure;
    }

    public void setDistanceMatrixCalculationMethod(String calculationMethod){
        if(calculationMethod.toLowerCase(Locale.ROOT).equals("default") || calculationMethod.toLowerCase(Locale.ROOT).equals("centerofmass"))
            this.distanceMatrixCalculationMethod = calculationMethod.toLowerCase(Locale.ROOT);
    }

}
