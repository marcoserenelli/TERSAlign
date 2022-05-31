package it.unicam.cs.bdslab.tersaling;

import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.GroupType;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.contact.Pair;
import org.biojava.nbio.structure.secstruc.SecStrucCalc;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Locale;

/**
 * Representation of an RNA/Protein structure, including its secondary structure and the methods to extract a bond list from its tertiary structure
 *
 * @author Filippo Lampa, Marco Serenelli
 */
public class TertiaryStructure {

    private final Structure structure;
    private double threshold; //Value between 4.5 and 12 ångström
    private SecStrucCalc secondaryStructure;
    private ArrayList<Pair<Integer>> bondList;
    private boolean[][] contactMatrix;
    private double[][] distanceMatrix;
    private String distanceMatrixCalculationMethod;
    private String sequence;

    /**
     * Creates a new TertiaryStructure from a PDB file's structure
     * @param structure the structure extracted from the PDB file
     */
    public TertiaryStructure(Structure structure) {
        this.structure = structure;
        this.sequence = null;
        this.threshold = 8;
        this.secondaryStructure = null;
        this.bondList = null;
        this.contactMatrix = null;
        this.distanceMatrix = null;
        this.distanceMatrixCalculationMethod = "default";
    }

    /**
     * Returns a list containing the indexes of bonded nucleotides/aminos, that is all nucleotides/aminos closer than the specified threshold as represented
     * within the contact map.
     * @return bond list
     */
    public ArrayList<Pair<Integer>> getBondList(){
        if(this.bondList == null)
            calculateBondList();
        return this.bondList;
    }

    private void calculateBondList() {
        boolean[][]contactMap = this.getContactMatrix();
        ArrayList<Pair<Integer>>bondList = new ArrayList<>();
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
    public boolean[][] getContactMatrix(){
        if(this.contactMatrix == null)
            calculateContactMatrix();
        return this.contactMatrix;
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
        if(this.distanceMatrixCalculationMethod.equals("default")){
            this.calculateDistanceMatrixDefault();
        }
        else if (this.distanceMatrixCalculationMethod.equals("centerofmass")){
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

    public String getSequence(){
        StringBuilder builder = new StringBuilder();
        for(Chain currentChain : this.structure.getChains())
            for(Group currentGroup : currentChain.getAtomGroups())
                if(currentGroup.isAminoAcid() || currentGroup.isNucleotide())
                    builder.append(StructureTools.get1LetterCode(currentGroup.getPDBName()));
        return builder.toString();
    }

    public void setDistanceMatrixCalculationMethod(String calculationMethod){
        if(calculationMethod.toLowerCase(Locale.ROOT).equals("default") || calculationMethod.toLowerCase(Locale.ROOT).equals("centerofmass"))
            this.distanceMatrixCalculationMethod = calculationMethod.toLowerCase(Locale.ROOT);
    }

    /**
     * Prints the distance matrix to a csv file
     */
    public void printDistanceMatrixToCSV(){
        try {
            double[][] distanceMatrix = this.getDistanceMatrix();
            FileWriter writer;
            if(this.distanceMatrixCalculationMethod.equals("default"))
                writer = new FileWriter("src/main/resources/DefaultDistanceMatrix.csv");
            else
                writer = new FileWriter("src/main/resources/DistanceMatrixCenterOfMass.csv");
            for (int i = 0; i < distanceMatrix.length; i++) {
                for (int j = 0; j < distanceMatrix.length; j++) {
                    writer.append(String.valueOf(i)).append(" ").append(String.valueOf(j)).append(":").append(" ").append(String.valueOf(distanceMatrix[i][j])).append("   ");
                }
                writer.append("\n");
                writer.append("\n");
            }
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * Prints the contact matrix to a csv file
     */
    public void printContactMatrixToCSV(){
        try {
            boolean[][] contactMatrix = this.getContactMatrix();
            FileWriter writer;
            writer = new FileWriter("src/main/resources/ContactMatrix.csv");
            for (int i = 0; i < contactMatrix.length; i++) {
                for (int j = 0; j < contactMatrix.length; j++) {
                    writer.append(String.valueOf(i)).append(" ").append(String.valueOf(j)).append(":").append(" ").append(String.valueOf(contactMatrix[i][j])).append("  ");
                }
                writer.append("\n");
            }
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * FOR TESTS PURPOSES
     * Replace current bonds list with a new one.
     * @param bondList new bond list
     */
    public void setBondList(ArrayList<Pair<Integer>> bondList) {
        this.bondList = bondList;
    }

}
