import org.biojava.nbio.structure.*;

public class TertiaryStructuresComparator {

    public static void main(String[] args) {
        try{
            Structure structure = StructureIO.getStructure("3mge");
            TertiaryStructure tertiaryStructure = new TertiaryStructure(structure);
            System.out.println(tertiaryStructure.getType());
            System.out.println(tertiaryStructure.getBondList());
            System.out.println(tertiaryStructure.getBondList().get(1));
            System.out.println(tertiaryStructure.getSecondaryStructure());
            System.out.println(tertiaryStructure.getThreshold());
            System.out.println(tertiaryStructure.getSequence());
            System.out.println(tertiaryStructure.getStructure());
            /*
            tertiaryStructure.setDistanceMatrixCalculationMethod("Default");
            tertiaryStructure.printDistanceMatrixToCSV();
            tertiaryStructure.setDistanceMatrixCalculationMethod("CenterOfMass");
            tertiaryStructure.printDistanceMatrixToCSV();
            tertiaryStructure.setDistanceMatrixCalculationMethod("Default");
            tertiaryStructure.printContactMatrixToCSV();
            tertiaryStructure.setDistanceMatrixCalculationMethod("CenterOfMass");
            tertiaryStructure.printContactMatrixToCSV(); */
        } catch (Exception e){
            e.printStackTrace();
        }

    }
}

