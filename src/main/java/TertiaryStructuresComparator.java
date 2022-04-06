import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.io.PDBFileReader;

public class TertiaryStructuresComparator {

    public static void main(String[] args) {
        PDBFileReader pdbReader = new PDBFileReader();
        pdbReader.setPath("../../main/resources");
        String filename = "4gxy";
        try{
            Structure structure = pdbReader.getStructureById(filename);

            TertiaryStructure tertiaryStructure = new TertiaryStructure(structure, 50);

            tertiaryStructure.setThreshold(12);
            double[][] distanceMatrix = tertiaryStructure.calculateDistanceMatrixCenterOfMass(structure);
            boolean[][] contactMatrix = tertiaryStructure.getContactMatrix();
            GroupType groupType = tertiaryStructure.getType();

        } catch (Exception e){
            e.printStackTrace();
        }

    }
}

