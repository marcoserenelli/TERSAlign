import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.io.PDBFileReader;

public class TertiaryStructuresComparator {

    public static void main(String[] args) {
        PDBFileReader pdbreader = new PDBFileReader();
        pdbreader.setPath("../../main/resources");
        String filename = "4gxy";
        try{
            Structure struc = pdbreader.getStructureById(filename);
            RepAtmsDistanceMatrixGenerator.printAmminoacidDistanceMatrix(RepAtmsDistanceMatrixGenerator.calculateDistanceMatrix(struc));
            //FirstDistanceMatrixGenerator.printAmminoacidDistanceMatrix(FirstDistanceMatrixGenerator.calcuateDistanceMatrix(struc));
        } catch (Exception e){
            e.printStackTrace();
        }
    }
}

