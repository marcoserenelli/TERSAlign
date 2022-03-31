import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.io.PDBFileReader;

public class TertiaryStructuresComparator {

    public static void main(String[] args) {
        String filename =  "/C:\\Users\\Filippo\\Downloads\\4gxy.pdb" ;
        PDBFileReader pdbreader = new PDBFileReader();
        try{
            Structure struc = pdbreader.getStructure(filename);
            DistanceMatrixGenerator.printAmminoacidDistanceMatrix(DistanceMatrixGenerator.calcuateDistanceMatrix(struc));
        } catch (Exception e){
            e.printStackTrace();
        }
    }
}

