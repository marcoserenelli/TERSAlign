import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.io.PDBFileReader;

public class TertiaryStructuresComparator {

    public static void main(String[] args) {
        PDBFileReader pdbReader = new PDBFileReader();
        pdbReader.setPath("../../main/resources");
        String filename = "3mge";
        try{
            Structure structure = pdbReader.getStructureById(filename);
            TertiaryStructure tertiaryStructure = new TertiaryStructure(structure);
            System.out.println(tertiaryStructure.getType());
        } catch (Exception e){
            e.printStackTrace();
        }

    }
}

