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
            tertiaryStructure.setThreshold(8.0);
            tertiaryStructure.getBondList().forEach(p -> System.out.println("first: " + p.getFirst() + ", second: " + p.getSecond()));
            System.out.println(tertiaryStructure.getType());
        } catch (Exception e){
            e.printStackTrace();
        }

    }
}

