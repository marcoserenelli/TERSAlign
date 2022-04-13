import org.biojava.nbio.structure.*;

public class TertiaryStructuresComparator {

    public static void main(String[] args) {
        try{
            Structure structure = StructureIO.getStructure("7jxh");
            TertiaryStructure tertiaryStructure = new TertiaryStructure(structure);
            System.out.println(tertiaryStructure.getType());
            System.out.println(tertiaryStructure.getBondList());
            System.out.println(tertiaryStructure.getSecondaryStructure());
            System.out.println(tertiaryStructure.getThreshold());
            System.out.println(tertiaryStructure.getSequence());
        } catch (Exception e){
            e.printStackTrace();
        }

    }
}

