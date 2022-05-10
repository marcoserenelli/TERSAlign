import org.biojava.nbio.structure.*;

public class TertiaryStructuresComparator {

    public static void main(String[] args) {
        try{
            Structure structure = StructureIO.getStructure("3mge");
            TertiaryStructure tertiaryStructure = new TertiaryStructure(structure);
            System.out.println(tertiaryStructure.getBondList());
            StructuralTree treeBuilder = new StructuralTree(tertiaryStructure);
            System.out.println(TreeOutputter.treeToString(treeBuilder.getStructuralRNATree()));
        } catch (Exception e){
            e.printStackTrace();
        }

    }
}

