import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.contact.AtomContactSet;
import java.util.ArrayList;

public class Test {

    public static void testMethod(Structure structure){
        System.out.println(structure.getChains(0).get(0).getAtomGroup(1).getAtom(0).getBonds());
        String[] atoms = {"P"};
        ArrayList<AtomContactSet>contactSetsList = new ArrayList<>();
        for(Chain currentChain : structure.getChains())
            contactSetsList.add(StructureTools.getAtomsInContact(currentChain, atoms, 4.0));
        int totalCounter = 0;
        for(AtomContactSet currentSet : contactSetsList)
            totalCounter+=currentSet.size();
        System.out.println(totalCounter);
        System.out.println("------------");
    }
}
