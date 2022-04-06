import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.contact.AtomContact;
import org.biojava.nbio.structure.contact.AtomContactSet;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.PDBFileReader;
import org.biojava.nbio.structure.secstruc.*;

import java.io.IOException;
import java.util.List;

public class TertiaryStructuresComparator {

    public static void main(String[] args) {
        PDBFileReader pdbreader = new PDBFileReader();
        pdbreader.setPath("../../main/resources");
        String filename = "3mge";
        try{
            Structure struc = pdbreader.getStructureById(filename);
            RepAtmsDistanceMatrixGenerator.printAmminoacidDistanceMatrix(RepAtmsDistanceMatrixGenerator.calculateDistanceMatrixCenterOfMass(struc));

            /*
            AtomCache cache = new AtomCache();
            StructureIO.setAtomCache(cache);

            Chain chain = struc.getChain("A");

            // we want contacts between Calpha atoms only
            String[] atoms = {"CA"};

            AtomContactSet contacts = StructureTools.getAtomsInContact(chain, atoms, 4.5);
            System.out.println("Total number of CA-CA contacts: "+contacts.size());
            */

            /*
            for(Chain c: struc.getChains()){
                AtomContactSet contacts = StructureTools.getAtomsInContact(c, atoms, 4.5);
                System.out.println("Total number of CA-CA contacts: "+contacts.size());
            }*/

            /*
            FOR SECONDARY STRUCT PRINTING
            SecStrucCalc secStrucCalc = new SecStrucCalc();
            secStrucCalc.calculate(struc, true);
            System.out.println(secStrucCalc.printDSSP());
             */



        } catch (Exception e){
            e.printStackTrace();
        }

    }
}

