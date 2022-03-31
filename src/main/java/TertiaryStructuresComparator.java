import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.io.PDBFileReader;

public class TertiaryStructuresComparator {

    public static void main(String[] args) {
        String filename =  "/C:\\Users\\Filippo\\Downloads\\4hhb.pdb" ;
        PDBFileReader pdbreader = new PDBFileReader();
        try{
            Structure struc = pdbreader.getStructure(filename);
            System.out.println(struc.getChains(0).get(0));
            printAmminoacidDistanceMatrix(calcuateAmminoacidDistanceMatrix(struc));
        } catch (Exception e){
            e.printStackTrace();
        }
    }

    private static Double[][] calcuateAmminoacidDistanceMatrix(Structure struc){
        int totalSequenceLength = 0;
        for(Chain currentChain : struc.getChains()) {
            for(Group currentGroup : currentChain.getAtomGroups()){
                if(!(currentGroup.getType() == GroupType.HETATM))
                    totalSequenceLength++;
            }
        }
        Double[][] amminoAcidDistancesMatrix = new Double[totalSequenceLength][totalSequenceLength];
        int amminoCount = 0;
        for(Chain currentChain: struc.getChains()) {
                Atom amminoacidAlphaCarbon;
                for (Group currentAmminoacid : currentChain.getAtomGroups()) {

                    if(!(currentAmminoacid.getType() == GroupType.HETATM)) {
                        //è sicuramente proteina o nucleo
                        if (currentChain.isNucleicAcid() || currentChain.isProtein()) {
                            if (currentChain.isProtein())
                                amminoacidAlphaCarbon = currentAmminoacid.getAtom("CA");
                            else {
                                amminoacidAlphaCarbon = currentAmminoacid.getAtom("P");
                            }
                            int comparedAmminoCount = 0;
                            for (Chain comparisonChain : struc.getChains()) {
                                for (Group comparisonAmminoacid : comparisonChain.getAtomGroups()) {
                                    if(!(comparisonAmminoacid.getType() == GroupType.HETATM)) {
                                        Atom comparedAmminoacidAlphaCarbon;
                                        if (currentChain.isNucleicAcid() || currentChain.isProtein()) {
                                            if (comparisonChain.isProtein()) {
                                                comparedAmminoacidAlphaCarbon = comparisonAmminoacid.getAtom("CA");
                                                if (comparedAmminoacidAlphaCarbon == null)
                                                    System.out.println(comparisonAmminoacid + " amino");
                                            } else {
                                                comparedAmminoacidAlphaCarbon = comparisonAmminoacid.getAtom("P");
                                                if (comparedAmminoacidAlphaCarbon == null)
                                                    System.out.println(comparisonAmminoacid + " nucleo");
                                            }
                                            amminoAcidDistancesMatrix[amminoCount][comparedAmminoCount] = Calc.getDistance(amminoacidAlphaCarbon, comparedAmminoacidAlphaCarbon);
                                            comparedAmminoCount++;
                                        }
                                    }
                                }
                            }
                        }
                        amminoCount++;
                    }
                }
            }
        return amminoAcidDistancesMatrix;
    }

    private static void printAmminoacidDistanceMatrix(Double[][] distanceMatrix){
        for (int i=0; i<distanceMatrix.length; i++) {
            for (int j=0; j<distanceMatrix.length; j++) {
                System.out.print("pos " + i + " " + j + ": " + distanceMatrix[i][j] + " ");
            }
            System.out.println();
        }
    }

}

