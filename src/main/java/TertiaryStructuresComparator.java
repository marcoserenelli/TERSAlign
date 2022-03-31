import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.io.PDBFileReader;

public class TertiaryStructuresComparator {

    public static void main(String[] args) {
        String filename =  "/C:\\Users\\Filippo\\Downloads\\4hhb.pdb" ;
        PDBFileReader pdbreader = new PDBFileReader();
        try{
            Structure struc = pdbreader.getStructure(filename);
            printAmminoacidDistanceMatrix(calcuateAmminoacidDistanceMatrix(struc));
        } catch (Exception e){
            e.printStackTrace();
        }
    }

    private static Double[][] calcuateAmminoacidDistanceMatrix(Structure struc){
        int totalSequenceLength = 0;
        for(Chain currentChain : struc.getChains()) 
           totalSequenceLength += currentChain.getSeqResLength();
        Double[][] amminoAcidDistancesMatrix = new Double[totalSequenceLength][totalSequenceLength];
        int amminoCount = 0;
        for(Chain currentChain: struc.getChains()){
            Atom amminoacidAlphaCarbon;
            for (Group currentAmminoacid : currentChain.getAtomGroups()) {
                int comparedAmminoCount = 0;
                for(Chain comparisonChain : struc.getChains()) {
                    amminoacidAlphaCarbon = currentAmminoacid.getAtom("CA");
                    Atom comparedAmminoacidAlphaCarbon;
                    for (Group comparisonAmminoacid : comparisonChain.getAtomGroups()) {
                        if (currentChain.isProtein() && comparisonChain.isProtein()) {
                            comparedAmminoacidAlphaCarbon = comparisonAmminoacid.getAtom("CA");
                            amminoAcidDistancesMatrix[amminoCount][comparedAmminoCount] = Calc.getDistance(amminoacidAlphaCarbon, comparedAmminoacidAlphaCarbon);
                            comparedAmminoCount++;
                        }
                    }
                }
                amminoCount++;
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

