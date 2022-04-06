import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.io.PDBFileReader;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

class RepAtmsDistanceMatrixGeneratorTest {

    @Test
    @DisplayName("Are values computed right?")
    void testCalculateDistanceMatrix(){
        //Load a PROTEIN file
        Structure struc = loadFile("3mge");

        double[][] resultMatrix = RepAtmsDistanceMatrixGenerator.calculateDistanceMatrix(struc);

        //Check matrix length
        assertEquals(204, resultMatrix.length, "Matrix length should be 204");
        //Test on the first two CA distance
        assertEquals(3.7951635010892466, resultMatrix[0][1], "Distance between the first two CA should be 3.7951...");

        //Test on two randoms CA distance in the middle
        assertEquals(5.029965009818657, resultMatrix[35][32], "Distance between two CA values in the middle should be 5.0299...");

        //Test on the last two CA distance
        assertEquals(3.799415876157805, resultMatrix[resultMatrix.length-2][resultMatrix.length-1], "Distance between the last two CA values should be 3.7994...");

        //Load an RNA file
        struc = loadFile("4gxy");
        resultMatrix = RepAtmsDistanceMatrixGenerator.calculateDistanceMatrix(struc);

        //Check matrix length
        assertEquals(162, resultMatrix.length, "Matrix length should be 162");

        //Test on the first two P distance
        assertEquals(6.002236000025325, resultMatrix[0][1], "Distance between the first two P should be 6.0022...");

        //Test on two randoms P distance in the middle
        assertEquals(17.053572382348516, resultMatrix[35][32], "Distance between two P values in the middle should be 17.0535...");

        //Test on the last two P distance
        assertEquals(6.141912731389141, resultMatrix[resultMatrix.length-2][resultMatrix.length-1], "Distance between the last two P values should be 6.1419...");
    }

    /**
     * Load a PDB files, returns a structure, file should be in ../resources
     * @param fileName name of the file to load
     * @return returns a structure
     */
    private Structure loadFile(String fileName){
        PDBFileReader pdbreader = new PDBFileReader();
        pdbreader.setPath("../resources");
        try{
            return pdbreader.getStructureById(fileName);
        } catch (Exception e){
            e.printStackTrace();
        }
        return null;
    }
}