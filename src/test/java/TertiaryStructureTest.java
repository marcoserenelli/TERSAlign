import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.io.PDBFileReader;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class TertiaryStructureTest {

    @Test
    @DisplayName("Distance matrix with default method")
    void testCalculateDistanceMatrix(){
        //Load a PROTEIN file
        Structure struc = loadFile("3mge");

        TertiaryStructure tertiaryStructure = new TertiaryStructure(struc);

        double[][] resultMatrix = tertiaryStructure.calculateDistanceMatrix(struc);

        //Check matrix length
        assertEquals(204, resultMatrix.length, "Matrix length should be 204");

        //Test the first value with himself
        assertEquals(0,resultMatrix[0][0], "Distance should be 0");

        //Test the last value with himself
        assertEquals(0,resultMatrix[resultMatrix.length-1][resultMatrix.length-1], "Distance should be 0");

        //Test on the first two CA distance
        assertEquals(3.7951635010892466, resultMatrix[0][1], "Distance between the first two CA should be 3.7951...");

        //Test on two randoms CA distance in the middle
        assertEquals(5.029965009818657, resultMatrix[35][32], "Distance between two CA values in the middle should be 5.0299...");

        //Test on the last two CA distance
        assertEquals(3.799415876157805, resultMatrix[resultMatrix.length-2][resultMatrix.length-1], "Distance between the last two CA values should be 3.7994...");

        //Load an RNA file
        struc = loadFile("4gxy");
        tertiaryStructure = new TertiaryStructure(struc);

        resultMatrix = tertiaryStructure.calculateDistanceMatrix(struc);

        //Check matrix length
        assertEquals(162, resultMatrix.length, "Matrix length should be 162");

        //Test the first value with himself
        assertEquals(0,resultMatrix[0][0], "Distance should be 0");

        //Test the last value with himself
        assertEquals(0,resultMatrix[resultMatrix.length-1][resultMatrix.length-1], "Distance should be 0");

        //Test on the first two P distance
        assertEquals(6.002236000025325, resultMatrix[0][1], "Distance between the first two P should be 6.0022...");

        //Test on two randoms P distance in the middle
        assertEquals(17.053572382348516, resultMatrix[35][32], "Distance between two P values in the middle should be 17.0535...");

        //Test on the last two P distance
        assertEquals(6.141912731389141, resultMatrix[resultMatrix.length-2][resultMatrix.length-1], "Distance between the last two P values should be 6.1419...");
    }

    @Test
    @DisplayName("Distance matrix with center of mass")
    void testCalculateDistanceMatrixCenterOfMass(){

        //Load a PROTEIN file
        Structure struc = loadFile("3mge");

        TertiaryStructure tertiaryStructure = new TertiaryStructure(struc);

        double[][] resultMatrix = tertiaryStructure.calculateDistanceMatrixCenterOfMass(struc);

        //Check matrix length
        assertEquals(346, resultMatrix.length, "Matrix length should be 346");

        //Test the first value with himself
        assertEquals(0,resultMatrix[0][0], "Distance should be 0");

        //Test the last value with himself
        assertEquals(Double.MAX_VALUE,resultMatrix[resultMatrix.length-1][resultMatrix.length-1], "Distance should be MAX VALUE");

        //Test on the first two distance
        assertEquals(4.620100901228067, resultMatrix[0][1], "Distance should be 4.6201...");

        //Test on two randoms distances in the middle
        assertEquals(5.236362921223696, resultMatrix[35][32], "Distance should be 5.2363...");

        //Test on the last two CA distance
        assertEquals(Double.MAX_VALUE, resultMatrix[resultMatrix.length-2][resultMatrix.length-1], "Distance should be MAX VALUE...");


        //Load an RNA file
        struc = loadFile("4gxy");
        tertiaryStructure = new TertiaryStructure(struc);

        resultMatrix = tertiaryStructure.calculateDistanceMatrixCenterOfMass(struc);

        //Check matrix length
        assertEquals(175, resultMatrix.length, "Matrix length should be 175");

        //Test the first value with himself
        assertEquals(Double.MAX_VALUE,resultMatrix[0][0], "Distance should be MAX VALUE");

        //Test the last value with himself
        assertEquals(Double.MAX_VALUE,resultMatrix[resultMatrix.length-1][resultMatrix.length-1], "Distance should be MAX VALUE");

        //Test on the first two distance
        assertEquals(Double.MAX_VALUE, resultMatrix[0][1], "Distance should be MAX VALUE");

        //Test on two randoms distance in the middle
        assertEquals(14.3414881620807, resultMatrix[35][32], "Distance should be 14.3414...");

        assertEquals(30.881515689238928, resultMatrix[160][32], "Distance should be 30.881515689238928");

        //Test on the last two P distance
        assertEquals(Double.MAX_VALUE, resultMatrix[resultMatrix.length-2][resultMatrix.length-1], "Distance should be MAX VALUE...");
    }

    @Test
    @DisplayName("Contact matrix using default distance matrix")
    void testGetContactMatrixDefault(){

        //Load a PROTEIN file
        Structure struc = loadFile("3mge");
        TertiaryStructure tertiaryStructure = new TertiaryStructure(struc);
        tertiaryStructure.setThreshold(12);

        boolean[][] actualContactMatrix = tertiaryStructure.getContactMatrixDefault();

        //USING 12 AS THRESHOLD

        //Check if matrix length is right
        assertEquals(actualContactMatrix.length, tertiaryStructure.getDistanceMatrixDefault().length, "contact matrix length should be the same as distance matrix length");

        //First value should be equal to himself.
        assertTrue(actualContactMatrix[0][0], "first value should always be equal to himself");

        //Last value should be equal to himself.
        assertTrue(actualContactMatrix[actualContactMatrix.length-1][actualContactMatrix.length-1], "last value should always be equal to himself");

        //Should be true, distance is 3.7951...
        assertTrue(actualContactMatrix[0][1], "first with second value should be true, distance is 3.7951...");

        //Should be true, distance is 3.7994
        assertTrue(actualContactMatrix[actualContactMatrix.length-1][actualContactMatrix.length-2], "last two values should be true, distance is 3.7994...");

        //Should be true, distance is 5.0299...
        assertTrue(actualContactMatrix[32][35], "32 - 35 value should be true, distance is 5.0299...");

        //Should be false, distance is 13.3787...
        assertFalse(actualContactMatrix[7][3], "7 - 3 value should be false, distance is 13.3787...");

        //Changing threshold to 14 should make 7 - 3 true
        tertiaryStructure.setThreshold(14);
        //Re-calculate contact matrix with new threshold
        actualContactMatrix = tertiaryStructure.getContactMatrixDefault();
        //Should be true
        assertTrue(actualContactMatrix[7][3], "7 - 3 should be true, with threshold 14, value is 13.3787...");

        //Changing threshold to 4 should make 32 - 35 false
        tertiaryStructure.setThreshold(4);
        //Re-calculate contact matrix with new threshold
        actualContactMatrix = tertiaryStructure.getContactMatrixDefault();
        //Should be false
        assertFalse(actualContactMatrix[32][35], "32 - 35 should be false with threshold 4, value is 5.0299...");

        //Changing threshold to 51
        tertiaryStructure.setThreshold(51);
        //Re-calculate contact matrix with new threshold
        actualContactMatrix = tertiaryStructure.getContactMatrixDefault();
        //Should be false
        assertFalse(actualContactMatrix[190][2], "190 - 2 should be false with threshold 51, value is 51.8775...");


        //Load an RNA file
        struc = loadFile("4gxy");
        tertiaryStructure = new TertiaryStructure(struc);
        tertiaryStructure.setThreshold(12);
        actualContactMatrix = tertiaryStructure.getContactMatrixDefault();

        tertiaryStructure.printContactMatrix(actualContactMatrix);
        System.out.println(actualContactMatrix.length);
        tertiaryStructure.printDistanceMatrix(tertiaryStructure.getDistanceMatrixDefault());


        //USING 12 AS THRESHOLD

        //Check if matrix length is right
        assertEquals(actualContactMatrix.length, tertiaryStructure.getContactMatrixDefault().length, "contact matrix length should be the same as distance matrix length");

        //First value should be equal to himself.
        assertTrue(actualContactMatrix[0][0], "first value should always be equal to himself");

        //Last value should be equal to himself.
        assertTrue(actualContactMatrix[actualContactMatrix.length-1][actualContactMatrix.length-1], "last value should always be equal to himself");

        //Should be true, distance is 4.6201...
        assertTrue(actualContactMatrix[0][1], "first with second value should be true, distance is 6.0022...");

        //Should be true, distance is 6.1419
        assertTrue(actualContactMatrix[actualContactMatrix.length-1][actualContactMatrix.length-2], "last two values should be true, distance is 6.1419...");

        //Should be false, distance is 17.0535...
        assertFalse(actualContactMatrix[32][35], "32 - 35 value should be true, distance is 17.0535...");

        //Should be false, distance is 19.8870...
        assertFalse(actualContactMatrix[7][3], "7 - 3 value should be false, distance is 19.8870...");

        //Changing threshold to 20 should make 7 - 3 true
        tertiaryStructure.setThreshold(20);
        //Re-calculate contact matrix with new threshold
        actualContactMatrix = tertiaryStructure.getContactMatrixDefault();
        //Should be true
        assertTrue(actualContactMatrix[7][3], "7 - 3 should be true, with threshold 20, value is 19.8870...");

        //Changing threshold to 4 should make 32 - 35 false
        tertiaryStructure.setThreshold(4);
        //Re-calculate contact matrix with new threshold
        actualContactMatrix = tertiaryStructure.getContactMatrixDefault();
        //Should be false
        assertFalse(actualContactMatrix[32][35], "32 - 35 should be false with threshold 4, value is 17.0535...");

        //Changing threshold to 41
        tertiaryStructure.setThreshold(41);
        //Re-calculate contact matrix with new threshold
        actualContactMatrix = tertiaryStructure.getContactMatrixDefault();
        //Should be false
        assertFalse(actualContactMatrix[132][35], "132 - 35 should be false with threshold 41, value is 41.3573...");
    }

    @Test
    @DisplayName("Contact matrix using center of mass distance matrix")
    void testGetContactMatrixDistanceCenterOfMass(){
        //Load a PROTEIN file
        Structure struc = loadFile("3mge");
        TertiaryStructure tertiaryStructure = new TertiaryStructure(struc);
        tertiaryStructure.setThreshold(12);

        boolean[][] actualContactMatrix = tertiaryStructure.getContactMatrixCenterOfMass();


        //USING 12 AS THRESHOLD

        //Check if matrix length is right
        assertEquals(actualContactMatrix.length, tertiaryStructure.getContactMatrixCenterOfMass().length, "contact matrix length should be the same as distance matrix length");

        //First value should be equal to himself.
        assertTrue(actualContactMatrix[0][0], "first value should be equal to himself");

        //Last value should not be to himself.
        assertFalse(actualContactMatrix[actualContactMatrix.length-1][actualContactMatrix.length-1], "last value should not be equal to himself");

        //Should be true, distance is 4.6201...
        assertTrue(actualContactMatrix[0][1], "first with second value should be true, distance is 4.6201...");

        //Should be true, distance is 3.7994
        assertFalse(actualContactMatrix[actualContactMatrix.length-1][actualContactMatrix.length-2], "last two values should be true, distance is MAX DOUBLE VALUE...");

        //Should be true, distance is 5.2363...
        assertTrue(actualContactMatrix[32][35], "32 - 35 value should be true, distance is 5.2363...");

        //Should be false, distance is 13.6747...
        assertFalse(actualContactMatrix[7][3], "7 - 3 value should be false, distance is 13.6747...");

        //Changing threshold to 14 should make 7 - 3 true
        tertiaryStructure.setThreshold(14);
        //Re-calculate contact matrix with new threshold
        actualContactMatrix = tertiaryStructure.getContactMatrixCenterOfMass();
        //Should be true
        assertTrue(actualContactMatrix[7][3], "7 - 3 should be true, with threshold 14, value is 13.6747...");

        //Changing threshold to 4 should make 32 - 35 false
        tertiaryStructure.setThreshold(4);
        //Re-calculate contact matrix with new threshold
        actualContactMatrix = tertiaryStructure.getContactMatrixCenterOfMass();
        //Should be false
        assertFalse(actualContactMatrix[32][35], "32 - 35 should be false with threshold 4, value is 5.2363...");

        //Changing threshold to 51
        tertiaryStructure.setThreshold(51);
        //Re-calculate contact matrix with new threshold
        actualContactMatrix = tertiaryStructure.getContactMatrixCenterOfMass();
        //Should be false
        assertFalse(actualContactMatrix[190][2], "190 - 2 should be false with threshold 51, value is 51.5449...");


        //Load an RNA file
        struc = loadFile("4gxy");
        tertiaryStructure = new TertiaryStructure(struc);
        tertiaryStructure.setThreshold(12);
        actualContactMatrix = tertiaryStructure.getContactMatrixCenterOfMass();

        tertiaryStructure.printContactMatrix(actualContactMatrix);
        tertiaryStructure.printDistanceMatrix(tertiaryStructure.getDistanceMatrixCenterOfMass());

        //USING 12 AS THRESHOLD

        //Check if matrix length is right
        assertEquals(actualContactMatrix.length, tertiaryStructure.getContactMatrixCenterOfMass().length, "contact matrix length should be the same as distance matrix length");

        //First value should be not be equal to himself.
        assertFalse(actualContactMatrix[0][0], "first value should not be equal to himself");

        //Last value should not be equal to himself.
        assertFalse(actualContactMatrix[actualContactMatrix.length-1][actualContactMatrix.length-1], "last value should not be equal to himself");

        //Should be false, distance is MAX DOBULE...
        assertFalse(actualContactMatrix[0][1], "first with second value should be false, distance is MAX DOUBLE...");

        //Should be false, distance is MAX DOUBLE
        assertFalse(actualContactMatrix[actualContactMatrix.length-1][actualContactMatrix.length-2], "last two values should be false, distance is MAX VALUE...");

        //Should be false, distance is 14.3414...
        assertFalse(actualContactMatrix[32][35], "32 - 35 value should be true, distance is 14.3414...");

        //Should be false, distance is 17.4826...
        assertFalse(actualContactMatrix[7][3], "7 - 3 value should be false, distance is 17.4826...");

        //Changing threshold to 20 should make 7 - 3 true
        tertiaryStructure.setThreshold(20);
        //Re-calculate contact matrix with new threshold
        actualContactMatrix = tertiaryStructure.getContactMatrixCenterOfMass();
        //Should be true
        assertTrue(actualContactMatrix[7][3], "7 - 3 should be true, with threshold 20, value is 17.4826...");

        //Changing threshold to 4 should make 32 - 35 false
        tertiaryStructure.setThreshold(4);
        //Re-calculate contact matrix with new threshold
        actualContactMatrix = tertiaryStructure.getContactMatrixCenterOfMass();
        //Should be false
        assertFalse(actualContactMatrix[32][35], "32 - 35 should be false with threshold 4, value is 14.3414...");

        //Changing threshold to 49
        tertiaryStructure.setThreshold(49);
        //Re-calculate contact matrix with new threshold
        actualContactMatrix = tertiaryStructure.getContactMatrixCenterOfMass();
        //Should be false
        assertTrue(actualContactMatrix[132][35], "132 - 35 should be true with threshold 49, value is 48.6214...");

    }

    @Test
    @DisplayName("Is type right?")
    void testGetType(){
        //Load a PROTEIN file
        Structure struc = loadFile("3mge");
        TertiaryStructure tertiaryStructure = new TertiaryStructure(struc);

        //Should be amminoacid
        assertEquals(GroupType.AMINOACID, tertiaryStructure.getType());

        //Load an RNA file
        struc = loadFile("4gxy");
        tertiaryStructure = new TertiaryStructure(struc);

        //Should be nucleotide
        assertEquals(GroupType.NUCLEOTIDE, tertiaryStructure.getType());
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