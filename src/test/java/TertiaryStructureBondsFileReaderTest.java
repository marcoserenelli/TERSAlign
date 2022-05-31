import static org.junit.jupiter.api.Assertions.*;

import it.unicam.cs.bdslab.tersaling.TertiaryStructureBondsFileReader;
import org.biojava.nbio.structure.contact.Pair;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

import java.io.IOException;
import java.util.ArrayList;

/**
 * Test class for the RNA/Protein's bonds file reader
 *
 * @author Marco Serenelli
 *
 */
class TertiaryStructureBondsFileReaderTest {


    @Test
    @DisplayName("Parsing bondlist from file")
    void testReadBondsList() throws IOException {
        //1
        ArrayList<Pair<Integer>> actualBonds = TertiaryStructureBondsFileReader.readBondsList("src/test/resources/bondsTest1.txt");
        //expectedBonds
        ArrayList<Pair<Integer>> expectedBonds = new ArrayList<>();
        expectedBonds.add(new Pair<>(0, 2));
        expectedBonds.add(new Pair<>(0, 3));
        expectedBonds.add(new Pair<>(1, 2));
        expectedBonds.add(new Pair<>(1, 3));
        expectedBonds.add(new Pair<>(2, 1));
        assertEquals(actualBonds, expectedBonds, "1st comparison failed");
        actualBonds.clear();
        expectedBonds.clear();

        //2
        actualBonds = TertiaryStructureBondsFileReader.readBondsList("src/test/resources/bondsTest2.txt");
        //expectedBonds
        expectedBonds.add(new Pair<>(1, 15));
        expectedBonds.add(new Pair<>(3, 13));
        expectedBonds.add(new Pair<>(4, 7));
        expectedBonds.add(new Pair<>(10, 18));
        expectedBonds.add(new Pair<>(20, 44));
        expectedBonds.add(new Pair<>(21, 26));
        expectedBonds.add(new Pair<>(22, 24));
        expectedBonds.add(new Pair<>(28, 35));
        expectedBonds.add(new Pair<>(29, 33));
        expectedBonds.add(new Pair<>(31, 34));
        expectedBonds.add(new Pair<>(36, 47));
        expectedBonds.add(new Pair<>(41, 50));
        expectedBonds.add(new Pair<>(53, 59));
        expectedBonds.add(new Pair<>(55, 61));
        assertEquals(actualBonds, expectedBonds, "2nd comparison failed");
        actualBonds.clear();
        expectedBonds.clear();

        //3
        //expected an error as you can't create two equals bonds
        assertThrows(IllegalArgumentException.class, () -> TertiaryStructureBondsFileReader.readBondsList("src/test/resources/bondsTest3.txt"), "3th comparison failed");
        actualBonds.clear();
        expectedBonds.clear();

        //4
        actualBonds = TertiaryStructureBondsFileReader.readBondsList("src/test/resources/bondsTest4.txt");
        //expectedBonds
        expectedBonds.add(new Pair<>(1,4));
        expectedBonds.add(new Pair<>(1,8));
        expectedBonds.add(new Pair<>(1,12));
        expectedBonds.add(new Pair<>(8,10));
        expectedBonds.add(new Pair<>(8,12));
        expectedBonds.add(new Pair<>(12,16));
        expectedBonds.add(new Pair<>(17,22));
        expectedBonds.add(new Pair<>(18,20));
        expectedBonds.add(new Pair<>(19,21));
        assertEquals(actualBonds, expectedBonds, "4th comparison failed");
        actualBonds.clear();
        expectedBonds.clear();

        //5
        actualBonds = TertiaryStructureBondsFileReader.readBondsList("src/test/resources/bondsTest5.txt");
        //expectedBonds
        expectedBonds.add(new Pair<>(1,4));
        expectedBonds.add(new Pair<>(1,7));
        expectedBonds.add(new Pair<>(7,10));
        expectedBonds.add(new Pair<>(7,12));
        expectedBonds.add(new Pair<>(12,17));
        expectedBonds.add(new Pair<>(19,26));
        expectedBonds.add(new Pair<>(20,23));
        expectedBonds.add(new Pair<>(22,25));
        assertEquals(actualBonds, expectedBonds, "5th comparison failed");
        actualBonds.clear();
        expectedBonds.clear();

        //6
        actualBonds = TertiaryStructureBondsFileReader.readBondsList("src/test/resources/bondsTest6.txt");
        //expectedBonds
        expectedBonds.add(new Pair<>(1,6));
        expectedBonds.add(new Pair<>(1,12));
        expectedBonds.add(new Pair<>(6,12));
        expectedBonds.add(new Pair<>(12,16));
        expectedBonds.add(new Pair<>(17,21));
        expectedBonds.add(new Pair<>(19,24));
        assertEquals(actualBonds, expectedBonds, "6th comparison failed");
        actualBonds.clear();
        expectedBonds.clear();

        //7
        actualBonds = TertiaryStructureBondsFileReader.readBondsList("src/test/resources/bondsTest7.txt");
        //expectedBonds
        expectedBonds.add(new Pair<>(1,6));
        expectedBonds.add(new Pair<>(1,12));
        expectedBonds.add(new Pair<>(6,12));
        expectedBonds.add(new Pair<>(7,17));
        expectedBonds.add(new Pair<>(12,17));
        expectedBonds.add(new Pair<>(18,22));
        expectedBonds.add(new Pair<>(18,24));
        expectedBonds.add(new Pair<>(20,24));
        assertEquals(actualBonds, expectedBonds, "7th comparison failed");
        actualBonds.clear();
        expectedBonds.clear();
    }
}