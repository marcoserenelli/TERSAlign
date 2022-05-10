import fr.orsay.lri.varna.models.treealign.*;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.contact.Pair;
import org.biojava.nbio.structure.*;
import static org.junit.jupiter.api.Assertions.*;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;

class StructuralTreeTest {

    @Test
    @DisplayName("Comparing expected and actual trees")
    void recBuildStructural() {
        ArrayList<Pair<Integer>> bonds = new ArrayList<>();
        //Simple tests
        // (1,4); (6,9); concat
        bonds.add(new Pair<>(1, 4));
        bonds.add(new Pair<>(6, 9));
        Tree<String> expectedTree = testConcat();
        //Check this expected tree with real tree using given bonds
        assertTrue(isEquals(expectedTree, bonds, 12), "Concat test failed");
        bonds.clear();


        // (1,5); (5,10); meet
        bonds.add(new Pair<>(1, 5));
        bonds.add(new Pair<>(5, 10));
        expectedTree = testMeet();
        //Check this expected tree with real tree using given bonds
        assertTrue(isEquals(expectedTree, bonds,11), "Meet test failed");
        bonds.clear();


        // (1,9); (3,6); nesting
        bonds.add(new Pair<>(1, 9));
        bonds.add(new Pair<>(3, 6));
        expectedTree = testNesting();
        //Check this expected tree with real tree using given bonds
        assertTrue(isEquals(expectedTree, bonds,10), "Nesting test failed");
        bonds.clear();

        // (1,5); (3,5); ending
        bonds.add(new Pair<>(1, 5));
        bonds.add(new Pair<>(3, 5));
        expectedTree = testEnding();
        //Check this expected tree with real tree using given bonds
        assertTrue(isEquals(expectedTree, bonds,7), "Ending test failed");
        bonds.clear();

        // (1,7); (1,5); starting
        bonds.add(new Pair<>(1, 7));
        bonds.add(new Pair<>(1, 5));
        expectedTree = testStarting();
        //Check this expected tree with real tree using given bonds
        assertTrue(isEquals(expectedTree, bonds,9), "Starting test failed");
        bonds.clear();

        // (1,6); (3,9); crossing
        bonds.add(new Pair<>(1, 6));
        bonds.add(new Pair<>(3, 9));
        expectedTree = testCrossing();
        //Check this expected tree with real tree using given bonds
        assertTrue(isEquals(expectedTree, bonds,10), "Crossing test failed");
        bonds.clear();

        // (1,6); (1,3); (3,6); diamond
        bonds.add(new Pair<>(1, 6));
        bonds.add(new Pair<>(1, 3));
        bonds.add(new Pair<>(3, 6));
        expectedTree = testDiamond();
        //Check this expected tree with real tree using given bonds
        assertTrue(isEquals(expectedTree, bonds,7), "Diamond test failed");
        bonds.clear();

        //More advanced tests

        // (1,5); (5,10); (12,15); testConcatMeet
        bonds.add(new Pair<>(1, 5));
        bonds.add(new Pair<>(5, 10));
        bonds.add(new Pair<>(12, 15));
        expectedTree = testConcatMeet();
        //Check this expected tree with real tree using given bonds
        assertTrue(isEquals(expectedTree, bonds,16), "ConcatMeet test failed");
        bonds.clear();

        //(1,3) (5,8) (10,13) (15,18) testFourConcat
        bonds.add(new Pair<>(1, 3));
        bonds.add(new Pair<>(5, 8));
        bonds.add(new Pair<>(10, 13));
        bonds.add(new Pair<>(15, 18));
        expectedTree = testThreeConcat();
        //Check this expected tree with real tree using given bonds
        assertTrue(isEquals(expectedTree, bonds,20), "FourConcat test failed");
        bonds.clear();

        //(1,13); (2,6); (4,8); (8,12); testNestingMeetCrossing
        bonds.add(new Pair<>(1, 13));
        bonds.add(new Pair<>(2, 6));
        bonds.add(new Pair<>(4, 8));
        bonds.add(new Pair<>(8, 12));
        expectedTree = testNestingMeetCrossing();
        //Check this expected tree with real tree using given bonds
        assertTrue(isEquals(expectedTree, bonds,14), "NestingMeetCrossing test failed");
        bonds.clear();

        // (1,5);(5,10); (15,20);(20;25); testMeetConcatMeet
        bonds.add(new Pair<>(1, 5));
        bonds.add(new Pair<>(5, 10));
        bonds.add(new Pair<>(15, 20));
        bonds.add(new Pair<>(20, 25));
        expectedTree = testMeetConcatMeet();
        //Check this expected tree with real tree using given bonds
        assertTrue(isEquals(expectedTree, bonds,26), "MeetConcatMeet test failed");
        bonds.clear();

        // (1,3); (3,6); (6,9); (12,15); (18,21); (21,24); testMeetConcatConcatMeetMeet
        bonds.add(new Pair<>(1, 3));
        bonds.add(new Pair<>(3, 6));
        bonds.add(new Pair<>(6, 9));
        bonds.add(new Pair<>(12, 15));
        bonds.add(new Pair<>(18, 21));
        bonds.add(new Pair<>(21, 24));
        expectedTree = testMeetConcatConcatMeetMeet();
        //Check this expected tree with real tree using given bonds
        assertTrue(isEquals(expectedTree, bonds,26), "MeetConcatConcatMeetMeet test failed");
        bonds.clear();

        //(1,7); (2,7); (4,8); (8,12); testMeetCrossEnding
        bonds.add(new Pair<>(1, 7));
        bonds.add(new Pair<>(2, 7));
        bonds.add(new Pair<>(4, 8));
        bonds.add(new Pair<>(8, 12));
        expectedTree = testMeetCrossEnding();
        //Check this expected tree with real tree using given bonds
        assertTrue(isEquals(expectedTree, bonds,13), "MeetCrossEnding test failed");
        bonds.clear();

        //(1,3); (1,12); (6,9); (7,9); (9,12); testMeetConcCrossNestDiamond
        bonds.add(new Pair<>(1, 3));
        bonds.add(new Pair<>(1, 12));
        bonds.add(new Pair<>(6, 9));
        bonds.add(new Pair<>(7, 9));
        bonds.add(new Pair<>(9, 12));
        expectedTree = testMeetEndConcatDiamond();
        //Check this expected tree with real tree using given bonds
        assertTrue(isEquals(expectedTree, bonds,13), "MeetConcCrossNestDiamond test failed");
        bonds.clear();
    }

    private Tree<String> testConcat() {
        // (1,4); (6,9); concat
        Tree<String> structuralTree = new Tree<>();

        structuralTree.setValue(Operators.CONCATENATION_LABEL);

        Tree<String> left = new Tree<>();
        Tree<String> right = new Tree<>();

        //update tree
        ArrayList<Tree<String>> concatChilds = new ArrayList<>();
        concatChilds.add(left);
        concatChilds.add(right);
        structuralTree.replaceChildrenListBy(concatChilds);

        right.setValue(Operators.HAIRPIN_LABEL + "(6,9)");
        left.setValue(Operators.HAIRPIN_LABEL + "(1,4)");

        return structuralTree;
    }

    private Tree<String> testMeet() {
        // (1,5); (5,10); meet
        Tree<String> structuralTree = new Tree<>();

        structuralTree.setValue(Operators.MEETING_LABEL);

        Tree<String> left = new Tree<>();
        Tree<String> right = new Tree<>();

        //update tree
        ArrayList<Tree<String>> meetChilds = new ArrayList<>();
        meetChilds.add(left);
        meetChilds.add(right);
        structuralTree.replaceChildrenListBy(meetChilds);

        left.setValue(Operators.HAIRPIN_LABEL + "(1,5)");
        right.setValue(Operators.HAIRPIN_LABEL + "(5,10)");

        return structuralTree;
    }

    private Tree<String> testNesting() {
        // (1,9); (3,6);
        Tree<String> structuralTree = new Tree<>();

        structuralTree.setValue(Operators.NESTING_LABEL);

        Tree<String> right = new Tree<>();
        Tree<String> left = new Tree<>();

        //update tree
        ArrayList<Tree<String>> nestingChilds = new ArrayList<>();
        nestingChilds.add(left);
        nestingChilds.add(right);
        structuralTree.replaceChildrenListBy(nestingChilds);

        right.setValue(Operators.HAIRPIN_LABEL +"(1,9)" );
        left.setValue(Operators.HAIRPIN_LABEL + "(3,6)" );
        return structuralTree;
    }

    private Tree<String> testEnding() {
        // (1,5); (3,5); ending
        Tree<String> structuralTree = new Tree<>();

        structuralTree.setValue(Operators.ENDING_LABEL);

        Tree<String> right = new Tree<>();
        Tree<String> left = new Tree<>();

        //update tree
        ArrayList<Tree<String>> endingChilds = new ArrayList<>();
        endingChilds.add(left);
        endingChilds.add(right);
        structuralTree.replaceChildrenListBy(endingChilds);

        right.setValue(Operators.HAIRPIN_LABEL + "(1,5)" );
        left.setValue(Operators.HAIRPIN_LABEL + "(3,5)");
        return structuralTree;
    }

    private Tree<String> testStarting() {
        // (1,7); (1,5); starting
        Tree<String> structuralTree = new Tree<>();

        structuralTree.setValue(Operators.STARTING_LABEL);

        Tree<String> right = new Tree<>();
        Tree<String> left = new Tree<>();

        //update tree
        ArrayList<Tree<String>> startingChilds = new ArrayList<>();
        startingChilds.add(left);
        startingChilds.add(right);
        structuralTree.replaceChildrenListBy(startingChilds);

        right.setValue(Operators.HAIRPIN_LABEL + "(1,7)");
        left.setValue(Operators.HAIRPIN_LABEL + "(1,5)");

        return structuralTree;
    }

    private Tree<String> testCrossing() {
        // (1,6); (3,9); crossing
        Tree<String> structuralTree = new Tree<>();

        structuralTree.setValue("(" + Operators.CROSSING_LABEL + ",1)");

        Tree<String> right = new Tree<>();
        Tree<String> left = new Tree<>();

        //update tree
        ArrayList<Tree<String>> crossingChilds = new ArrayList<>();
        crossingChilds.add(left);
        crossingChilds.add(right);
        structuralTree.replaceChildrenListBy(crossingChilds);

        right.setValue(Operators.HAIRPIN_LABEL + "(3,9)");
        left.setValue(Operators.HAIRPIN_LABEL + "(1,6)");

        return structuralTree;
    }

    private Tree<String> testDiamond() {
        // (1,6); (1,3); (3,6); diamond
        Tree<String> structuralTree = new Tree<>();

        structuralTree.setValue(Operators.DIAMOND_LABEL);

        Tree<String> left = new Tree<>();
        Tree<String> right = new Tree<>();

        //update tree
        ArrayList<Tree<String>> diamondChilds = new ArrayList<>();
        diamondChilds.add(left);
        diamondChilds.add(right);
        structuralTree.replaceChildrenListBy(diamondChilds);

        left.setValue(Operators.MEETING_LABEL);
        right.setValue(Operators.HAIRPIN_LABEL + "(1,6)");

        //2nd level

        Tree<String> left2 = new Tree<>();
        Tree<String> right2 = new Tree<>();

        //update tree
        ArrayList<Tree<String>> concatChilds2 = new ArrayList<>();
        concatChilds2.add(left2);
        concatChilds2.add(right2);
        left.replaceChildrenListBy(concatChilds2);

        right2.setValue(Operators.HAIRPIN_LABEL + "(3,6)");
        left2.setValue(Operators.HAIRPIN_LABEL + "(1,3)");

        return structuralTree;
    }

    private Tree<String> testNestingMeetCrossing() {
        //(1,13); (2,6); (4,8); (8,12);
        Tree<String> structuralTree = new Tree<>();
        structuralTree.setValue(Operators.NESTING_LABEL);

        Tree<String> left = new Tree<>();
        Tree<String> right = new Tree<>();

        //update tree
        ArrayList<Tree<String>> nestingChilds = new ArrayList<>();
        nestingChilds.add(left);
        nestingChilds.add(right);
        structuralTree.replaceChildrenListBy(nestingChilds);

        right.setValue(Operators.HAIRPIN_LABEL + "(1,13)");
        left.setValue(Operators.MEETING_LABEL);

        //2nd level
        Tree<String> left2 = new Tree<>();
        Tree<String> right2 = new Tree<>();

        //update tree
        ArrayList<Tree<String>> meetingChilds = new ArrayList<>();
        meetingChilds.add(left2);
        meetingChilds.add(right2);
        left.replaceChildrenListBy(meetingChilds);

        right2.setValue(Operators.HAIRPIN_LABEL + "(8,12)");
        left2.setValue("(" + Operators.CROSSING_LABEL + ",1)");

        //3rd level
        Tree<String> left3 = new Tree<>();
        Tree<String> right3 = new Tree<>();

        //update tree
        ArrayList<Tree<String>> meetingChilds2 = new ArrayList<>();
        meetingChilds2.add(left3);
        meetingChilds2.add(right3);
        left2.replaceChildrenListBy(meetingChilds2);

        right3.setValue(Operators.HAIRPIN_LABEL + "(4,8)");
        left3.setValue(Operators.HAIRPIN_LABEL + "(2,6)");

        return structuralTree;
    }

    private Tree<String> testMeetConcatConcatMeetMeet() {

        // (1,3); (3,6); (6,9); (12,15); (18,21); (21,24);
        Tree<String> structuralTree = new Tree<>();
        structuralTree.setValue(Operators.MEETING_LABEL);

        // create the node for building the left part
        Tree<String> left = new Tree<>();
        // create node for building the right part
        Tree<String> right = new Tree<>();

        // update tree
        ArrayList<Tree<String>> meetingChilds = new ArrayList<>();
        meetingChilds.add(left);
        meetingChilds.add(right);
        structuralTree.replaceChildrenListBy(meetingChilds);

        left.setValue(Operators.CONCATENATION_LABEL);
        right.setValue(Operators.HAIRPIN_LABEL + "(21,24)");

        // create the node for building the left part
        Tree<String> left2 = new Tree<>();
        // create node for building the right part
        Tree<String> right2 = new Tree<>();

        // update tree
        ArrayList<Tree<String>> concatChilds = new ArrayList<>();
        concatChilds.add(left2);
        concatChilds.add(right2);
        left.replaceChildrenListBy(concatChilds);

        left2.setValue(Operators.CONCATENATION_LABEL);
        right2.setValue(Operators.HAIRPIN_LABEL + "(18,21)");

        // create the node for building the left part
        Tree<String> left3 = new Tree<>();
        // create node for building the right part
        Tree<String> right3 = new Tree<>();

        // update tree
        ArrayList<Tree<String>> concatChilds2 = new ArrayList<>();
        concatChilds2.add(left3);
        concatChilds2.add(right3);
        left2.replaceChildrenListBy(concatChilds2);

        left3.setValue(Operators.MEETING_LABEL);
        right3.setValue(Operators.HAIRPIN_LABEL + "(12,15)");

        // create the node for building the left part
        Tree<String> left4 = new Tree<>();
        // create node for building the right part
        Tree<String> right4 = new Tree<>();

        // update tree
        ArrayList<Tree<String>> meetingChilds2 = new ArrayList<>();
        meetingChilds2.add(left4);
        meetingChilds2.add(right4);
        left3.replaceChildrenListBy(meetingChilds2);

        left4.setValue(Operators.MEETING_LABEL);
        right4.setValue(Operators.HAIRPIN_LABEL + "(6,9)");

        // create the node for building the left part
        Tree<String> left5 = new Tree<>();
        // create node for building the right part
        Tree<String> right5 = new Tree<>();

        // update tree
        ArrayList<Tree<String>> meetingChilds3 = new ArrayList<>();
        meetingChilds3.add(left5);
        meetingChilds3.add(right5);
        left4.replaceChildrenListBy(meetingChilds3);

        left5.setValue(Operators.HAIRPIN_LABEL + "(1,3)");
        right5.setValue(Operators.HAIRPIN_LABEL + "(3,6)");

        return structuralTree;
    }

    private Tree<String> testMeetConcatMeet() {
        Tree<String> structuralTree = new Tree<>();

        // (1,5);(5,10); (15,20);(20;25); first we encounter a meet, then a concat, then a meet
        structuralTree.setValue(Operators.MEETING_LABEL);

        // create the node for building the left part
        Tree<String> left = new Tree<>();
        // create node for building the right part
        Tree<String> right = new Tree<>();

        // update tree
        ArrayList<Tree<String>> meetingChilds = new ArrayList<>();
        meetingChilds.add(left);
        meetingChilds.add(right);
        structuralTree.replaceChildrenListBy(meetingChilds);

        right.setValue(Operators.HAIRPIN_LABEL + "(20,25)");
        left.setValue(Operators.CONCATENATION_LABEL);

        //2nd level
        // create the node for building the left part
        Tree<String> left2 = new Tree<>();
        // create node for building the right part
        Tree<String> right2 = new Tree<>();

        // update tree
        ArrayList<Tree<String>> concatChilds = new ArrayList<>();
        concatChilds.add(left2);
        concatChilds.add(right2);
        left.replaceChildrenListBy(concatChilds);

        right2.setValue(Operators.HAIRPIN_LABEL + "(15,20)");
        left2.setValue(Operators.MEETING_LABEL);

        //3rd level
        // create the node for building the left part
        Tree<String> left3 = new Tree<>();
        // create node for building the right part
        Tree<String> right3 = new Tree<>();

        // update tree
        ArrayList<Tree<String>> meetingChilds2 = new ArrayList<>();
        meetingChilds2.add(left3);
        meetingChilds2.add(right3);
        left2.replaceChildrenListBy(meetingChilds2);

        left3.setValue(Operators.HAIRPIN_LABEL + "(1,5)");
        right3.setValue(Operators.HAIRPIN_LABEL + "(5,10)");

        return structuralTree;
    }

    private Tree<String> testConcatMeet() {
        Tree<String> structuralTree = new Tree<>();

        // (1,5); (5,10); (12,15); first we encounter a concat, then a zero interval, then a meet
        structuralTree.setValue(Operators.CONCATENATION_LABEL);

        // create the node for building the left part
        Tree<String> left = new Tree<>();

        // create node for building the right part
        Tree<String> right = new Tree<>();

        // update tree
        ArrayList<Tree<String>> concatChilds = new ArrayList<>();
        concatChilds.add(left);
        concatChilds.add(right);
        structuralTree.replaceChildrenListBy(concatChilds);

        right.setValue(Operators.HAIRPIN_LABEL + "(12,15)");
        left.setValue(Operators.MEETING_LABEL);

        //2nd level
        Tree<String> left2 = new Tree<>();
        Tree<String> right2 = new Tree<>();

        // update tree
        ArrayList<Tree<String>> meetingChilds = new ArrayList<>();
        meetingChilds.add(left2);
        meetingChilds.add(right2);
        left.replaceChildrenListBy(meetingChilds);

        right2.setValue(Operators.HAIRPIN_LABEL + "(5,10)");
        left2.setValue(Operators.HAIRPIN_LABEL + "(1,5)");

        return structuralTree;
    }

    private Tree<String> testThreeConcat() {
        //(1,3) (5,8) (10,13) (15,18)
        Tree<String> structuralTree = new Tree<>();
        structuralTree.setValue(Operators.CONCATENATION_LABEL);
        Tree<String> left = new Tree<>();
        Tree<String> right = new Tree<>();
        ArrayList<Tree<String>> concatChilds = new ArrayList<>();
        concatChilds.add(left);
        concatChilds.add(right);
        structuralTree.replaceChildrenListBy(concatChilds);
        right.setValue(Operators.HAIRPIN_LABEL + "(15,18)");
        left.setValue(Operators.CONCATENATION_LABEL);
        //2nd level
        Tree<String> left2 = new Tree<>();
        Tree<String> right2 = new Tree<>();
        ArrayList<Tree<String>> concatChilds2 = new ArrayList<>();
        concatChilds2.add(left2);
        concatChilds2.add(right2);
        left.replaceChildrenListBy(concatChilds2);
        right2.setValue(Operators.HAIRPIN_LABEL + "(10,13)");
        left2.setValue(Operators.CONCATENATION_LABEL);
        //3rd level
        Tree<String> left3 = new Tree<>();
        Tree<String> right3 = new Tree<>();
        ArrayList<Tree<String>> concatChilds3 = new ArrayList<>();
        concatChilds3.add(left3);
        concatChilds3.add(right3);
        left2.replaceChildrenListBy(concatChilds3);
        right3.setValue(Operators.HAIRPIN_LABEL + "(5,8)");
        left3.setValue(Operators.HAIRPIN_LABEL + "(1,3)");

        return structuralTree;
    }

    private Tree<String> testMeetCrossEnding() {
        //(1,7); (2,7); (4,8); (8,12);
        Tree<String> structuralTree = new Tree<>();

        structuralTree.setValue(Operators.MEETING_LABEL);

        Tree<String> left = new Tree<>();
        Tree<String> right = new Tree<>();

        ArrayList<Tree<String>> meetingChilds = new ArrayList<>();
        meetingChilds.add(left);
        meetingChilds.add(right);
        structuralTree.replaceChildrenListBy(meetingChilds);

        right.setValue(Operators.HAIRPIN_LABEL + "(8,12)");
        left.setValue("(" + Operators.CROSSING_LABEL + ",2)");

        //2nd level
        Tree<String> left2 = new Tree<>();
        Tree<String> right2 = new Tree<>();

        ArrayList<Tree<String>> crossingChilds = new ArrayList<>();
        crossingChilds.add(left2);
        crossingChilds.add(right2);
        left.replaceChildrenListBy(crossingChilds);

        right2.setValue(Operators.HAIRPIN_LABEL + "(4,8)");
        left2.setValue(Operators.ENDING_LABEL);

        //3rd level
        Tree<String> left3 = new Tree<>();
        Tree<String> right3 = new Tree<>();

        ArrayList<Tree<String>> endingChilds = new ArrayList<>();
        endingChilds.add(left3);
        endingChilds.add(right3);
        left2.replaceChildrenListBy(endingChilds);

        right3.setValue(Operators.HAIRPIN_LABEL + "(1,7)");
        left3.setValue(Operators.HAIRPIN_LABEL + "(2,7)");

        return structuralTree;
    }

    private Tree<String> testMeetEndConcatDiamond() {
        //(1,2); (1,12); (6,9); (7,9); (9,12);
        Tree<String> structuralTree = new Tree<>();

        structuralTree.setValue(Operators.DIAMOND_LABEL);

        Tree<String> left = new Tree<>();
        Tree<String> right = new Tree<>();

        ArrayList<Tree<String>> diamondChilds = new ArrayList<>();
        diamondChilds.add(left);
        diamondChilds.add(right);
        structuralTree.replaceChildrenListBy(diamondChilds);

        right.setValue(Operators.HAIRPIN_LABEL + "(1,12)");
        left.setValue(Operators.MEETING_LABEL);

        //2nd level
        Tree<String> left2 = new Tree<>();
        Tree<String> right2 = new Tree<>();

        ArrayList<Tree<String>> meetChilds = new ArrayList<>();
        meetChilds.add(left2);
        meetChilds.add(right2);
        left.replaceChildrenListBy(meetChilds);

        right2.setValue(Operators.HAIRPIN_LABEL + "(9,12)");
        left2.setValue(Operators.CONCATENATION_LABEL);

        //3rd level
        Tree<String> left3 = new Tree<>();
        Tree<String> right3 = new Tree<>();

        ArrayList<Tree<String>> concatChilds = new ArrayList<>();
        concatChilds.add(left3);
        concatChilds.add(right3);
        left2.replaceChildrenListBy(concatChilds);

        right3.setValue(Operators.ENDING_LABEL);
        left3.setValue(Operators.HAIRPIN_LABEL + "(1,3)");

        //4th level
        Tree<String> left4 = new Tree<>();
        Tree<String> right4 = new Tree<>();

        ArrayList<Tree<String>> endingChilds = new ArrayList<>();
        endingChilds.add(left4);
        endingChilds.add(right4);
        right3.replaceChildrenListBy(endingChilds);

        right4.setValue(Operators.HAIRPIN_LABEL + "(6,9)");
        left4.setValue(Operators.HAIRPIN_LABEL + "(7,9)");


        return structuralTree;
    }


    /**
     * Iterate the trees and checks if they have the same values in each nodes
     * @param expectedTree the expected tree
     * @param bonds bonds for the actual tree
     * @return true if the trees are equals false otherwise
     */
    private boolean isEquals(Tree<String> expectedTree, ArrayList<Pair<Integer>> bonds, int sequenceLenght ) {

        StructuralTree structuralTree = CreateStructuralTreeWithCustomBonds(bonds,sequenceLenght);
        Tree<String> actualTree = structuralTree.getStructuralRNATree();

        //If the trees don't have the same number of nodes return false
        if(expectedTree.countNodes() != actualTree.countNodes()) {
            System.out.println("Trees length is different: expected " + expectedTree.countNodes() + " actual " + actualTree.countNodes());
            return false;
        }

        Iterator<Tree<String>> expectedIterator = expectedTree.iterator();
        Iterator<Tree<String>> actualIterator = actualTree.iterator();

        while(expectedIterator.hasNext() && actualIterator.hasNext()) {
            Tree<String> expectedNode = expectedIterator.next();
            Tree<String> actualNode = actualIterator.next();
            System.out.println(expectedNode.getValue() + " vs " + actualNode.getValue());
            if(!expectedNode.getValue().equals(actualNode.getValue())) {
                return false;
            }
        }

        return true;
    }

    /**
     * Creates a new structual tree by loading a random structure but replacing the current bonds with
     * the bonds passed as parameter
     * @param bonds bonds for the actual tree
     * @param sequenceLength sequence length of the actual tree
     * @return the new structual tree
     */
    private StructuralTree CreateStructuralTreeWithCustomBonds(ArrayList<Pair<Integer>> bonds, int sequenceLength) {
        Structure structure = null;
        try {
            structure = StructureIO.getStructure("3mge");
        } catch (IOException | StructureException e) {
            e.printStackTrace();
        }
        TertiaryStructure tertiaryStructure = new TertiaryStructure(structure);
        tertiaryStructure.setBondList(bonds);
        StructuralTree treeGenerator = new StructuralTree(tertiaryStructure);
        treeGenerator.setSequenceLength(sequenceLength);
        return treeGenerator;
    }


}