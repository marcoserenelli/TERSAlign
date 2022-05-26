package it.unicam.cs.bdslab.tersaling;

import org.antlr.v4.runtime.CharStream;
import org.antlr.v4.runtime.CharStreams;
import org.antlr.v4.runtime.CommonTokenStream;
import org.antlr.v4.runtime.tree.ParseTreeWalker;
import org.biojava.nbio.structure.contact.Pair;

import java.io.IOException;
import java.util.ArrayList;


public class TertiaryStructureBondsFileReader {

    /**
     * Use ANTLR 4 and the grammar defined for TertiaryStructure to read a
     * bonds list from a file.
     *
     * @param filename       the name of the file to read
     */
    public static ArrayList<Pair<Integer>> readBondsList(String filename) throws IOException {
        CharStream input = CharStreams.fromFileName(filename);
        // create a lexer that feeds off of input CharStream
        TertiaryStructureBondsLexer lexer = new TertiaryStructureBondsLexer(input);
        // create a buffer of tokens pulled from the lexer
        CommonTokenStream tokens = new CommonTokenStream(lexer);
        // create a parser that feeds off the tokens buffer
        TertiaryStructureBondsParser bondsParser = new TertiaryStructureBondsParser(tokens);
        // remove default error listeners
        bondsParser.removeErrorListeners();
        // begin parsing at bonds rule
        TertiaryStructureBondsParser.BondsContext tree = bondsParser.bonds();
        // Create a generic parse tree walker that can trigger callbacks
        ParseTreeWalker walker = new ParseTreeWalker();
        // Create the specialised listener for bonds
        TertiaryStructureBondsConstructor constructor = new TertiaryStructureBondsConstructor();
        // Walk the tree created during the parse, trigger callbacks
        walker.walk(constructor, tree);
        // Get the parsed secondary structure
        return constructor.getBondsList();
    }

}
