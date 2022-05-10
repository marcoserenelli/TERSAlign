import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;

import fr.orsay.lri.varna.models.treealign.*;
import org.apache.commons.cli.*;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.io.PDBFileReader;

/**
 * MainComparator class interacting with the user through command line options.
 *
 * @author Luca Tesei
 *
 */
public class MainComparator {

    public static void main(String[] args) {
        // Use Apache Commons CLI 1.4
        // create Options object for Command Line Definition
        Options options = new Options();
        // define command line options
        Option o1 = new Option("sc","structcode",true,"Produce the structural RNA/Protein tree corresponding to the given structure by PDB code");
        o1.setArgName("input-file");
        options.addOption(o1);
        Option o2 = new Option("sf","structfile",true,"Produce the structural RNA/Protein tree corresponding to the given structure by PDB file");
        o2.setArgName("input-file");
        options.addOption(o2);
        Option o3 = new Option("ac","aligncode",true,"Align two given structures by PDB code producing alignment tree and distance");
        o3.setArgs(2);
        o3.setArgName("input-file1 input-file2");
        options.addOption(o3);
        Option o11 = new Option("af","alignfile",true,"Align two given structures by PDB file producing alignment tree and distance");
        o11.setArgs(2);
        o11.setArgName("input-file1 input-file2");
        options.addOption(o11);
        Option o4 = new Option("o","out",true,"Output result on the given file instead of standard output");
        o4.setArgName("output-file");
        options.addOption(o4);
        Option o5 = new Option("l","latexout",false,"Output in LaTeX format instead of linearised tree");
        options.addOption(o5);
        Option o6 = new Option("i","info",false,"Show license and other info");
        options.addOption(o6);
        Option o7 = new Option("h","help",false,"Show usage information");
        options.addOption(o7);
        Option o8 = new Option("d","outdist",false,"Output only distance, no alignment tree (works only with option -a)");
        options.addOption(o8);
        Option o9 = new Option("e","showscores",false,"Show current values of edit scores used for alignment");
        options.addOption(o9);
        Option o10 = new Option("n","useconffile",true,"Use the specified configuration file instead of the default one");
        o10.setArgName("conf-file");
        options.addOption(o10);

        // Parse command line
        HelpFormatter formatter = new HelpFormatter();
        CommandLineParser commandLineParser = new DefaultParser();
        CommandLine cmd = null;
        try {
            cmd = commandLineParser.parse(options, args);
        } catch (ParseException e) {
            // oops, something went wrong
            System.err.println("ERROR: Command Line parsing failed.  Reason: " + e.getMessage() + "\n");
            formatter.printHelp(CommandLineMessages.LAUNCH_COMMAND, CommandLineMessages.HEADER, options,
                    CommandLineMessages.USAGE_EXAMPLES + CommandLineMessages.COPYRIGHT
                            + CommandLineMessages.SHORT_NOTICE + CommandLineMessages.REPORT_TO,
                    true);
            System.exit(1);
        }

        // Manage Option h
        if (cmd.hasOption("h")) {
            formatter.printHelp(CommandLineMessages.LAUNCH_COMMAND, CommandLineMessages.HEADER, options,
                    CommandLineMessages.USAGE_EXAMPLES + CommandLineMessages.COPYRIGHT
                            + CommandLineMessages.SHORT_NOTICE + CommandLineMessages.REPORT_TO,
                    true);
            return;
        }

        // Manage Option i
        if (cmd.hasOption("i")) {
            Options optionsEmpty = new Options();
            formatter
                    .printHelp(CommandLineMessages.LAUNCH_COMMAND, "", optionsEmpty,
                            CommandLineMessages.COPYRIGHT + CommandLineMessages.LONG_NOTICE
                                    + CommandLineMessages.REPORT_TO + "\n\nUse option -h for full usage information",
                            true);
            return;
        }

        // Manage Option n
        String configurationFileName = ScoringFunction.DEFAULT_PROPERTY_FILE;
        if (cmd.hasOption("n")) {
            configurationFileName = cmd.getOptionValue("n");
        }

        // Manage Option e
        if (cmd.hasOption("e")) {
            ScoringFunction f = new ScoringFunction(configurationFileName);
            String scores = "ASPRAlign current costs:\n\n" + "Cost for Operator Insertion = "
                    + f.getInsertOperatorCost() + "\nCost for Operator Deletion = " + f.getDeleteOperatorCost()
                    + "\nCost for Operator Replacement with Operator = " + f.getReplaceOperatorCost()
                    + "\nCost for Hairpin Insertion = " + f.getInsertHairpinCost() + "\nCost for Hairpin Deletion = "
                    + f.getDeleteHairpinCost() + "\nCost for One Crossing Mismatch (Local Cost) = "
                    + f.getCrossingMismatchCost();
            System.out.println(scores);
            return;
        }

        // Manage option s
        if (cmd.hasOption("sc") || cmd.hasOption("sf")) {
            boolean filepath = cmd.hasOption("sf");
            TertiaryStructure tertiaryStructure = null;
            Structure struc;
            try {
                if(filepath){
                    String filename =  cmd.getOptionValue("sf");
                    PDBFileReader pdbreader = new PDBFileReader();
                    struc = pdbreader.getStructure(filename);
                } else {
                    struc = StructureIO.getStructure(cmd.getOptionValue("sc"));
                }
                tertiaryStructure = new TertiaryStructure(struc);
            } catch (Exception e) {
                System.err.println("ERROR:" + e.getMessage());
                System.exit(3);
            }
            // Construct the ASPRATtree
            Tree<String> t;
            StructuralTree tree = new StructuralTree(tertiaryStructure);
            // get the structural RNA/Protein tree
            t = tree.getStructuralRNATree();
            // Produce Output
            String output;
            if (cmd.hasOption("l"))
                // produce LaTeX
                output = TreeOutputter.toLatex(t);
            else
                // produce linearised tree
                output = TreeOutputter.treeToString(t);

            // Write Output on proper file or on standard output
            if (cmd.hasOption("o")) {
                String outputFile = cmd.getOptionValue("o");
                try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile, false))) {
                    writer.write(output);
                } catch (FileNotFoundException e) {
                    System.err.println("ERROR: Output file " + outputFile + " cannot be created.");
                    System.exit(3);
                } catch (IOException e) {
                    System.err.println("ERROR: Output file" + outputFile + " cannot be written.");
                    System.exit(3);
                }
            } else
                System.out.println(output);
            return;
        }

        // Manage Option a
        if (cmd.hasOption("ac") || cmd.hasOption("af")) {
            boolean filePath = cmd.hasOption("af");
            // variables for structural RNA trees
            Tree<String> t1;
            Tree<String> t2;
            // Parse the first input file for the secondary structure
            Structure struc;
            TertiaryStructure tertiaryStructure = null;
            try {
                if(filePath){
                    PDBFileReader pdbreader = new PDBFileReader();
                    struc = pdbreader.getStructure(cmd.getOptionValues("af")[0]);
                } else {
                    struc = StructureIO.getStructure(cmd.getOptionValues("ac")[0]);
                }
                tertiaryStructure = new TertiaryStructure(struc);
            } catch (Exception e) {
                System.err.println("ERROR:" + e.getMessage());
                System.exit(3);
            }
            // Construct structural RNA/Protein tree 1
            StructuralTree s1 = new StructuralTree(tertiaryStructure);
            t1 = s1.getStructuralRNATree();

            // Parse the second input file for the secondary structure
            Structure struc2;
            TertiaryStructure tertiaryStructure2 = null;
            try {
                if(filePath){
                    PDBFileReader pdbreader = new PDBFileReader();
                    struc2 = pdbreader.getStructure(cmd.getOptionValues("af")[1]);
                } else {
                    struc2 = StructureIO.getStructure(cmd.getOptionValues("ac")[1]);
                }
                tertiaryStructure2 = new TertiaryStructure(struc2);
            } catch (Exception e) {
                System.err.println("ERROR:" + e.getMessage());
                System.exit(3);
            }
            // Construct structural RNA/Protein tree 2
            StructuralTree s2 = new StructuralTree(tertiaryStructure2);
            t2 = s2.getStructuralRNATree();

            // Align t1 and t2, which contain two structural RNA trees
            AlignmentResult r = null;
            ScoringFunction f = new ScoringFunction(configurationFileName);
            try {
                r = new AlignmentResult(t1, t2, f);
            } catch (TreeAlignException e) {
                System.err.println("ERROR: Alignment Exception: " + e.getMessage());
                System.exit(4);
            }

            // Produce Output
            if (!cmd.hasOption("d")) {
                Tree<AlignedNode<String, String>> t = r.getAlignedTree();
                String output;
                if (cmd.hasOption("l"))
                    // produce LaTeX
                    output = TreeOutputter.toLatexAligned(t);
                else
                    // produce linearised tree
                    output = TreeOutputter.treeToStringAligned(t);
                // Write Output on proper file or on standard output
                if (cmd.hasOption("o")) {
                    String outputFile = cmd.getOptionValue("o");
                    try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile, false))) {
                        writer.write(output);
                    } catch (FileNotFoundException e) {
                        System.err.println("ERROR: Output file " + outputFile + " cannot be created.");
                        System.exit(3);
                    } catch (IOException e) {
                        System.err.println("ERROR: Output file" + outputFile + " cannot be written.");
                        System.exit(3);
                    }
                } else
                    System.out.println(output + "\n");
            }
            // Output distance
            System.out.println("Distance = " + r.getDistance());
            return;
        }

        // If no option is given, output usage
        System.err.println("No operation specified...");
        formatter.printHelp(CommandLineMessages.LAUNCH_COMMAND, CommandLineMessages.HEADER, options,
                CommandLineMessages.USAGE_EXAMPLES + CommandLineMessages.COPYRIGHT + CommandLineMessages.SHORT_NOTICE
                        + CommandLineMessages.REPORT_TO,
                true);
    }
}
