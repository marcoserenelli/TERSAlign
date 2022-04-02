import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.StructureAlignment;
import org.biojava.nbio.structure.align.StructureAlignmentFactory;
import org.biojava.nbio.structure.align.ce.CeMain;
import org.biojava.nbio.structure.align.ce.CeParameters;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.model.AfpChainWriter;

public class Aligner {

    public static AFPChain alignStructures(Structure struc1, Structure struc2){

      try {

          StructureAlignment algorithm  = StructureAlignmentFactory.getAlgorithm(CeMain.algorithmName);

          Atom[] ca1 = StructureTools.getAtomCAArray(struc1);
          Atom[] ca2 = StructureTools.getAtomCAArray(struc2);

          // get default parameters
          CeParameters params = new CeParameters();

          // add more print
          params.setShowAFPRanges(true);

          // set the maximum gap size to unlimited
          params.setMaxGapSize(-1);

          // The results are stored in an AFPChain object
          AFPChain afpChain = algorithm.align(ca1,ca2,params);

          afpChain.setName1(struc1.getName());
          afpChain.setName2(struc2.getName());

           // show a nice summary print`
          System.out.println(AfpChainWriter.toWebSiteDisplay(afpChain, ca1, ca2));

          // print rotation matrices`
          System.out.println(afpChain.toRotMat());
          //System.out.println(afpChain.toCE(ca1, ca2));

          // print XML representation`
          //System.out.println(AFPChainXMLConverter.toXML(afpChain,ca1,ca2));
          return afpChain;

      } catch (Exception e) {
          e.printStackTrace();
          return null;
      }
  }

}
