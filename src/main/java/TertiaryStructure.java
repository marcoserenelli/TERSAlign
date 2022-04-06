import org.biojava.nbio.structure.GroupType;
import org.biojava.nbio.structure.Structure;

public class TertiaryStructure {

    private Structure structure;

    public TertiaryStructure(Structure structure) {
        this.structure = structure;
    }

    public double[][] getDistanceMatrix(){
        return RepAtmsDistanceMatrixGenerator.calculateDistanceMatrix(this.structure);
    }

}
