import fr.orsay.lri.varna.models.treealign.*;
import org.biojava.nbio.structure.contact.Pair;

import java.util.ArrayList;
import java.util.Arrays;

public class StructuralTree {

    private final TertiaryStructure tertiaryStructure;
    private Tree<String> structuralTree;
    private int[] c;
    private ArrayList<Integer>[] p;
    private int[] m;

    public StructuralTree(TertiaryStructure tertiaryStructure){
        this.tertiaryStructure = tertiaryStructure;
        this.structuralTree = null;
    }

    /**
     * @return the original RNA secondary structure
     */
    public TertiaryStructure getTertiaryStructure() {
        return this.tertiaryStructure;
    }

    public Tree<String> getStructuralRNATree() {
        if (this.structuralTree == null)
            buildStructural();
        return this.structuralTree;
    }

    private void buildStructural() {
        // create array for counting openings and closings of loops, by default it is
        // initialized to zeroes, starts at 1, position 0 not used
        int[] c = new int[this.tertiaryStructure.getSequence().length()];

        // initialize the counting array
        initc(c);

        // init indexes for later recursion call
        int l = 1; // left index, starts at 1, position 0 not used
        int r = this.tertiaryStructure.getSequence().length(); // right index

        // move l to the start of the structure
        while (c[l] == 0)
            l++;

        // move r to the tail of the structure
        while (c[r] == 0)
            r--;
        r++; // last closing loop has 0 count, but belongs to the loop, so it's not part of
        // the tail

        // the largest pseudoloop is now identified by the interval [l,r]
        assert c[l] == 1 && c[r] == 0 : "Largest pseudoloop at [" + l + "," + r + "]\nCounting array: " + Arrays.toString(c);

        // create an empty list of zero intervals to detect concatenations
        ArrayList<Interval> zi = new ArrayList<>();

        // find zero intervals in the outermost pseudoloop, if any
        detectZeroIntervals(c, zi, l, r);

        // create a copy of the list of weak bonds of the structure, ensuring that the
        // copied list conserves ordering
        ArrayList <Pair<Integer>> bonds = new ArrayList<>(this.tertiaryStructure.getBondList());

        // create the root node of the structural RNA tree
        Tree<String> t = new Tree<>();

        // start the recursive construction of the structural RNA Tree on the node ct
        recBuildStructural(t, bonds, c, this.p, zi, l, r);

        // assign to the root of this tree
        this.structuralTree = t;
    }

    private void recBuildStructural(Tree<String> ct, ArrayList<Pair<Integer>> bonds, int[] c, ArrayList<Integer>[] p ,ArrayList<Interval> zi, int l, int r){

        // checks
        assert c[l] == 1 && c[r] == 0 : "Pseudoloop bounds error while parsing at [" + l + "," + r
                + "]\nCounting array: " + Arrays.toString(c);
        Pair<Integer> lastBond = bonds.get(bonds.size() - 1);
        assert lastBond.getFirst() == r : "Mismatch among indexes of secondary structure and of "
                + "determined loop in WeakBond: (" + l + "," + r + ") vs (" + lastBond.getSecond() + ","
                + lastBond.getFirst() + ")";

        //If there are no zero intervals, the hairpin can't be a concat
        if(zi.isEmpty()){
            if(!detectMeet(l,r)){
                //If the C array has a value greater than 1, in the arriving point of the hairpin it means that the hairpin is a cross
                if(c[p[r].get(0)] > 1 ){
                    //cross case
                    ct.setValue(Operators.CROSSING_LABEL);
                }
                //If the P size in the r position, is more than one then it can be an ending or a diamond
                if(p[r].size() > 1){
                    Pair<Integer> longestBond = new Pair<Integer>(Integer.MAX_VALUE,Integer.MAX_VALUE);
                    for(int currentBond : p[r]){
                        if(p[currentBond].size() > 1){
                            longestBond = currentBond - r < longestBond.getFirst() ? new Pair<Integer>(p[currentBond].get(0),r) : longestBond;

                            if(longestBond.getFirst() >= l){
                                //ending case
                                ct.setValue(Operators.ENDING_LABEL);
                            }
                        } else {
                            //diamond case
                            ct.setValue(Operators.DIAMOND_LABEL);
                        }
                    }
                } else {
                    //If the P size in the arriving point of the hairpin is greater than 1, then the hairpin it's a starting
                    if(p[p[r].get(0)].size() > 1){
                        //starting case
                        ct.setValue(Operators.STARTING_LABEL);
                    } else {
                        //If the P in the arriving point of the hairpin is equals to l, it means that the hairpin is a nest
                        if(p[r].get(0) == l){
                            //nest case
                            ct.setValue(Operators.NESTING_LABEL);
                        }
                    }
                }
            }
            //If meets are detected than it's a meet
            else {
                //meet case
                ct.setValue(Operators.MEETING_LABEL);
            }
        }
        //If there are zero intervals, the hairpin is a concat
        else {
            //concat case
            ct.setValue(Operators.CONCATENATION_LABEL);
        }
    }

    /**
     * Returns true if a meet is find in the array between l and r indexes
     * @param l left index
     * @param r right index
     * @return true if a meet is founded, false otherwise
     */
    private boolean detectMeet(int l, int r ){
        for(int i = l; i <= r; i++){
            if(this.m[i] == 1){
                return true;
            }
        }
        return false;
    }

    private void splitBonds(ArrayList<Pair<Integer>> bonds, int ll, int lr, int rl, int rr, ArrayList<Pair<Integer>> lbonds,
                            ArrayList<Pair<Integer>> rbonds) {
        for (Pair<Integer> b : bonds) {
            if (ll <= b.getSecond() && b.getFirst() <= lr)
                lbonds.add(b);
            else if (rl <= b.getSecond() && b.getFirst() <= rr)
                rbonds.add(b);
            else
                assert false : "Error in splitting bonds: weak bond (" + b.getSecond() + "," + b.getFirst()
                        + ") not in range [" + ll + "," + lr + "] and not in range [" + rl + "," + rr + "]";
        }
        assert lbonds.size() > 0 : "Error in splitting bonds: left bonds list empty";
        assert rbonds.size() > 0 : "Error in splitting bonds: right bonds list empty";
    }

    private int determineNumberOfCrossings(ArrayList<Pair<Integer>> bonds) {
        int n = 0;
        int lastBondLeftIndex = bonds.get(bonds.size() - 1).getSecond();
        for (int i = bonds.size() - 2; i >= 0; i--) {
            Pair<Integer> b = bonds.get(i);
            if (b.getSecond() < lastBondLeftIndex && lastBondLeftIndex < b.getFirst())
                n++;
        }
        assert n > 0 : "Crossing number equal to zero!";
        return n;
    }

    /*
     * Finds all the (possibly empty) zero intervals inside the pseudoloop [l,r]. A
     * zero interval is a section of the primary sequence that separates two
     * concatenated pseudoloops. It is called zero interval because in the counting
     * array the count goes to zero before the end of the pseudoloop. A zero
     * interval always starts at the first position after the count went to zero. It
     * stops at the first position in which the counting raises to one again. If
     * these positions coincide the zero interval is empty. All the found intervals
     * are put into the list zi, which is assumed to be empty at the calling time.
     */
    private void detectZeroIntervals(int[] c, ArrayList<Interval> zi, int l, int r) {
        assert l < r : "Empty pseudoloop while detecting zero intervals at [" + l + "," + r + "]";
        assert c[l] == 1 && c[r] == 0 : "Pseudoloop bounds error while detecting zero intervals at [" + l + "," + r
                + "]\nCounting array: " + Arrays.toString(c);
        assert zi.isEmpty() : "Not empty zero interval list while detecting zero intervals: " + zi;
        int i = l;
        do {
            // search for the next zero interval
            while (c[i] != 0)
                i++;
            if (i == r)
                break; // reached end of the interval, stop searching more zero intervals
            i++;
            // determine start of the zero interval
            int start = i;
            // search for the end of the zero interval
            while (c[i] == 0)
                i++;
            // determine the stop of the zero interva, if start == stop, the zero interval
            // is empty
            int stop = i;
            // create the interval and add it to the list
            Interval interval = new Interval(start, stop);
            zi.add(interval);
        } while (true);
    }

    private void initc(int[] c) {
        int count = 0;
        for (int i = 1; i <= this.tertiaryStructure.getSequence().length(); i++) {
            if (loop_start(i)) {
                // a loop has just started in the arc annotated sequence, so count is
                // incremented
                count++;
                c[i] = count;
            } else if (loop_stop(i)) {
                // a loop has just stopped in the arc annotated sequence, so count is
                // decremented
                count--;
                c[i] = count;
            } else {// no star/stop of a loop, the value of count is kept
                if(loop_start(i) && loop_stop(i))
                    this.m[i] = this.p[i].size();
                c[i] = count;
            }
        }
        assert count == 0 : "Value of count after initialization of counting array: " + count + "\nCounting array: "
                + Arrays.toString(c);
    }

    private boolean isPartner(int indexToCheck, int partnerToCheck){
        for(Integer currentPartner : this.p[indexToCheck])
            if(currentPartner == partnerToCheck)
                return true;
        return false;
    }

    @SuppressWarnings("unchecked")
    private void initp(){
        this.p = new ArrayList[500];
        for(Pair<Integer> currentBond : this.tertiaryStructure.getBondList()) {
            this.p[currentBond.getFirst()] = new ArrayList<>();
            this.p[currentBond.getFirst()].add(currentBond.getSecond());
            this.p[currentBond.getSecond()] = new ArrayList<>();
            this.p[currentBond.getSecond()].add(currentBond.getFirst());
        }
    }

    /*
     * Tells if at position i there is the starting of an hairpin loop of the
     * secondary structure represented by the original arc annotated sequence.
     */
    private boolean loop_start(int i){
        if(!(this.p[i] == null))
            for(Integer partner : this.p[i])
                if(partner > i)
                    return true;
        return false;
    }

    /*
     * Tells if at position i there is the ending of an hairpin loop of the
     * secondary structure represented by the original arc annotated sequence.
     */
    private boolean loop_stop(int i) {
        if(!(this.p[i] == null))
            for(Integer partner : this.p[i])
                if(partner < i)
                    return true;
        return false;
    }

    /*
     * Service class for holding zero intervals.
     */
    protected class Interval {
        private int i;
        private int j;

        protected Interval(int i, int j) {
            assert i <= j : "Interval [" + i + "," + j + "]";
            this.i = i;
            this.j = j;
        }

        /*
         * (non-Javadoc)
         *
         * @see java.lang.Object#hashCode()
         */
        @Override
        public int hashCode() {
            final int prime = 31;
            int result = 1;
            result = prime * result + getOuterType().hashCode();
            result = prime * result + i;
            result = prime * result + j;
            return result;
        }

        /*
         * A zero interval is empty if i == j
         */
        protected boolean isEmpty() {
            return this.i == this.j;
        }

        /*
         * (non-Javadoc)
         *
         * @see java.lang.Object#equals(java.lang.Object)
         */
        @Override
        public boolean equals(Object obj) {
            if (this == obj)
                return true;
            if (obj == null)
                return false;
            if (getClass() != obj.getClass())
                return false;
            Interval other = (Interval) obj;
            if (!getOuterType().equals(other.getOuterType()))
                return false;
            if (i != other.i)
                return false;
            return j == other.j;
        }

        private StructuralTree getOuterType() {
            return StructuralTree.this;
        }

        /*
         * (non-Javadoc)
         *
         * @see java.lang.Object#toString()
         */
        @Override
        public String toString() {
            return "Interval [i=" + i + ", j=" + j + "]";
        }

    }


}

