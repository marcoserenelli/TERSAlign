import fr.orsay.lri.varna.models.treealign.*;
import org.biojava.nbio.structure.contact.Pair;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Objects;

public class StructuralTree {

    private final TertiaryStructure tertiaryStructure;
    private Tree<String> structuralTree;
    private int sequenceLength;

    public StructuralTree(TertiaryStructure tertiaryStructure){
        this.tertiaryStructure = tertiaryStructure;
        this.sequenceLength = this.tertiaryStructure.getSequence().length();
        this.structuralTree = null;
    }

    public TertiaryStructure getTertiaryStructure(){
        return this.tertiaryStructure;
    }

    public Tree<String> getStructuralRNATree() {
        if (this.structuralTree == null)
            buildStructural();
        return this.structuralTree;
    }

    private void buildStructural() {

        int[] c = new int[this.sequenceLength];
        ArrayList<Integer>[] p = new ArrayList[this.sequenceLength];
        int[] m = new int[this.sequenceLength];

        //initialize the pointers array
        initp(p);

        // initialize counting and meets array
        initmc(m,c,p);

        // init indexes for later recursion call
        int l = 0; // left index
        int r = this.sequenceLength - 1; // right index

        // move l to the start of the structure
        while (c[l] == 0)
            l++;

        // move r to the tail of the structure
        while (c[r] == 0)
            r--;
        r++; // last closing loop has 0 count, but belongs to the loop, so it's not part of
        // the tail

        // the largest pseudoloop is now identified by the interval [l,r]
        assert c[l] >= 1 && c[r] == 0 : "Largest pseudoloop at [" + l + "," + r + "]\nCounting array: " + Arrays.toString(c);

        // create an empty list of zero intervals to detect concatenations
        ArrayList<Interval> zi = new ArrayList<>();

        // find zero intervals in the outermost pseudoloop, if any
        detectZeroIntervals(c, zi, l, r);

        // create an empty list of meets indexes to detect meetings
        ArrayList<Integer> meetIndexesList = new ArrayList<>();

        //find meets in the new pseudoloop, if any
        getMeetsInInterval(meetIndexesList,l,r,c,m);

        // create the root node of the structural RNA tree
        Tree<String> t = new Tree<>();

        // start the recursive construction of the structural RNA Tree on the node ct
        recBuildStructural(t, meetIndexesList, zi, c, p, m, l, r);

        // assign to the root of this tree
        this.structuralTree = t;
    }

    private void recBuildStructural(Tree<String> ct, ArrayList<Integer> meetsInInterval, ArrayList<Interval> zi, int[] c, ArrayList<Integer>[] p, int[] m, int l, int r) {

        assert c[l] >= 1 && c[r] == 0 : "Pseudoloop bounds error while parsing at [" + l + "," + r
                + "]\nCounting array: " + Arrays.toString(c);

        if((!meetsInInterval.isEmpty() && Objects.equals(p[r].get(0), meetsInInterval.get(meetsInInterval.size() - 1))) || !zi.isEmpty()) {
            int rl = 0;
            int lr = 0;
            int[] meetConcatC = c.clone();
            ArrayList<Integer>[] leftP = p;
            ArrayList<Integer>[] rightP = p;
            if(!meetsInInterval.isEmpty() && Objects.equals(p[r].get(0), meetsInInterval.get(meetsInInterval.size() - 1))){
                //meet case
                ct.setValue(Operators.MEETING_LABEL);
                //filter parents on indexes l, m and r
                leftP = p.clone();
                leftP[l] = this.filterPartnerList(l,l,meetsInInterval.get(meetsInInterval.size() - 1),p);
                leftP[meetsInInterval.get(meetsInInterval.size() - 1)] = this.filterPartnerList(meetsInInterval.get(meetsInInterval.size() - 1),l,meetsInInterval.get(meetsInInterval.size() - 1),p);

                rightP = p.clone();
                rightP[r] = this.filterPartnerList(r,meetsInInterval.get(meetsInInterval.size() - 1),r,p);
                rightP[meetsInInterval.get(meetsInInterval.size() - 1)] = this.filterPartnerList(meetsInInterval.get(meetsInInterval.size() - 1),meetsInInterval.get(meetsInInterval.size() - 1),r,p);

                //find boundaries of the right pseudoloop and of the left pseudoloop
                lr = meetsInInterval.get(meetsInInterval.size() - 1);
                rl = meetsInInterval.get(meetsInInterval.size() - 1);

                //When two pseudoloops are split by a meet, the c array changes
                meetConcatC[meetsInInterval.get(meetsInInterval.size() - 1)] = 0;

                meetsInInterval.remove(meetsInInterval.size() - 1);
            } else if(!zi.isEmpty()){
                //concat case
                ct.setValue(Operators.CONCATENATION_LABEL);

                // get rightmost zero interval
                Interval rmzi = zi.get(zi.size() - 1);
                zi.remove(zi.size() - 1);

                // find boundaries of the right pseudoloop and of the left pseudoloop
                lr = rmzi.i - 1;
                rl = rmzi.j;
            }
            // the new right pseudoloop to consider has bounds [rl,rr]
            assert c[rl] >= 1 && c[r] == 0 : "Determined wrong pseudoloop at [" + rl + "," + r
                    + "]\nCounting array: " + Arrays.toString(c);

            // create the node for building the left part
            Tree<String> left = new Tree<>();

            // create node for building the right part
            Tree<String> right = new Tree<>();

            // update tree
            ArrayList<Tree<String>> meetConcChilds = new ArrayList<>();
            meetConcChilds.add(left);
            meetConcChilds.add(right);
            ct.replaceChildrenListBy(meetConcChilds);

            // create an empty list of zero intervals to detect concatenations in the right
            // pseudoloop
            ArrayList<Interval> zir = new ArrayList<>();

            // find zero intervals in the right pseudoloop, if any
            detectZeroIntervals(c, zir, rl, r);

            // create an empty list of meets indexes to detect meetings
            ArrayList<Integer> meetIndexesListRight = new ArrayList<>();
            ArrayList<Integer> meetIndexesListLeft = new ArrayList<>();

            //find meets in the new pseudoloop, if any
            getMeetsInInterval(meetIndexesListRight, rl, r, c, m);
            getMeetsInInterval(meetIndexesListLeft,l,lr,c,m);

            // recursive construction of the algebraic RNA subTree on the right
            recBuildStructural(right, meetIndexesListRight, zir, c, rightP, m, rl, r);
            // recursive construction of the algebraic RNA subTree on the left
            recBuildStructural(left, meetIndexesListLeft, zi, meetConcatC, leftP, m, l, lr);
        } else {
            if(p[r].get(0) > l){
                //cross case
                // determine number of crossings and set label
                int numberOfCrossings = determineNumberOfCrossings(p, p[r].get(0));
                ct.setValue("(" + Operators.CROSSING_LABEL + "," + numberOfCrossings + ")");

                // left end of the rightmost crossing hairpin
                int lpp = p[r].get(0);
                // index for the right end of the new pseudoloop
                int rp = r;

                // decrease counting array according to the elimination of this hairpin
                for (int i = lpp; i < r; i++)
                    c[i]--;

                // determine the ending of the last loop on the right
                while (c[rp] == 0)
                    rp--;
                rp++; // last closing loop has 0 count, but belongs to the loop

                // the new pseudoloop to consider has bounds [l,rp]
                assert c[l] >= 1 && c[rp] == 0 : "Determined wrong pseudoloop at [" + l + "," + rp
                        + "]\nCounting array: " + Arrays.toString(c);

                // create the empty node for building the rest of the tree on the left
                Tree<String> rest = new Tree<>();

                // create hairpin subtree
                Tree<String> h = new Tree<>();
                h.setValue(Operators.HAIRPIN_LABEL + "(" + p[r].get(0) + "," + r + ")");

                // update tree
                ArrayList<Tree<String>> crossChilds = new ArrayList<>();
                crossChilds.add(rest);
                crossChilds.add(h);
                ct.replaceChildrenListBy(crossChilds);

                // create an empty list of zero intervals to detect concatenations
                ArrayList<Interval> zip = new ArrayList<>();

                // find zero intervals in the new pseudoloop, if any
                detectZeroIntervals(c, zip, l, rp);    //errore, dopo un meet r non è uguale a 0, archi multipli sballano tutto

                // create an empty list of meets indexes to detect meetings
                ArrayList<Integer> meetIndexesList = new ArrayList<>();

                //find meets in the new pseudoloop, if any
                getMeetsInInterval(meetIndexesList, l, rp, c, m);

                //remove partner indexes from both hairpin ending's p arrays
                p[p[r].get(0)].remove(Integer.valueOf(r));
                p[r].remove(p[r].get(0));

                // recursive construction of the algebraic RNA subTree on the node rest
                recBuildStructural(rest, meetIndexesList, zip, c, p, m, l, rp);

            } else {
                if(p[r].size() > 1 && p[p[r].get(0)].size() == 1 && p[r].get(0) == l){
                    //ending case
                    ct.setValue(Operators.ENDING_LABEL);
                    // decrease counting array according to the elimination of this hairpin
                    for (int i = l; i < r; i++)
                        c[i]--;

                    // init indexes for later recursion call
                    int lp = l + 1;

                    // determine the starting of the next loop on the left
                    while (c[lp] == 0 && lp < r)
                        lp++;
                    if (lp == r) {
                        // no subloops of this ending, there will be no more complex subtrees
                        // revert to just a single hairpin
                        ct.setValue(Operators.HAIRPIN_LABEL + "(" + p[r].get(0) + "," + r + ")");

                        // end recursion
                        return;
                    }

                    // the new pseudoloop to consider has bounds [lp,rp]
                    assert c[lp] >= 1 && c[r] == 0 : "Determined wrong pseudoloop at [" + lp + "," + r + "]\nCounting array: "
                            + Arrays.toString(c);

                    // create the empty node for building the rest of the tree on the left
                    Tree<String> rest = new Tree<>();

                    // create hairpin subtree
                    Tree<String> h = new Tree<>();
                    h.setValue(Operators.HAIRPIN_LABEL + "(" + p[r].get(0) + "," + r + ")");

                    // update tree
                    ArrayList<Tree<String>> endChild = new ArrayList<>();
                    endChild.add(rest);
                    endChild.add(h);
                    ct.replaceChildrenListBy(endChild);

                    // create an empty list of zero intervals to detect concatenations
                    ArrayList<Interval> zip = new ArrayList<>();

                    // find zero intervals in the new pseudoloop, if any
                    detectZeroIntervals(c, zip, lp, r);

                    // create an empty list of meets indexes to detect meetings
                    ArrayList<Integer> meetIndexesList = new ArrayList<>();

                    //find meets in the new pseudoloop, if any
                    getMeetsInInterval(meetIndexesList, lp, r, c, m);

                    //remove partner indexes from both hairpin ending's p arrays
                    p[p[r].get(0)].remove(Integer.valueOf(r));
                    p[r].remove(p[r].get(0));

                    // recursive construction of the algebraic RNA subTree on the node rest
                    recBuildStructural(rest, meetIndexesList, zip, c, p, m, lp, r);

                } else if(p[r].size() == 1 && p[p[r].get(0)].size() > 1 && p[r].get(0) == l){
                    //starting case
                    ct.setValue(Operators.STARTING_LABEL);

                    // decrease counting array according to the elimination of this hairpin
                    for (int i = l; i < r; i++)
                        c[i]--;

                    // init indexes for later recursion call
                    int rp = r - 1;

                    // determine the ending of the last loop on the right
                    while (c[rp] == 0)
                        rp--;
                    rp++; // last closing loop has 0 count, but belongs to the loop

                    // the new pseudoloop to consider has bounds [lp,rp]
                    assert c[l] >= 1 && c[rp] == 0 : "Determined wrong pseudoloop at [" + l + "," + rp + "]\nCounting array: "
                            + Arrays.toString(c);

                    // create the empty node for building the rest of the tree on the left
                    Tree<String> rest = new Tree<>();

                    // create hairpin subtree
                    Tree<String> h = new Tree<>();
                    h.setValue(Operators.HAIRPIN_LABEL + "(" + p[r].get(0) + "," + r + ")");

                    // update tree
                    ArrayList<Tree<String>> startChids = new ArrayList<>();
                    startChids.add(rest);
                    startChids.add(h);
                    ct.replaceChildrenListBy(startChids);

                    // create an empty list of zero intervals to detect concatenations
                    ArrayList<Interval> zip = new ArrayList<>();

                    // find zero intervals in the new pseudoloop, if any
                    detectZeroIntervals(c, zip, l, rp);

                    // create an empty list of meets indexes to detect meetings
                    ArrayList<Integer> meetIndexesList = new ArrayList<>();

                    //find meets in the new pseudoloop, if any
                    getMeetsInInterval(meetIndexesList, l, rp, c, m);

                    //remove partner indexes from both hairpin ending's p arrays
                    p[p[r].get(0)].remove(Integer.valueOf(r));
                    p[r].remove(p[r].get(0));

                    // recursive construction of the structural RNA subTree on the node rest
                    recBuildStructural(rest, meetIndexesList, zip, c, p, m, l, rp);

                } else if(p[r].size() > 1 && p[p[r].get(0)].size() > 1 && p[r].get(0) == l){
                    //diamond case
                    ct.setValue(Operators.DIAMOND_LABEL);

                    // decrease counting array according to the elimination of this hairpin
                    for (int i = l; i < r; i++)
                        c[i]--;

                    // indexes for later recursion call are l and r

                    // Starting and ending of the loop are l and r

                    // the new pseudoloop to consider has bounds [lp,rp]
                    assert c[l] >= 1 && c[r] == 0 : "Determined wrong pseudoloop at [" + l + "," + r + "]\nCounting array: "
                            + Arrays.toString(c);

                    // create the empty node for building the rest of the tree on the left
                    Tree<String> rest = new Tree<>();

                    // create hairpin subtree
                    Tree<String> h = new Tree<>();
                    h.setValue(Operators.HAIRPIN_LABEL + "(" + p[r].get(0) + "," + r + ")");

                    // update tree
                    ArrayList<Tree<String>> diamondChilds = new ArrayList<>();
                    diamondChilds.add(rest);
                    diamondChilds.add(h);
                    ct.replaceChildrenListBy(diamondChilds);

                    // create an empty list of zero intervals to detect concatenations
                    ArrayList<Interval> zip = new ArrayList<>();

                    // find zero intervals in the new pseudoloop, if any
                    detectZeroIntervals(c, zip, l, r);

                    // create an empty list of meets indexes to detect meetings
                    ArrayList<Integer> meetIndexesList = new ArrayList<>();

                    //find meets in the new pseudoloop, if any
                    getMeetsInInterval(meetIndexesList, l, r, c, m);

                    //remove partner indexes from both hairpin ending's p arrays
                    p[p[r].get(0)].remove(Integer.valueOf(r));
                    p[r].remove(p[r].get(0));

                    // recursive construction of the structural RNA subTree on the node rest
                    recBuildStructural(rest, meetIndexesList, zip, c, p, m, l, r);

                } else {
                    //nest case
                    ct.setValue(Operators.NESTING_LABEL);

                    // decrease counting array according to the elimination of this hairpin
                    for (int i = l; i < r; i++)
                        c[i]--;

                    // init indexes for later recursion call
                    int lp = l + 1;
                    int rp = r;

                    // determine the starting of the next loop on the left and increment k
                    while (c[lp] == 0 && lp < rp)
                        lp++;
                    if (lp == rp) {
                        // no subloops of this nesting, there will be no more complex subtrees

                        // revert to just a single hairpin
                        ct.setValue(Operators.HAIRPIN_LABEL + "(" + p[r].get(0) + "," + r + ")");

                        // end recursion
                        return;
                    }

                    // determine the ending of the last loop on the right
                    while (c[rp] == 0)
                        rp--;
                    rp++; // last closing loop has 0 count, but belongs to the loop

                    // the new pseudoloop to consider has bounds [lp,rp]
                    assert c[lp] >= 1 && c[rp] == 0 : "Determined wrong pseudoloop at [" + lp + "," + rp + "]\nCounting array: "
                            + Arrays.toString(c);

                    // create the empty node for building the rest of the tree on the left
                    Tree<String> rest = new Tree<>();

                    // create hairpin subtree
                    Tree<String> h = new Tree<>();
                    h.setValue(Operators.HAIRPIN_LABEL + "(" + p[r].get(0) + "," + r + ")");

                    // update tree
                    ArrayList<Tree<String>> nestChilds = new ArrayList<>();
                    nestChilds.add(rest);
                    nestChilds.add(h);
                    ct.replaceChildrenListBy(nestChilds);

                    // create an empty list of zero intervals to detect concatenations
                    ArrayList<Interval> zip = new ArrayList<>();

                    // find zero intervals in the new pseudoloop, if any
                    detectZeroIntervals(c, zip, lp, rp);

                    // create an empty list of meets indexes to detect meetings
                    ArrayList<Integer> meetIndexesList = new ArrayList<>();

                    //find meets in the new pseudoloop, if any
                    getMeetsInInterval(meetIndexesList, lp, rp, c, m);

                    //remove partner indexes from both hairpin ending's p arrays
                    p[p[r].get(0)].remove(Integer.valueOf(r));
                    p[r].remove(p[r].get(0));

                    // recursive construction of the structural RNA subTree on the node rest
                    recBuildStructural(rest, meetIndexesList, zip, c, p, m, lp, rp);

                }
            }
        }
    }

    private ArrayList<Integer> filterPartnerList(int index, int l, int r, ArrayList<Integer>[] p) {
        ArrayList<Integer> filteredPartnerList = new ArrayList<>();
        p[index].forEach(partnerIndex -> {
            if(partnerIndex >= l && partnerIndex <= r)
                filteredPartnerList.add(partnerIndex);
        });
        return filteredPartnerList;
    }

    /**
     * Returns true if a meet is found in the array between l and r indexes
     * @param l left index
     * @param r right index
     */
    private void getMeetsInInterval(ArrayList<Integer> meetList, int l, int r, int[] c, int[] m){
        for(int i = l + 1; i < r; i++){
            if(c[i] < m[i]){
                meetList.add(i);
            }
        }
    }

    private int determineNumberOfCrossings(ArrayList<Integer>[] p, int bondStart) {
        int n = 0;
        for(int i=0; i<p.length; i++){
            if(p[i] != null) {
                for (Integer currentPartner : p[i]) {
                    if (i < bondStart && bondStart < currentPartner)
                        n++;
                }
            }
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
        assert c[l] >= 1 && c[r] == 0 : "Pseudoloop bounds error while detecting zero intervals at [" + l + "," + r
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

    private void initmc(int[] m, int[] c, ArrayList<Integer>[] p) {
        int count = 0;
        int currentIndexStartingLoops;
        int currentIndexStoppingLoops;
        for (int i = 0; i < this.sequenceLength; i++) {
            if (!(p[i] == null)) {
                currentIndexStartingLoops = getStartingLoopsNumber(i,p);
                currentIndexStoppingLoops = getStoppingLoopsNumber(i,p);
                count = count + (currentIndexStartingLoops - currentIndexStoppingLoops);
                if (currentIndexStartingLoops > 0 && currentIndexStoppingLoops > 0) {
                    // a loop stops and another starts in the current index, so this is a meet
                    m[i] = p[i].size();
                }
            }
            c[i] = count;
        }
        assert count == 0 : "Value of count after initialization of counting array: " + count + "\nCounting array: "
                + Arrays.toString(c);
    }

    private void initp(ArrayList<Integer>[] p){
        for(Pair<Integer> currentBond : this.tertiaryStructure.getBondList()) {
            if (!currentBond.getFirst().equals(currentBond.getSecond())) {
                if (p[currentBond.getFirst()] == null) {
                    p[currentBond.getFirst()] = new ArrayList<>();
                }
                p[currentBond.getFirst()].add(currentBond.getSecond());
                if (p[currentBond.getSecond()] == null) {
                    p[currentBond.getSecond()] = new ArrayList<>();
                }
                p[currentBond.getSecond()].add(currentBond.getFirst());
            }
        }
    }

    /*
     * Tells if at position i there is the starting of an hairpin loop of the
     * secondary structure represented by the original arc annotated sequence.
     */
    private int getStartingLoopsNumber(int i, ArrayList<Integer>[] p){
        int startingLoopsNumber = 0;
        if(!(p[i] == null))
            for(Integer partner : p[i])
                if(partner > i)
                    startingLoopsNumber ++;
        return startingLoopsNumber;
    }

    /*
     * Tells if at position i there is the ending of an hairpin loop of the
     * secondary structure represented by the original arc annotated sequence.
     */
    private int getStoppingLoopsNumber(int i,  ArrayList<Integer>[] p) {
        int stoppingLoopsNumber = 0;
        if(!(p[i] == null))
            for(Integer partner : p[i])
                if(partner < i)
                    stoppingLoopsNumber ++;
        return stoppingLoopsNumber;
    }

    /**
     * FOR TESTS PURPOSES
     * Replace current sequence length with a new one.
     * @param sequenceLength new sequence length
     */
    public void setSequenceLength(int sequenceLength){
        this.sequenceLength = sequenceLength;
    }

    /*
     * Service class for holding zero intervals.
     */
    protected class Interval {
        private final int i;
        private final int j;

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

