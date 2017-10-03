import java.util.ArrayList;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Created by sebastiansoberg on 10/09/17.
 */
public class DFS {

    double accumulatedScore;
    int numLeaves;
    String fileEnding;
    DFS(NewickTree.Node root, String fileEnding ){
        this.accumulatedScore = 0.0;
        this.numLeaves = 0;
        this.fileEnding = fileEnding;
        this.depthFirstSearch(root);
        this.printResults();
    }

    void printResults(){
        double res = (this.accumulatedScore / (double)(this.numLeaves - 1)) * 100;
        //System.out.println("numleaves:" + this.numLeaves);
        System.out.println("Score:\t"+ this.accumulatedScore);
        System.out.println("Similarity in %:\t "+ res);
    }


     void depthFirstSearch(NewickTree.Node n){
        for(NewickTree.Node child: n.children){
            if(child.children.size() == 0) { // child is leaf, find similarity score
                child.setTraversedNodes();
                //System.out.println("Calc sim for "+ child.getName()+ " with parent: "+n.getName());
                HashSet<NewickTree.Node> internalSameDistanceFromLeafNodes = new HashSet<>();
                calc_child_sim_score(n, child, 0, internalSameDistanceFromLeafNodes);
                this.accumulatedScore += child.getSimilarity();
                this.numLeaves++;
            }else{ //child is non leaf, traverse.
                depthFirstSearch(child);
            }
        }
    }

    void calc_child_sim_score(NewickTree.Node parent, NewickTree.Node leaf, int level, HashSet<NewickTree.Node> internalSameDistanceFromLeafNodes){
        if (parent.getNumChildrenLeafsExceptCurr(leaf) < 2 ){
            if(parent.getNumChildrenLeafsExceptCurr(leaf) == 1) {
                for (NewickTree.Node leaves : parent.children) {
                    if (!leaf.equals(leaves) && leaves.children.size() == 0) { // only compare leafs
                        leaf.increaseNumCompares();
                        compare(leaves, leaf); // COMPARE nodes
                    }
                }
            }
            internalSameDistanceFromLeafNodes.addAll(parent.getInternalChildrenNodes());
            if(parent.parent != null){
                internalSameDistanceFromLeafNodes.add(parent.parent);
            }
            helper(parent, leaf, level, internalSameDistanceFromLeafNodes);
        }else{
            iterateChildrenAndCompare(parent, leaf);
        }
        //fixandPrintStats(leaf);
    }

    void helper(NewickTree.Node initialParent, NewickTree.Node leaf, int level, HashSet<NewickTree.Node> internalSameDistanceFromLeafNodes){
        boolean nyeNoder = true;
        boolean firstIteration = true;

        ArrayList<NewickTree.Node> newNodes = new ArrayList<>();
        HashSet<NewickTree.Node> loopedNodes = new HashSet<>();
        loopedNodes.addAll(internalSameDistanceFromLeafNodes);
        /*
            IF internalSameDistance Fromleafasnodes
            we dont need two methods for it

         */
        while(nyeNoder) {
            loopedNodes.addAll(internalSameDistanceFromLeafNodes);
            /*System.out.println("presize: "+internalSameDistanceFromLeafNodes.size()+" Nodes pre:");
            for(NewickTree.Node n: internalSameDistanceFromLeafNodes){
                System.out.println(n.getName());
            }*/
            int loops = 0;
            for (NewickTree.Node currentInternalNodes : internalSameDistanceFromLeafNodes) {
                //System.out.println("iterating for : "+currentInternalNodes.getName());
                ArrayList<NewickTree.Node> downNodes = recursiveDownStep(currentInternalNodes, leaf,0, firstIteration);
                if(downNodes!= null)
                    newNodes.addAll(downNodes);
                /*System.out.println("downstepnodes: ");
                for(NewickTree.Node tmp2: newNodes){
                    System.out.println(tmp2.getName());
                }*/

               // TODO: also add parent somewhere to new internal nodes.
                // now ONLY internals are added
                /*if(currentInternalNodes.parent != null)
                    newNodes.add(currentInternalNodes.parent); // or similar*/

                loops++;
            }
            internalSameDistanceFromLeafNodes.removeAll(loopedNodes);
            internalSameDistanceFromLeafNodes.addAll(newNodes);
            newNodes.clear();
            /*System.out.println("looped: "+loops);
            System.out.println(" Nodes post:");
            for(NewickTree.Node n: internalSameDistanceFromLeafNodes){
                System.out.println(n.getName());
            }
            System.out.println("leafcompares: "+leaf.get_num_compares());*/

            if (leaf.get_num_compares() >= 2) {
                //System.out.println("leafcompares: "+leaf.get_num_compares());
                return;
            }
            if(leaf.getCycle()){
                return;
            }
            if(internalSameDistanceFromLeafNodes.size() > 0){
                nyeNoder = true;
                level++;
            }else{
                nyeNoder = false;
            }
            if(firstIteration) {
               // System.out.println("resets size. "+internalSameDistanceFromLeafNodes.size());
                leaf.traversedNodes.addAll(internalSameDistanceFromLeafNodes);
                firstIteration = false;
            }
            leaf.increaseLevels();
        }
    }

     NewickTree.Node recursiveUpStep(NewickTree.Node parent, NewickTree.Node leaf, int level, boolean firstIter){
        if(firstIter) {
            if (parent.parent != null) {
                //System.out.println("recursice up comparing againtst: "+parent.parent.getName());
                calc_child_sim_score_one_step(parent.parent, leaf, (level + 1), firstIter); // check parent parents
            } else {
                return null;
            }
            return parent.parent;
        }
        return null;
    }

     ArrayList<NewickTree.Node> recursiveDownStep(NewickTree.Node internalNode, NewickTree.Node leaf, int level, boolean firstIter){
         if (leaf.traversedNodes.contains(internalNode)) {
             leaf.setCycle();
             //System.out.println("setcycle for :" + leaf.getName() + " multiple of node: " + internalNode.getName());
             return null;
         }else{
             leaf.traversedNodes.add(internalNode);
             for (NewickTree.Node child : internalNode.children) {
                 if (!leaf.equals(child) && child.children.size() == 0) {
                     //System.out.println("\t recursive down comparing against: " + child.getName());
                     leaf.increaseNumCompares();
                     compare(child, leaf);
                 }
             }
         }
        return internalNode.getInternalChildrenNodes();
    }

    void calc_child_sim_score_one_step(NewickTree.Node internalNode, NewickTree.Node leaf, int level, boolean firstIter){
        if (!firstIter) {
            if (leaf.traversedNodes.contains(internalNode)) { // avoid cycles
                leaf.setCycle();
                //System.out.println("setcycle for :" + leaf.getName() + " multiple of " + internalNode.getName());
                return;
            } else {
                leaf.traversedNodes.add(internalNode);
            }
        }
        iterateChildrenAndCompare(internalNode, leaf);
    }

    void iterateChildrenAndCompare(NewickTree.Node parent, NewickTree.Node leaf){
        for(NewickTree.Node child: parent.children){
            if(!leaf.equals(child) && child.children.size() == 0) { // only compare leafs
                leaf.increaseNumCompares();
                compare(child, leaf); // COMPARE nodes
            }
        }
    }

    void compare(NewickTree.Node neighbour, NewickTree.Node leaf){
        String pattern = "";
        if(this.fileEnding != null){
            pattern = fileEnding+"(.*)";
        }else{
            pattern = "fasta(.*)";
        }

        Pattern p = Pattern.compile( pattern ); // genome simulator appends number at end of filename
        Matcher neigh = p.matcher( neighbour.getName() );
        Matcher leaff = p.matcher( leaf.getName() );
        if ( neigh.find() && leaff.find() ) {
            String neighh = neigh.group(1);
            String leafff = leaff.group(1);
            int diff = Math.abs( Integer.parseInt(neighh) - Integer.parseInt(leafff));
            if(diff == 1){
                leaf.increaseSimilarity(0.5);
            }
        }
        leaf.comparedTo.add(neighbour);
    }

    void fixandPrintStats(NewickTree.Node n){
        n.printcomparedTo();
        n.clearTraversedNodes();
    }

}
