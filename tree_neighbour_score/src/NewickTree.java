/**
 * Created by sebastiansoberg on 05/09/17.
 */
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.stream.Collectors;

public class NewickTree {

    private static int node_uuid = 0;
    public ArrayList<Node> nodeList = new ArrayList<>();
    public Node root;

    static NewickTree readNewickFormat(String newick) {
        return new NewickTree().innerReadNewickFormat(newick);
    }

    private static String[] split(String s) {

        ArrayList<Integer> splitIndices = new ArrayList<>();

        int rightParenCount = 0;
        int leftParenCount = 0;
        for (int i = 0; i < s.length(); i++) {
            switch (s.charAt(i)) {
                case '(':
                    leftParenCount++;
                    break;
                case ')':
                    rightParenCount++;
                    break;
                case ',':
                    if (leftParenCount == rightParenCount) splitIndices.add(i);
                    break;
            }
        }

        int numSplits = splitIndices.size() + 1;
        String[] splits = new String[numSplits];

        if (numSplits == 1) {
            splits[0] = s;
        } else {

            splits[0] = s.substring(0, splitIndices.get(0));

            for (int i = 1; i < splitIndices.size(); i++) {
                splits[i] = s.substring(splitIndices.get(i - 1) + 1, splitIndices.get(i));
            }

            splits[numSplits - 1] = s.substring(splitIndices.get(splitIndices.size() - 1) + 1);
        }

        return splits;
    }

    private NewickTree innerReadNewickFormat(String newick) {

        // single branch = subtree (?)
        this.root = readSubtree(newick.substring(0, newick.length() - 1));

        return this;
    }

    private Node readSubtree(String s) {

        int leftParen = s.indexOf('(');
        int rightParen = s.lastIndexOf(')');

        if (leftParen != -1 && rightParen != -1) {

            String name = s.substring(rightParen + 1);
            String[] childrenString = split(s.substring(leftParen + 1, rightParen));

            Node node = new Node(name);
            //node.children = new ArrayList<>();
            for (String sub : childrenString) {
                Node child = readSubtree(sub);
                node.children.add(child);
                child.parent = node;
            }

            nodeList.add(node);
            return node;
        } else if (leftParen == rightParen) {

            Node node = new Node(s);
            nodeList.add(node);
            return node;

        } else throw new RuntimeException("unbalanced ()'s");
    }

    static class Node {
        final String name;
        final String weight;
        boolean realName = false;
        ArrayList<Node> children = new ArrayList<>();
        Node parent;
        public int num_compares = 0;
        public int levels = 0;
        ArrayList<Node> comparedTo = new ArrayList<>();
        HashSet<Node> traversedNodes = new HashSet<>();
        public boolean cycle = false;

        public void setCycle(){
            this.cycle = true;
        }

        public boolean getCycle(){
            return this.cycle;
        }

        public void setTraversedNodes(){
            traversedNodes = new HashSet<>();
        }

        public void clearTraversedNodes(){
            this.traversedNodes.clear();
        }

        public ArrayList<Node> getInternalChildrenNodes(){
            List<Node> tmp = children.stream().filter(c -> c.children.size() != 0 )
                    .collect(Collectors.toList());;
            return (ArrayList<Node>) tmp;
        }

        public int getNumChildrenLeafsExceptCurr(Node curr){
            List<Node> tmp = children.stream().filter(c -> c.children.size() == 0 && !curr.equals(c))
            .collect(Collectors.toList());;
            return tmp.size();
        }

        public void printcomparedTo(){
            System.out.println(getName()+ " is compared to: ");
            for (Node n: comparedTo){
                System.out.println("\t"+n.getName());
            }
            System.out.println("Num comparisons: "+this.comparedTo.size());
            System.out.println("Num levels: "+this.levels);
            System.out.println("Score: "+this.getSimilarity());
        }
        public void increaseLevels(){
            this.levels++;
        }

        public void increaseNumCompares(){
            this.num_compares++;
        }

        public int get_num_compares(){
            return this.num_compares;
        }

        public double similarity = 0;

        public void increaseSimilarity(double score){
            this.similarity = score + this.getSimilarity();
            if(this.similarity > 1){
                System.out.println("Sim score to high for: "+getName());
            }
        }

        public double getSimilarity(){
            return this.similarity;
        }

        /**
         * @param name name in "actualName:weight" format, weight defaults to zero if colon absent
         */
        Node(String name) {

            int colonIndex = name.indexOf(':');
            String actualNameText;
            if (colonIndex == -1) {
                actualNameText = name;
                weight = "0";
            } else {
                actualNameText = name.substring(0, colonIndex);
                weight = "1";
            }

            if (actualNameText.equals("")) {
                this.realName = false;
                this.name = Integer.toString(node_uuid);
                node_uuid++;
            } else {
                this.realName = true;
                this.name = actualNameText;
            }
        }

        @Override
        public int hashCode() {
            return name.hashCode();
        }

        @Override
        public boolean equals(Object o) {
            if (!(o instanceof Node)) return false;
            Node other = (Node) o;
            return this.name.equals(other.name);
        }

        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder();
            if (children != null && children.size() > 0) {
                sb.append("(");
                for (int i = 0; i < children.size() - 1; i++) {
                    sb.append(children.get(i).toString());
                    sb.append(",");
                }
                sb.append(children.get(children.size() - 1).toString());
                sb.append(")");
            }
            if (name != null) sb.append(this.getName());
            return sb.toString();
        }

        String getName() {
            if (realName)
                return name;
            else
                return "";
        }
    }

    @Override
    public String toString() {
        return root.toString() + ";";
    }

}