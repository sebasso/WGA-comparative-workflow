import java.io.*;
import java.nio.file.Files;

public class Main {

    public static void main(String[] args) throws IOException{
        //System.out.println(System.getProperty("user.dir"));
        // resolves into /Users/sebastiansoberg/IdeaProjects/testnewickparser
        /*File f = new File(System.getProperty("user.dir")+"/src/testclass.fasta");
        String contents = new String(Files.readAllBytes(f.toPath()));
        System.out.println("Contents (Java 7) : " + contents);*/

        /*NewickTree newroot = new NewickTree();
        newroot = test.readNewickFormat(contents);
        DFS search2 = new DFS(newroot.root);*/

        /*
        NewickTree test = new NewickTree();
        test = test.readNewickFormat("(EEE_Florida91-4697.fasta12:0.00055,(EEE_Florida91-4697.fasta11:0.00055,(EEE_Florida91-4697.fasta10:0.00055,(EEE_Florida91-4697.fasta9:0.00055,(EEE_Florida91-4697.fasta8:0.00055,(EEE_Florida91-4697.fasta7:0.00055,(EEE_Florida91-4697.fasta6:0.00055,(EEE_Florida91-4697.fasta5:0.00055,(EEE_Florida91-4697.fasta4:0.00055,(EEE_Florida91-4697.fasta3:0.00055,(EEE_Florida91-4697.fasta2:0.00055,(EEE_Florida91-4697.fasta0:0.00055,EEE_Florida91-4697.fasta1:0.00055)0.951:0.00055)0.952:0.00055)0.996:0.00078)0.964:0.00055)0.972:0.00055)0.946:0.00055)1:0.00103)0.919:0.00055)0.991:0.00061)0.824:0.00055)0.999:0.00095,(EEE_Florida91-4697.fasta13:0.00055,(EEE_Florida91-4697.fasta14:0.00055,(EEE_Florida91-4697.fasta15:0.00055,(EEE_Florida91-4697.fasta16:0.00055,(EEE_Florida91-4697.fasta17:0.00055,(EEE_Florida91-4697.fasta18:0.00055,EEE_Florida91-4697.fasta19:0.00055)0.972:0.00055)0.911:0.00055)1:0.00129)0.908:0.00055)0.97:0.00055)0.907:0.00055);");
        DFS search = new DFS(test.root);
        */

        File f = new File(args[0]);
        String fileEnding = null;
        if (args.length > 1)
            fileEnding = args[1];
        String fasttree = new String(Files.readAllBytes(f.toPath()));
        //System.out.println("Contents fasttree: " + fasttree);
        NewickTree test = new NewickTree();
        test = test.readNewickFormat(fasttree);
        DFS search = new DFS(test.root, fileEnding);

    }

}
