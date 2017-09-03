


# Structure
```
 WGA-comparative-workflow/
├── README.md
├── common-sw
├── docs
├── evolve_genome
├── flows
├── ksnp
├── parsnp
├── results
├── snp_comparator.py
├── timemodule.sh
├── timings
├── tree_comp
├── tree_comparator.py
└── workflowmanager.sh ```

### Description
```
The \textit{WGA-standalone-workflow} overview consists of several modules. The elements marked in green and black are files while those in blue are directories. The \textbf{workflowmanager.sh} script is the entry-point for the project, it's the main file of the project and controls the execution flow. Workflowmanager uses the modules inside the \textbf{flows} directory, extra functionality like cleanup, abel modules, program execution and definitions for workflowmanager is here.
    The tool executable are in their named folders \textbf{ksnp, parsnp} respectively. Software that was used across modules are in the \textbf{common-sw} folder. The \textbf{tree\_comp} folder contains the tree comparator executables from \url{http://kaims.eti.pg.gda.pl/~dambo/treecmp/}. The \textbf{snp\_comparator.py} analyzes the snp's from the compared tools or against a reference produced by the evolve\_genome module, this is also the case for the \textbf{tree\_comparator.py} module but it compares the phylogenetic trees directly.
    The \textbf{timemodule.sh} is a module that wraps the workflowmanager with  timing abilities.
    The \textbf{results} folder contains all the output from comparators and tools that are relevant for a run, likewise the \textbf{contains detailed timing output.}
    The \textbf{docs} folder contains documentation of some of the modules used in this project.
    ```
