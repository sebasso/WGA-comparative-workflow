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
The WGA-standalone-workflow overview consists of several modules. The elements marked in green and black are files while those in blue are direc- tories. The workflowmanager.sh script is the entry-point for the project, it’s the main file of the project and controls the execution flow. Workflow- manager uses the modules inside the flows directory, extra functionality like cleanup, abel modules, program execution and definitions for workflowmanager is here. The tool executable are in their named folders ksnp, parsnp respectively. Software that was used across modules are in the common-sw folder. The tree comp folder contains the tree comparator executables from http://kaims.eti.pg.gda.pl/~dambo/treecmp/. The snp comparator.py analyzes the snp’s from the compared tools or against a reference produced by the evolve genome module, this is also the case for the tree comparator.py module but it compares the phylogenetic trees directly. The timemodule.sh is a module that wraps the workflowmanager with timing abilities. The results folder contains all the output from comparators and tools that are relevant for a run, likewise the contains detailed timing output. The docs folder contains documentation of some of the modules used in this project.
```
