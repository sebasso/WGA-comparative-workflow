#!/usr/bin/env python
# -*- coding: utf-8 -*-

from ete3 import Tree
import sys
import json
import os

"""
NB: robinson_fould doesnt compare leafs that ARENT in the other tree, skips them.
so only similar leafs are part of the stats.
"compare topologies of different size and content. When two trees contain a different set of labels, only shared leaves will be used."
"""
# IF trees contain different leaves, the tree with extra leafes are PRUNED/deleted.
def robinson_fould(outputdir, input_trees):

    ete_trees = []
    treenames = []
    for tree in input_trees:
        treenames.append(os.path.basename(tree).split(".")[0]+"pruned.tree")
        ete_trees.append(Tree(tree))

    rf, max_rf, common_leaves, parts_t1, parts_t2,u1,u2 = ete_trees[0].robinson_foulds(ete_trees[1],unrooted_trees=True)

    stats = {}
    stats["Robinson foulds"] = rf
    stats["Splits in tree1 but not tree2"] = list(parts_t1 - parts_t2)
    stats["Splits in tree2 but not tree1"] = list(parts_t2 - parts_t1)
    if len(common_leaves) > 0:
        stats["Common leaves"] = list(common_leaves)
        stats["num_common_leaves"] = len(common_leaves)

    for i in range(0,len(input_trees)):
        cached_content = ete_trees[i].get_cached_content()
        num_leafes = len(cached_content[ete_trees[i]])
        #total_nodes = len(cached_content)
        stats[os.path.basename(input_trees[i][:-5])+" stats"] = {"num genomes": num_leafes}

    #TODO: check cached_content vs len(treenode)
    #num_genomes = common_leaves always when comparing
    r = json.dumps(stats, indent=4, encoding=" utf-8", sort_keys=True)


    with open(outputdir+"/tree_stats.json","w") as f: #warning will write this relative to exection path - sys.executables
        f.write(r)

    if len(ete_trees[0]) > len(common_leaves) or len(ete_trees[1]) > len(common_leaves):
        print "num leaves before pruning:", len(ete_trees[0]), len(ete_trees[1])

        for leave in ete_trees[0].iter_leaves():
            if not leave.name in common_leaves:
                leave.delete()

        for leave in ete_trees[1].iter_leaves():
            if not leave.name in common_leaves:
                leave.delete()

        print "num leaves after pruning:", len(ete_trees[0]), len(ete_trees[1])
        for trees in ete_trees:
            for i in xrange(0,len(ete_trees)):
                with open(outputdir+"/"+treenames[i],"w") as f:
                    f.write(ete_trees[i].write())



if __name__ == '__main__':
    print "tree_comparator args:"
    print str(sys.argv)+"\n"
    files = sys.argv[1:]
    if len(files) <= 2:
        sys.stderr.write("tree_comparator requires a minimum of 2 trees")
        exit(1)

    robinson_fould(files[0],files[1:])
