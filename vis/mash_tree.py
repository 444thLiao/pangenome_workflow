import pandas as pd
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as ssd


def getNewick(node, newick, parentdist, leaf_names):
    if node.is_leaf():
        return "%s:%.2f%s" % (leaf_names[node.id], parentdist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = "):%.2f%s" % (parentdist - node.dist, newick)
        else:
            newick = ");"
        newick = getNewick(node.get_left(), newick, node.dist, leaf_names)
        newick = getNewick(node.get_right(), ",%s" % (newick), node.dist, leaf_names)
        newick = "(%s" % (newick)
        return newick

def dist2newick(dist_path):
    dis = pd.read_csv(dist_path, sep=',', index_col=0)
    schlink = sch.linkage(ssd.squareform(dis))
    tree = sch.to_tree(schlink, False)

    newick_text = getNewick(tree, "", tree.dist, dis.index)
    with open(dist_path.replace('.dist', '.newick'), 'w') as f1:
        f1.write(newick_text)
