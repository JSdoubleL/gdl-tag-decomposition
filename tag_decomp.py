import treeswift
import argparse

def unroot(tree):
    """
    Unroots treeswift tree. Adapted from treeswift 'deroot' function.
    This one doesn't contract (A,B); to A;

    Parameters
    ----------
    tree: treeswift tree

    Returns unrooted treeswift tree
    """
    if tree.root.num_children() == 2:
        [left, right] = tree.root.child_nodes()
        if not right.is_leaf():
            right.contract()
        elif not left.is_leaf():
            left.contract()
    tree.is_rooted = False
    return tree


def get_min_root(tree, delimiter=None):
    """
    Calculates the root with the minimum score.

    Parameters
    ----------
    tree: treeswift tree
    delimiter: delimiter separating species name from rest of leaf label

    Returns vertex corresponding to best edge to root tree on
    """
    # root tree if not rooted
    if tree.root.num_children() != 2:
        tree.reroot(tree.root)
    tree.resolve_polytomies()

    # Get down scores pass
    for node in tree.traverse_postorder():
        if node.is_leaf():
            node.down = set([node.get_label().split(delimiter)[0]])
            node.d_score = 0
        else:
            if node.num_children() != 2:
                raise Exception("Vertex has more than 2 children")

            [left, right] = node.child_nodes()
            node.down = left.down.union(right.down)
            node.d_score = left.d_score + right.d_score

            if not len(left.down.intersection(right.down)) == 0:
                if node.down == left.down or node.down == right.down:
                    if left.down == right.down:
                        node.d_score += 1
                    else:
                        node.d_score += 2
                else:
                    node.d_score += 3

    min_score, best_root = float("inf"), None

    # Get scores above edge pass
     # Find root node, also label rest of nodes 'skip = False'
    for node in tree.traverse_preorder():
        if node.is_root():
            root = node
            root.skip = True
        else:
            node.skip = False

    # Get 'up' set for children of root
    [left, right] = root.child_nodes()
    left.up = right.down
    left.u_score = right.d_score
    right.up = left.down
    right.u_score = left.d_score
    left.skip = True
    right.skip = True

    min_score = right.d_score + left.d_score
    best_root = left 

    for node in tree.traverse_preorder():
        if not node.skip:
            
            parent = node.get_parent()
            if parent.child_nodes()[0] != node:
                other = parent.child_nodes()[0]
            else: 
                other = parent.child_nodes()[1]
            # calculate up score
            node.u_score = parent.u_score + other.d_score
            
            node.up = parent.up.union(other.down)
            # check for duplications
            if not len(parent.up.intersection(other.down)) == 0:
                if node.up == parent.up or node.up == other.down:
                    if parent.up == other.down:
                        node.u_score += 1
                    else:
                        node.u_score += 2
                else:
                    node.u_score += 3 
            total_score = node.u_score + node.d_score
            # check for possible min
            if total_score < min_score:
                min_score = total_score
                best_root = node
    return best_root


def tag(tree, delimiter=None):
    """
    Tags tree according to its current rooting.

    Parameters
    ----------
    tree: treeswift tree
    delimiter: delimiter separating species name from rest of leaf label
    """
    tree.suppress_unifurcations()
    tree.resolve_polytomies()
    for node in tree.traverse_postorder():
        if node.is_leaf():
            node.s = set([node.get_label().split(delimiter)[0]])
        else:
            [left, right] = node.child_nodes()

            node.s = left.s.union(right.s)

            if len(left.s.intersection(right.s)) == 0:
                node.tag = 'S'
            else: 
                node.tag = 'D'


def decompose(tree, max_only=False, no_subsets=False):
    """
    Decomposes a tagged tree, by separating clades at duplication vetices

    NOTE: must be run after 'tag()'

    Parameters
    ----------
    tree: tagged treeswift tree
    max_only: return only the single large
    no_subsets: filters out duplicate clades that have leafsets which are subsets of their counterpart

    Returns result of the decomposition as a list of trees
    """
    out = []
    root = tree.root
    for node in tree.traverse_postorder(leaves=False):
        if node.tag == 'D':
            # trim off smallest subtree (i.e. subtree with least species)
            [left, right] = node.child_nodes()
            delete = left if len(left.s) < len(right.s) else right
            if not max_only and not (right.s.issubset(left.s) and no_subsets):
                out.append(tree.extract_subtree(delete))
            node.remove_child(delete)
    tree.suppress_unifurcations() # all the duplication nodes will be unifurcations
    out.append(tree)
    return out


def trim(tree):
    """
    Trims duplicate leaves under ever duplication vertex.

    NOTE: This may be buggy. I haven't tested it enough.

    Parameters:
    tree: tagged treeswift tree

    Returns a single tree with removed leaves
    """
    
    root = tree.root
    for node in tree.traverse_postorder(leaves=False):
        if node.tag == 'D':
            # trim only duplicate leaves from smallest subsection
            [left, right] = node.child_nodes()
            to_delete = left.s.intersection(right.s)
            smallest = left if len(left.s) < len(right.s) else right
            for v in smallest.traverse_postorder():
                if v.s.issubset(to_delete) or (v.is_leaf() and v.label is None):
                    parent = v.get_parent()
                    parent.remove_child(v)
    return tree


def trivial(newick_str):
    """
    Determines if a newick string represents a trivial tree (tree containing no quartets).

    Parameters
    ----------
    newick_str: newick string

    Returns True if tree contains less than two '('
    """
    count = 0
    for c in newick_str:
        if c == '(':
            count += 1
        if count > 1:
            return False
    return True


def main(args):

    if args.output is None:
        split = args.input.rsplit('.', 1)
        output = split[0] + '-decomp.' + split[1]
    else:
        output = args.output

    with open(args.input, 'r') as fi, open(output, 'w') as fo:
        for line in fi:
            tree = treeswift.read_tree_newick(line)
            tree.reroot(get_min_root(tree, args.delimiter))
            tag(tree, args.delimiter)
            if args.trim:
                out = [trim(tree)]
            else:
                out = decompose(tree, args.max_only, args.no_subsets)
            for t in out:
                unroot(t)
                t.suppress_unifurcations()
                nwck = t.newick()
                if not trivial(nwck):
                    fo.write(nwck + '\n')
            

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", type=str,
                        help="Input tree list file", required=True)
    parser.add_argument("-o", "--output", type=str,
                        help="Output tree list file")
    parser.add_argument('-d', "--delimiter", type=str, 
                        help="Delimiter separating species name from rest of leaf label")
    parser.add_argument('-m', '--max_only', action='store_true',
                        help="Only output maximum tree")
    parser.add_argument('-s', '--no_subsets', action='store_true',
                        help="Do not include sections of tree that are subsets")
    parser.add_argument('-t', '--trim', action='store_true',
                        help="Trim duplicate leaves--otherwise decompose")

    main(parser.parse_args())
