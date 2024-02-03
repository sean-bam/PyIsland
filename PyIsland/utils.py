from pathlib import Path
import shutil


class ete3_utils(object):
    """
    functions for helping manipulate ete3 tree objects
    """

    def mean(array):
        return sum(array) / float(len(array))

    def cache_distances(tree):
        """precalculate distances of all nodes to the root"""
        node2rootdist = {tree: 0}
        for node in tree.iter_descendants("preorder"):
            node2rootdist[node] = node.dist + node2rootdist[node.up]
        return node2rootdist

    def collapse(self, tree, min_dist):
        """
        calculates the average node-tip branch distance
        if under the given threshold, places leafs in a group
        returns a dictionary

        Modified from https://www.biostars.org/p/97409/
        """
        groups = {}
        i = 0
        # cache the tip content of each node to reduce the number of times the tree is traversed
        node2tips = get_cached_content()
        root_distance = ete3_utils.cache_distances(tree)

        for node in tree.get_descendants("preorder"):
            if not node.is_leaf():
                avg_distance_to_tips = ete3_utils.mean(
                    [
                        root_distance[tip] - root_distance[node]
                        for tip in node2tips[node]
                    ]
                )

                if avg_distance_to_tips < min_dist:
                    # set a new group number
                    i += 1
                    for tip in node2tips[node]:
                        # if it isn't already in a group
                        if tip.name not in groups.keys():
                            # assign it to the group
                            groups[tip.name] = i

        # for the singletons
        for leaf in tree.iter_leaves():
            if leaf.name not in groups.keys():
                i += 1
                groups[leaf.name] = i

        return groups


def concatenate(list_of_filepaths, output):
    """
    A wrapper function to concatenate files
    """

    output = Path(output)

    with output.open("wb") as outfile:
        for file in list_of_filepaths:
            # get the file
            filepath = Path(file)
            if not filepath.is_file():
                print(f"cant find {file}")

            # copy
            shutil.copyfileobj(filepath.open("rb"), outfile)
