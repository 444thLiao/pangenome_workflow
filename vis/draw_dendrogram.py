from Bio import Phylo
import plotly.graph_objs as go

def get_tree(tree_pth,rooted=False):
    t = Phylo.read(tree_pth, 'newick')
    if rooted is False:
        return t
    elif rooted == 'midpoint':
        t.root_at_midpoint()
        return t
    elif rooted in [_.name for _ in t.get_terminals()]:
        rooted_node = [_ for _ in t.get_terminals() if _.name == rooted]
        t.root_with_outgroup(rooted_node[0])
        return t
    return t

def main(tree_pth,rooted =False):
    def get_x_positions(tree):
        """Create a mapping of each clade to its horizontal position.

        Dict of {clade: x-coord}
        """
        depths = tree.depths()
        # If there are no branch lengths, assume unit branch lengths
        if not max(depths.values()):
            depths = tree.depths(unit_branch_lengths=True)
        return depths

    def get_y_positions(tree):
        """Create a mapping of each clade to its vertical position.

        Dict of {clade: y-coord}.
        Coordinates are negative, and integers for tips.
        """
        maxheight = tree.count_terminals()
        # Rows are defined by the tips
        heights = dict((tip, maxheight - i)
                       for i, tip in enumerate(reversed(tree.get_terminals())))

        # Internal nodes: place at midpoint of children
        def calc_row(clade):
            for subclade in clade:
                if subclade not in heights:
                    calc_row(subclade)
            # Closure over heights
            heights[clade] = (heights[clade.clades[0]] +
                              heights[clade.clades[-1]]) / 2.0

        if tree.root.clades:
            calc_row(tree.root)
        return heights

    def draw_clade_lines(orientation='horizontal',
                         y_here=0, x_start=0, x_here=0, y_bot=0, y_top=0):
        """Create a line with or without a line collection object.

        Graphical formatting of the lines representing clades in the plot can be
        customized by altering this function.
        """
        if orientation == 'horizontal':
            horizontal_linecollections.append(
                [(x_start, y_here), (x_here, y_here)])
        if orientation == 'vertical':
            vertical_linecollections.append(
                [(x_here, y_bot), (x_here, y_top)])

    def draw_clade(clade, x_start, color, lw,x_posns,y_posns):
        """Recursively draw a tree, down from the given clade."""
        x_here = x_posns[clade]
        y_here = y_posns[clade]
        # Draw a horizontal line from start to here
        draw_clade_lines(orientation='horizontal',
                         y_here=y_here, x_start=x_start, x_here=x_here)
        label = str(clade)
        labels.append(label)
        if clade.clades:
            # for those point non-terminals.
            # Draw a vertical line connecting all children
            y_top = y_posns[clade.clades[0]]
            y_bot = y_posns[clade.clades[-1]]
            # Only apply widths to horizontal lines, like Archaeopteryx

            draw_clade_lines(orientation='vertical',
                             x_here=x_here, y_bot=y_bot, y_top=y_top)
            # Draw descendents
            for child in clade:
                draw_clade(child, x_here, color, lw,x_posns,y_posns)

    if type(tree_pth) == str:
        tree = get_tree(tree_pth,rooted=rooted)
    else:
        tree = tree_pth

    x_posns = get_x_positions(tree)
    y_posns = get_y_positions(tree)
    horizontal_linecollections = []
    vertical_linecollections = []
    labels = []

    draw_clade(tree.root, 0, 'k', '', x_posns, y_posns)

    fig = go.Figure()
    labels_draw_x = []
    labels_draw_y = []
    labels = [None if _i == 'Clade' else _i for _i in labels]

    for _o, _t in horizontal_linecollections:

        horizontal_draws_x = []
        horizontal_draws_y = []
        horizontal_draws_x.append(_o[0])
        horizontal_draws_y.append(_o[1] - 1)

        if labels[horizontal_linecollections.index([_o, _t])]:
            horizontal_draws_x.append(_t[0])
            horizontal_draws_y.append(_t[1] - 1)
            labels_draw_x.append(_t[0])
            labels_draw_y.append(_t[1] - 1)
        else:
            horizontal_draws_x.append(_t[0])
            horizontal_draws_y.append(_t[1] - 1)

        trace1 = go.Scatter(x=horizontal_draws_x, y=horizontal_draws_y,
                            mode='lines',
                            line=go.scatter.Line(color='#444', width=1),
                            hoverinfo='none', xaxis='x1', yaxis='y1')
        fig.add_trace(trace1)
    for _o, _t in vertical_linecollections:
        vertical_draws_x = []
        vertical_draws_y = []
        vertical_draws_x.append(_o[0])
        vertical_draws_y.append(_o[1] - 1)
        vertical_draws_x.append(_t[0])
        vertical_draws_y.append(_t[1] - 1)
        trace1 = go.Scatter(x=vertical_draws_x, y=vertical_draws_y,
                            mode='lines',
                            line=go.scatter.Line(color='#444', width=1),
                            hoverinfo='none', xaxis='x1', yaxis='y1')
        fig.add_trace(trace1)
    fig.layout.yaxis.tickvals = labels_draw_y
    fig.layout.yaxis.ticktext = [_ for _ in labels if _ is not None]
    return fig

if __name__ == '__main__':
    pass