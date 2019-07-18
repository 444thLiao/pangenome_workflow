# -*- coding: utf-8 -*-
import sys
from os.path import dirname, join

sys.path.insert(0, dirname(dirname(__file__)))
import warnings
from luigi.cmdline_parser import CmdlineParser
from luigi.task import flatten
from toolkit.utils import run_cmd
from io import StringIO
from collections import Counter

warnings.filterwarnings("ignore")


def get_graph(task):
    children = flatten(task.requires())
    count_c = Counter([_.__class__.__name__
                       for _ in children])
    graph = {"nodes": [task],
             "edges": [(task, c, count_c[c.__class__.__name__])
                       for c in children]}

    for index, child in enumerate(children):
        result_g = get_graph(child)
        graph["nodes"] += result_g["nodes"]
        graph["edges"] += result_g["edges"]
    return graph


def construct_dot_output(graph):
    # graph nodes and edges is a list of luigi.Task
    stream = StringIO()

    stream.write('digraph G {\n')
    # stream.write('size="8,11";\n')
    stream.write('splines=true;\n')
    stream.write('fontsize="30";\n')
    stream.write('ranksep = 0.3;\n')
    stream.write('node[shape=box,fontsize="20"];\n')
    stream.write('graph[clusterrank="local"];\n')
    edges = [(p1.__class__.__name__,
              p2.__class__.__name__,
              num)
             for (p1, p2, num) in graph["edges"]]

    for p1, p2, num in set(edges):
        if num > 1:
            stream.write(f"{p1} -> {p2}[style=bold, color=red, arrowtype=normal];\n")
        else:
            stream.write(f"{p1} -> {p2};\n")
    stream.write('}')
    return stream.getvalue()


def set_graph_attributes():
    # todo
    pass


def main():
    cmdline_args = sys.argv[1:]
    if "--tab" not in cmdline_args:
        cmdline_args += ["--tab", join(dirname(__file__),
                                       "static_tree_for.tab")]
    with CmdlineParser.global_instance(cmdline_args) as cp:
        task = cp.get_task_obj()
        graph = get_graph(task)
        dot_graph = construct_dot_output(graph)
        with open('/tmp/tmp.dot', 'w') as f1:
            f1.write(dot_graph)
        dot_graph = '/tmp/tmp.dot'
        ofile = join(task.odir, "pipelines.png")
        run_cmd(f"dot -Tpng < {dot_graph} > {ofile}", )
        # run_cmd(f"rm {dot_graph}")


if __name__ == '__main__':
    main()
