#!/usr/bin/env python
# # -*- coding: utf-8 -*-
# import sys
# sys.path.append('/home/monica/Documents/ano_1/2_semestre/projeto/pydna/pydna/')

import pytest
from dseqrecord import Dseqrecord
from golden_gate import *


def test_graph_assembly():
    seqs = [Dseqrecord("ATCG"), Dseqrecord("CGTA")]
    graph = graph_assembly(seqs)
    assert len(graph.nodes) == 2
    assert len(graph.edges) == 2

def test_graph_assembly_error():
    seqs = [Dseqrecord("ATCG"), Dseqrecord("GCTA")]
    with pytest.raises(ValueError):
        graph_assembly(seqs)

def test_find_all_paths():
    graph = nx.DiGraph()
    graph.add_edge(1, 2)
    graph.add_edge(2, 3)
    paths = find_all_paths(graph)
    assert paths == [[1, 2, 3]]

def test_find_all_paths_no_paths():
    graph = nx.DiGraph()
    graph.add_node(1, dseq=Dseqrecord("ATCG"))
    graph.add_node(2, dseq=Dseqrecord("CGTA"))
    paths = find_all_paths(graph)
    assert paths == []

def test_find_paths_seqs():
    graph = nx.DiGraph()
    graph.add_node(1, dseq=Dseqrecord("ATCG"))
    graph.add_node(2, dseq=Dseqrecord("CGTA"))
    paths = [[1, 2], [2, 1]]
    sequences = find_paths_seqs(paths, graph)
    assert sequences == {(1, 2): Dseqrecord("ATCGCGTA")}


def test_GoldenGateAssembler():
    seqs = [Dseqrecord("ATCG"), Dseqrecord("CGTA")]
    sequences = GoldenGateAssembler(seqs)
    assert sequences == {(1, 2): Dseqrecord("ATCGCGTA")}

def test_GoldenGateAssembler_circular():
    seqs = [Dseqrecord("ATCG"), Dseqrecord("CGAT", circular=True)]
    sequences = GoldenGateAssembler(seqs)
    assert sequences == {(1, 2): Dseqrecord("ATCGCGAT", circular=True)}

