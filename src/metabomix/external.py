import networkx as nx
from itertools import islice

def filter_component(G, max_component_size):
    """From GNPS! drop lowest edges from cluster if it exceeds max size """
    if max_component_size == 0:
        return

    big_components_present = True

    while big_components_present == True:
        big_components_present = False
        components = nx.connected_components(G)
        for component in components:
            if len(component) > max_component_size:
                prune_component(G, component)
                big_components_present = True
        #print("After Round of Component Pruning", len(G.edges()))
        
def prune_component(G, component, cosine_delta=0.02):
    component_edges = get_edges_of_component(G, component)

    min_score = 1000
    for edge in component_edges:
        min_score = min(min_score, edge[2]["weight"])

    cosine_threshold = cosine_delta + min_score
    for edge in component_edges:
        if edge[2]["weight"] < cosine_threshold:
            #print(edge)
            G.remove_edge(edge[0], edge[1])

def get_edges_of_component(G, component):
    component_edges = {}
    for node in component:
        node_edges = G.edges((node), data=True)
        for edge in node_edges:
            if edge[0] < edge[1]:
                key = edge[0] + "-" + edge[1]
            else:
                key = edge[1] + "-" + edge[0]
            component_edges[key] = edge

    component_edges = component_edges.values()
    return component_edges

#from 
def batched(iterable, n, *, strict=False):
    """from itertools python 3.12, split iterable in batches of size n"""
    # batched('ABCDEFG', 3) → ABC DEF G
    if n < 1:
        raise ValueError('n must be at least one')
    iterator = iter(iterable)
    while batch := tuple(islice(iterator, n)):
        if strict and len(batch) != n:
            raise ValueError('batched(): incomplete batch')
        yield batch

