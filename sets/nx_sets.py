"""
sets_nx.py
----------
Rewrite of the inversion-set finder using NetworkX.

The original code manually walked up and downstream through SWORD reach
topology using while-loops and manual bookkeeping.  networkx lets us
do the same after building the SWORD data into a graph. Once turned into
an nx object, these topological operations are super efficient.

2/24/2026 -Ted L

"""

import json
import os

import networkx as nx


def is_river_reach(reach_id: int) -> bool:
    return str(reach_id)[-1] == "1"


def has_consistent_orbits(G: nx.DiGraph, node_ids: list) -> bool:
    """Return True if every node in the chain/graph shares the same orbits."""
    orbit_sets = [G.nodes[n]["swot_orbits"] for n in node_ids]
    return len(set(orbit_sets)) == 1


def has_area_jumps(G: nx.DiGraph, node_ids: list, cutoff_pct: float) -> bool:
    """
    Return True if any consecutive pair of nodes in the (ordered) chain has a
    fractional accumulation-area change that exceeds *cutoff_pct*.
    """
    for a, b in zip(node_ids, node_ids[1:]):
        facc_a = G.nodes[a]["facc"]
        facc_b = G.nodes[b]["facc"]
        if facc_a == 0:
            continue
        pct_diff = abs(facc_b - facc_a) / facc_a * 100
        if pct_diff > cutoff_pct:
            return True
    return False


def calc_overlap_frac(set_a: frozenset, set_b: frozenset) -> float:
    n_overlap = len(set_a & set_b)
    mean_size = (len(set_a) + len(set_b)) / 2
    return n_overlap / mean_size if mean_size else 0.0


class Sets:
    """Divide a list of SWORD reaches into inversion sets."""

    def __init__(self, params: dict, reaches: list, sword_dataset):
        self.params = params
        self.reaches = reaches 
        self.reach_ids = {int(r["reach_id"]) for r in reaches}
        self.sword_dataset = sword_dataset

        

    def _build_graph(self) -> nx.DiGraph:
        """
        Parse the SWORD continent file and return a directed graph where
        every node is a reach and edges point *downstream*
        (upstream_reach -> downstream_reach).

        Node attributes
        ---------------
        facc, swot_obs, swot_orbits, n_rch_up, n_rch_down
        """
        ds = self.sword_dataset
        reach_ids = ds["reaches/reach_id"][:]
        facc = ds["reaches/facc"][:]
        n_rch_up = ds["reaches/n_rch_up"][:]
        n_rch_down = ds["reaches/n_rch_down"][:]
        rch_id_dn = ds["reaches/rch_id_dn"][:]
        swot_obs = ds["reaches/swot_obs"][:]
        swot_orbits = ds["reaches/swot_orbits"][:]

        G = nx.DiGraph()

        for i, rid in enumerate(reach_ids):
            rid = int(rid)
            n_up = int(n_rch_up[i])
            n_dn = int(n_rch_down[i])
            obs = int(swot_obs[i])
            orbits = tuple(int(o) for o in swot_orbits[:obs, i])

            G.add_node(
                rid,
                facc=float(facc[i]),
                swot_obs=obs,
                swot_orbits=orbits,
                n_rch_up=n_up,
                n_rch_down=n_dn,
            )

            for dn_id in rch_id_dn[:n_dn, i]:
                dn_id = int(dn_id)
                if dn_id != 0:
                    G.add_edge(rid, dn_id)

        return G

    def _valid_river_nodes(self, G: nx.DiGraph) -> list[int]:
        """
        Return reach IDs that pass the basic structural requirements:
          - last digit == '1'  (river reach in SWORD)
          - exactly one upstream neighbour (in-degree 1)
          - exactly one downstream neighbour (out-degree 1)
          - restricted to self.reach_ids when RequireSetReachesInput is True
        """
        valid = []
        for n in G.nodes:
            if not is_river_reach(n):
                continue
            if self.params.get("RequireSetReachesInput") and n not in self.reach_ids:
                continue
            if G.in_degree(n) == 1 and G.out_degree(n) == 1:
                valid.append(n)
        return valid

    def _extract_chains(self, G: nx.DiGraph, valid_nodes: list[int]) -> list[list]:
        """
        Restrict the graph to *valid_nodes* and return the weakly-connected
        components as topologically-sorted lists (upstream -> downstream).

        Each component is a linear chain because all nodes have in/out-degree 1
        in the full graph, so within the subgraph the structure is a simple path.
        """
        sub = G.subgraph(valid_nodes)
        chains = []
        for component in nx.weakly_connected_components(sub):
            ordered = list(nx.topological_sort(sub.subgraph(component)))
            chains.append(ordered)
        return chains

    def _filter_chains(self, G: nx.DiGraph, chains: list[list]) -> list[list]:
        """
        Apply per-chain validity checks controlled by the params.

        Checks
        ------
        RequireIdenticalOrbits      - all reaches share the same SWOT orbit set
        DrainageAreaPctCutoff       - no consecutive pair exceeds the pct threshold
        MaximumReachesEachDirection - chain trimmed to this radius around midpoint
        MinimumReaches              - chains shorter than the floor are discarded
        """
        max_radius = self.params.get("MaximumReachesEachDirection", 999)
        min_reaches = self.params.get("MinimumReaches", 1)
        da_cutoff = self.params.get("DrainageAreaPctCutoff", float("inf"))
        req_orbits = self.params.get("RequireIdenticalOrbits", False)

        valid = []
        for chain in chains:
            # Trim to MaximumReachesEachDirection around the midpoint
            if len(chain) > 2 * max_radius + 1:
                mid = len(chain) // 2
                chain = chain[max(0, mid - max_radius) : mid + max_radius + 1]

            if len(chain) < min_reaches:
                continue
            if req_orbits and not has_consistent_orbits(G, chain):
                continue
            if has_area_jumps(G, chain, da_cutoff):
                continue

            valid.append(chain)

        return valid

    def _deduplicate(self, chains: list[list]) -> list[list]:
        """
        Remove duplicate and high-overlap chains.

        Two chains are near-duplicates when their overlap fraction exceeds
        AllowedReachOverlap.  Longer chains are preferred when there is a
        conflict.
        """
        overlap_threshold = self.params.get("AllowedReachOverlap", 1.0)

        # Longest-first so conflicts keep the more informative set
        chains_sorted = sorted(chains, key=len, reverse=True)

        kept: list[frozenset] = []
        kept_ordered: list[list] = []

        for chain in chains_sorted:
            fs = frozenset(chain)
            if overlap_threshold == -1:
                if fs in kept:
                    continue
            elif any(calc_overlap_frac(fs, k) > overlap_threshold for k in kept):
                continue
            kept.append(fs)
            kept_ordered.append(chain)

        return kept_ordered

    def _add_singleton_sets(self, chains: list[list], G: nx.DiGraph) -> list[list]:
        """
        For every input river reach that does not appear in any chain, add a
        length-1 chain containing just that reach.  Only applied when
        MinimumReaches == 1.
        """
        covered = {r for chain in chains for r in chain}
        for rid in self.reach_ids:
            if rid not in covered and is_river_reach(rid) and rid in G:
                chains.append([rid])
        return chains

    def _chains_to_inversion_sets(self, chains: list[list], G: nx.DiGraph) -> dict:
        """
        Produce the same {reach_id: {ReachList, numReaches, Reaches}}
        dictionary that the original getsets returned.
        """
        inversion_sets = {}
        for chain in chains:
            origin = chain[len(chain) // 2]
            reaches_dict = {
                rid: {**dict(G.nodes[rid]), "reach_id": rid}
                for rid in chain
                if rid in G
            }
            inversion_sets[origin] = {
                "ReachList": chain,
                "numReaches": len(chain),
                "Reaches": reaches_dict,
            }
        return inversion_sets

    # def print_stats(self, inversion_sets: dict) -> None:
    #     sizes = [v["numReaches"] for v in inversion_sets.values()]
    #     covered = {r for v in inversion_sets.values() for r in v["ReachList"]}

    #     print(f"    Total input reaches      : {len(self.reach_ids)}")
    #     print(f"    Inversion sets found     : {len(inversion_sets)}")
    #     print(f"    Total reaches in sets    : {sum(sizes)}")
    #     print(f"    Input reaches covered    : {len(covered & self.reach_ids)}")
    #     if sizes:
    #         print(f"    Mean set size            : {np.mean(sizes):.1f}")
    #         print(f"    Max set size             : {max(sizes)}")
    #     else:
    #         print(f"    Mean set size            : n/a (no sets found)")
    #         print(f"    Max set size             : n/a (no sets found)")

    def print_stats(self,InversionSets):
        # output some stats
        numReaches=[]
        for IS in InversionSets:
            numReaches.append(InversionSets[IS]['numReaches'])
            
        reaches_in_sets=[]
        for IS in InversionSets:
            for reachid in InversionSets[IS]['ReachList']:
                reaches_in_sets.append(reachid)        
                
        reaches_specified=[]
        for reach in self.reaches:
            reaches_specified.append(reach['reach_id'])  

        reaches_specified = [r["reach_id"] for r in self.reaches]                                              
        
        numOverlap=len( list( set(reaches_in_sets) & set(reaches_specified) ) )                            
        
        print('    total number of reaches:',len(self.reaches))
        print('    A total of', len(InversionSets.keys()),'sets were identified.')
        print('    Total reaches included in sets:',sum(numReaches))
        print('    Of the ',len(self.reaches),' reaches input to SetFinder, ',numOverlap,' are included in the sets')
        

    def write_inversion_set_data(
        self,
        inversion_sets: dict,
        output_dir: str,
        expanded: bool = False,
    ) -> list[int]:
        """Serialize inversion_sets to a JSON file and return all reach IDs."""
        prefix = "expanded_" if expanded else ""
        out_path = os.path.join(output_dir, prefix + self.params["Filename"])

        # Extract required file paths from the reaches attribute
        swordfile = self.reaches[0]['sword']
        sosfile = self.reaches[0]['sos']
        
        # Utilize helper method for data structure generation
        InversionSetsWrite, all_reaches = self.get_IS_list(inversion_sets, swordfile, sosfile, 'writing')

        with open(out_path, 'w') as json_file:
            json.dump(InversionSetsWrite, json_file, indent=2)
            print(f"    File written: {out_path}")

        return all_reaches
    
    def get_IS_list(self, InversionSets, swordfile, sosfile, mode):
        """
        Constructs the nested list structure for inversion set serialization.
        Restores the 'remove_dupes' mode functionality and origin reach tracking.
        """
        InversionSetsList = []
        all_reaches = []

        for IS in InversionSets:
            InversionSetWrite = []
            for reach in InversionSets[IS]['ReachList']:
                all_reaches.append(reach)
                reachdict = {
                    'reach_id': int(reach),
                    'sword': swordfile,
                    'swot': f"{reach}_SWOT.nc",
                    'sos': sosfile
                }
                
                # Restore conditional data inclusion based on mode
                if mode == 'remove_dupes':
                    reachdict['origin'] = InversionSets[IS]['OriginReach']['reach_id']
                    
                InversionSetWrite.append(reachdict)
            InversionSetsList.append(InversionSetWrite)

        return InversionSetsList, all_reaches


    def getsets(self) -> dict:
        print("Building graph from SWORD...")
        G = self._build_graph()

        print("Identifying valid chain nodes...")
        if self.params.get("AllowRiverJunction", False):
            # Relax the degree-1 constraint — every river reach is a candidate
            candidates = (
                self.reach_ids
                if self.params.get("RequireSetReachesInput")
                else set(G.nodes)
            )
            valid_nodes = [n for n in candidates if n in G and is_river_reach(n)]
        else:
            valid_nodes = self._valid_river_nodes(G)
        print(len(valid_nodes))

        print("Extracting chains...")
        chains = self._extract_chains(G, valid_nodes)
        print(len(chains))

        print("Filtering chains...")
        chains = self._filter_chains(G, chains)
        print(len(chains))

        print("Deduplicating...")
        chains = self._deduplicate(chains)
        print(len(chains))

        if self.params.get("MinimumReaches", 1) == 1:
            print("Adding singleton sets for uncovered reaches...")
            chains = self._add_singleton_sets(chains, G)

        print("Building output structure...")
        inversion_sets = self._chains_to_inversion_sets(chains, G)

        print("Stats:")
        self.print_stats(inversion_sets)

        return inversion_sets
