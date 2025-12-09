import numpy as np
from scipy.optimize import linear_sum_assignment
from collections import defaultdict


class LabeledGraph:
    def __init__(self, graph, labels):

        self.graph = graph
        self.labels = labels
        self.label2idxs = defaultdict(list)
        self.props = {}
        self.build_label2idxs()

        self.ini_labels = self.labels.copy()
        self.ini_label2idxs = self.label2idxs.copy()
        self.ini_wl_labels = self.get_WL_labels()

        self.irred_labels = {}

    def build_label2idxs(self):
        self.label2idxs.clear()
        for idx, l in enumerate(self.labels):
            self.label2idxs[l].append(idx)

    def copy(self):
        obj = type(self).__new__(type(self))
        obj.graph = {k: v.copy() for k, v in self.graph.items()}
        obj.labels = self.labels.copy()
        obj.label2idxs = defaultdict(
            list, {k: v.copy() for k, v in self.label2idxs.items()}
        )
        obj.props = self.props

        obj.ini_labels = self.ini_labels
        obj.ini_label2idxs = self.ini_label2idxs
        obj.ini_wl_labels = self.ini_wl_labels

        obj.irred_labels = self.irred_labels.copy()

        return obj

    def set_prop(self,name,prop):
        self.props[name] = prop

    def binarize_graph(self):
        for i in self.graph:
            for j in self.graph[i]:
                self.graph[i][j] = 1

    def get_WL_labels(self):

        wl_labels = self.labels.copy()

        for _ in range(2 * len(wl_labels)):
            temp_wl_labels = []
            for i in range(len(wl_labels)):
                nbrs = defaultdict(list)
                if i in self.graph:
                    for j, w in self.graph[i].items():
                        nbrs[wl_labels[j]].append(w)

                hash_input = []
                for k, ws in nbrs.items():
                    ws.sort(reverse=True)
                    hash_input.append((k, tuple(ws)))

                next_label = abs(hash((wl_labels[i], tuple(sorted(hash_input)))))
                temp_wl_labels.append(next_label)

            if len(set(temp_wl_labels)) == len(set(wl_labels)):
                return wl_labels
            else:
                wl_labels = temp_wl_labels


class SlapMapper:
    INF_INT = 100000

    def __init__(self, binary=False):

        self.stack = []
        self.visited = []
        self.results = []
        self.minval = SlapMapper.INF_INT

        self.binary = binary
        self.valfactor = 2

    def reset(self):

        self.stack.clear()
        self.visited.clear()
        self.results.clear()
        self.minval = SlapMapper.INF_INT

    def solve(self):
        while self.stack:
            lgp = self.stack.pop()
            self._solve_all_lap(lgp)

    def remove_non_min_results(self):

        min_val = min(r["val"] for r in self.results)
        self.results = [r for r in self.results if r["val"] == min_val]

    def remove_isomorphic_results(self, ini_mlg, break_sym_targets):

        if len(self.results) < 2:
            return

        visited = []
        new_results = []
        for r in self.results:
            nshift = len(r["lgp"][0].labels)
            mlg = ini_mlg.copy()
            for idx0 in break_sym_targets:
                l = r["lgp"][0].labels[idx0]
                idx1 = r["lgp"][1].label2idxs[l][0]
                mlg.graph[idx0][idx1 + nshift] = 1
                mlg.graph[idx1 + nshift][idx0] = 1

            key = hash(tuple(sorted(mlg.get_WL_labels())))
            if key not in visited:
                visited.append(key)
                new_results.append(r)

        self.results = new_results

    def get_maps(self, lgp, break_sym_targets=None, interactive=False, base=None):

        self.reset()

        if break_sym_targets is None:
            interactive = False

        ini_lgp = [lgp[0].copy(), lgp[1].copy()]
        temp_lgps = [[lgp[0].copy(), lgp[1].copy()]]
        temp_results = []

        if self.binary:
            for temp_lgp in temp_lgps:
                for lg in temp_lgp:
                    lg.binarize_graph()

        ini_mlg = merge_lgp(temp_lgps[0])

        if interactive:
            given = [[], []]
            while base not in [0, 1]:
                while True:
                    base = input("select 0 or 1-based indexing [0/1]:")
                    try:
                        base = int(base)
                        break
                    except ValueError:
                        pass

        while temp_lgps:
            self.stack.extend(temp_lgps)
            self.results.clear()
            self.minval = SlapMapper.INF_INT
            self.solve()
            self.remove_non_min_results()

            if break_sym_targets is None:
                break

            if interactive:
                temp_lgps, given = self._break_sym_interactive(
                    ini_lgp, base, given, break_sym_targets
                )
                self.visited.clear()
            else:
                temp_lgps, temp_results = self._break_sym(
                    temp_results, break_sym_targets
                )

        if not interactive and break_sym_targets is not None:
            self.results = temp_results

        self.remove_non_min_results()
        if break_sym_targets is not None:
            self.remove_isomorphic_results(ini_mlg, break_sym_targets)

        for r in self.results:
            simplify_labels(r['lgp'])

            if r['val']%self.valfactor==0:
                r['cd'] = int(r['val']//self.valfactor)
            else:
                r['cd'] = r['val']/self.valfactor

            if interactive:
                r["base"] = base
                r["label_string"] = ";".join(
                    [f"{i+base}>>{j+base}" for i, j in zip(given[0], given[1])]
                )

    def _break_sym(self, temp_results, break_sym_targets):
        temp_lgps = []
        while self.results:
            r = self.results.pop()
            lgp = r["lgp"]

            sorted_labels = self._get_sorted_labels(lgp)
            l = None
            for temp_l in reversed(sorted_labels):
                if (
                    len(lgp[0].label2idxs[temp_l]) > 1
                    and lgp[0].label2idxs[temp_l][0] in break_sym_targets
                ):
                    l = temp_l
                    break


            if l is None:
                temp_results.append(r)
                continue

            for new_l in range(1000, 2000):
                if new_l not in lgp[0].label2idxs.keys():
                    break

            sols, next_info_pair = r["lap_sols"][l]
            for sol in sols:
                if len(set(sol["groups_pair"][0])) > 1:
                    continue

                info00 = list(next_info_pair[0].values())[0]
                idx0 = min(info00["idxs"])

                if sol["fingerprint"] is None:
                    j_nonzero = 0
                else:
                    j_nonzero = np.where(sol["fingerprint"][0] > 0)[0][0]
                info1j = list(next_info_pair[1].values())[j_nonzero]

                wl_labels1 = lgp[1].ini_wl_labels
                visited = []
                for idx1 in info1j["idxs"]:

                    if wl_labels1[idx1] not in visited:
                        visited.append(wl_labels1[idx1])

                        new_lgp = [lgp[0].copy(), lgp[1].copy()]
                        new_lgp[0].labels[idx0] = new_l
                        new_lgp[0].label2idxs[l].remove(idx0)
                        new_lgp[0].label2idxs[new_l].append(idx0)
                        new_lgp[1].labels[idx1] = new_l
                        new_lgp[1].label2idxs[l].remove(idx1)
                        new_lgp[1].label2idxs[new_l].append(idx1)

                        key = hash_current_labels(new_lgp)
                        if key not in self.visited:
                            self._del_irred_labels(new_lgp,[l,new_l])
                            temp_lgps.append(new_lgp)
                            self.visited.append(key)

        return temp_lgps, temp_results


    def _break_sym_interactive(self, ini_lgp, base, given, break_sym_targets):
        idx2idxs = defaultdict(set)
        for r in self.results:
            lgp = r["lgp"]
            for l, idxs in lgp[0].label2idxs.items():
                if idxs[0] not in break_sym_targets:
                    continue
                sols, next_info_pair = r["lap_sols"][l]
                info0 = list(next_info_pair[0].values())
                ninfo0 = len(next_info_pair[0])
                info1 = list(next_info_pair[1].values())
                ninfo1 = len(next_info_pair[1])
                if ninfo0 == 1 or ninfo1 == 1:
                    fng = np.ones([ninfo0, ninfo1], dtype=bool)
                else:
                    fng = np.zeros([ninfo0, ninfo1], dtype=bool)
                    for sol in sols:
                        if len(set(sol["groups_pair"][0])) > 1:
                            continue
                        fng = fng | sol["fingerprint"] > 0
                list_idxs0 = [set(info["idxs"]) for info in info0]
                list_idxs1 = [set(info["idxs"]) for info in info1]
                for i in range(ninfo0):
                    for j in range(ninfo1):
                        if fng[i, j]:
                            for idx0 in list_idxs0[i]:
                                idx2idxs[idx0].update(list_idxs1[j])

        min_item = None
        min_len = SlapMapper.INF_INT
        for idx0 in break_sym_targets:
            idxs1 = idx2idxs[idx0]
            if len(idxs1) > 1 and len(idxs1) < min_len:
                min_item = (idx0, idxs1)
                min_len = len(idxs1)

        temp_lgps = []
        if min_item is not None:
            idx0, idxs1 = min_item
            str_idxs1 = ", ".join(str(i + base) for i in sorted(idxs1))

            while True:
                idx1 = input(f"{idx0+base} >> ? (probably in {str_idxs1}):")
                try:
                    idx1 = int(idx1) - base
                    break
                except ValueError:
                    pass

            given[0].append(idx0)
            given[1].append(idx1)

            new_lgp = [ini_lgp[0].copy(), ini_lgp[1].copy()]
            new_l = 1000
            for idx0, idx1 in zip(given[0], given[1]):
                for _ in range(SlapMapper.INF_INT):
                    new_l += 1
                    if new_l not in new_lgp[0].label2idxs.keys():
                        break
                new_lgp[0].labels[idx0] = new_l
                new_lgp[1].labels[idx1] = new_l
            new_lgp[0].build_label2idxs()
            new_lgp[1].build_label2idxs()

            temp_lgps.append(new_lgp)

        return temp_lgps, given

    def _solve_all_lap(self, lgp, avoid_multi_sols=[1]):

        sorted_labels = self._get_sorted_labels(lgp)

        tot_val = 0
        lap_sols = dict()
        for l in sorted_labels:

            is_irred = False
            while True:

                if l in lgp[0].irred_labels:
                    sols, next_info_pair = lgp[0].irred_labels[l]
                    is_irred = True
                else:
                    sols, next_info_pair = self._solve_lap_with_label(lgp, l)

                tot_val += sols[0]["val"]
                if tot_val > self.minval:
                    return

                if is_irred:
                    break

                ini_label = lgp[0].ini_labels[lgp[0].label2idxs[l][0]]
                if ini_label in avoid_multi_sols:
                    cnt = 0
                    for sol in sols:
                        if sol["proper"]:
                            cnt += 1
                    if cnt>1:
                        is_irred = True

                if is_irred:
                    break

                for sol in sols:
                    if sol["proper"]:
                        if len(set(sol["groups_pair"][0])) == 1:
                            is_irred = True
                            break
                        else:
                            new_lgp = self._update_labels(
                                lgp, l, sol, next_info_pair
                            )
                            key = hash_current_labels(new_lgp)
                            if key not in self.visited:
                                self.stack.append(new_lgp)
                                self.visited.append(key)

                break

            lap_sols[l] = (sols, next_info_pair)

            if is_irred:
                if l not in lgp[0].irred_labels:
                    lgp[0].irred_labels[l] = (sols, next_info_pair)
                continue
            else:
                return

        self.results.append({"lgp": lgp, "val": tot_val, "lap_sols": lap_sols})
        self.minval = min(tot_val, self.minval)

    def _diff(self, x, y):

        if self.binary:
            return abs(len(x) - len(y))

        l = min(len(x), len(y))
        d = 0

        for i in range(l):
            d += abs(x[i] - y[i])

        if len(x) > l:
            d += sum(x[l:])
        elif len(y) > l:
            d += sum(y[l:])

        return d

    def _diff_nbrs(self, nbrs_x, nbrs_y):

        all_keys = set(nbrs_x.keys()) | set(nbrs_y.keys())

        d = 0
        for key in all_keys:
            d += self._diff(nbrs_x.get(key, []), nbrs_y.get(key, []))

        return d

    def _get_next_info(self, lg, label):

        next_info = {}
        for i in lg.label2idxs[label]:
            nbrs = defaultdict(list)
            if i in lg.graph:
                for j, w in lg.graph[i].items():
                    nbrs[lg.labels[j]].append(w)

            hash_input = []
            for k, ws in nbrs.items():
                ws.sort(reverse=True)
                hash_input.append((k, tuple(ws)))

            next_label = abs(hash((label, tuple(sorted(hash_input)))))

            if next_label not in next_info:
                next_info[next_label] = {"idxs": [i], "nbrs": nbrs}
            else:
                next_info[next_label]["idxs"].append(i)

        return next_info

    def _get_cost_matrix(self, lgp, label):

        next_info_pair = []
        slices_pair = []
        mask_mat_pair = []
        for lg in lgp:
            next_info = self._get_next_info(lg, label)
            nlabels = len(next_info)
            slices = [0]
            for l in next_info.keys():
                slices.append(slices[-1] + len(next_info[l]["idxs"]))

            mask_mat = np.zeros([nlabels, slices[-1]], dtype=int)
            for i in range(nlabels):
                mask_mat[i, slices[i] : slices[i + 1]] = 1

            next_info_pair.append(next_info)
            slices_pair.append(slices)
            mask_mat_pair.append(mask_mat)

        cost_matrix = np.empty([slices_pair[0][-1], slices_pair[1][-1]], dtype=int)

        for s00, s01, v0 in zip(
            slices_pair[0][:-1], slices_pair[0][1:], next_info_pair[0].values()
        ):
            for s10, s11, v1 in zip(
                slices_pair[1][:-1], slices_pair[1][1:], next_info_pair[1].values()
            ):
                cost_matrix[s00:s01, s10:s11] = self._diff_nbrs(v0["nbrs"], v1["nbrs"])

        return cost_matrix, mask_mat_pair, next_info_pair

    def _solve_lap(self, cost, mask_mat_pair):

        nlabels = []
        nlabels.append(len(mask_mat_pair[0]))
        nlabels.append(len(mask_mat_pair[1]))

        sols = []

        if nlabels[0] == 1 or nlabels[1] == 1:
            row = np.arange(len(cost), dtype=int)
            col = np.arange(len(cost), dtype=int)
            val = np.sum(cost[row, col])

            sol = {}
            sol["row"] = row
            sol["col"] = col
            sol["val"] = val
            sol["fingerprint"] = None
            sol["groups_pair"] = [[0] * nlabels[0], [0] * nlabels[1]]

            sols.append(sol)

        else:
            min_val = 1e10
            pert_cost = 100 * cost
            for cnt in range(10000):

                row, col = linear_sum_assignment(pert_cost)
                val = np.sum(cost[row, col])

                if cnt == 0:
                    min_val = val

                if min_val != val:
                    break

                fingerprint = mask_mat_pair[0][:, row] @ (mask_mat_pair[1][:, col].T)

                is_obtained = False
                for s in sols:
                    is_obtained = np.all(fingerprint == s["fingerprint"])
                    if is_obtained:
                        break
                if is_obtained:
                    break

                grps0 = [-1] * nlabels[0]
                grps1 = [-1] * nlabels[1]
                cnt_g = 0
                for i in range(nlabels[0]):
                    if grps0[i] != -1:
                        continue

                    set0 = set()
                    set1 = set()
                    next_set0 = {i}

                    for _ in range(nlabels[0]):
                        diff0 = next_set0 - set0
                        if diff0:
                            next_set1 = set()
                            for i0 in diff0:
                                grps0[i0] = cnt_g
                                next_set1.update(np.where(fingerprint[i0] > 0)[0])
                            set0.update(diff0)
                        else:
                            cnt_g += 1
                            break

                        diff1 = next_set1 - set1
                        if diff1:
                            next_set0 = set()
                            for i1 in diff1:
                                grps1[i1] = cnt_g
                                next_set0.update(np.where(fingerprint.T[i1] > 0)[0])
                            set1.update(diff1)
                        else:
                            cnt_g += 1
                            break

                sol = {}
                sol["row"] = row
                sol["col"] = col
                sol["val"] = val
                sol["fingerprint"] = fingerprint
                sol["groups_pair"] = [grps0, grps1]

                sols.append(sol)

                pert_cost += (
                    (mask_mat_pair[0].T @ fingerprint @ mask_mat_pair[1]) > 0
                ).astype(int)

        grps = []
        for sol in sols:
            sol["proper"] = True
            grps.append(sol["groups_pair"][0] + sol["groups_pair"][1])

        nsols = len(sols)
        nls = len(grps[0])

        for ig in range(nsols):
            for jg in range(nsols):
                if ig == jg or not sols[jg]["proper"]:
                    continue
                d = defaultdict(set)
                for idx in range(nls):
                    d[grps[ig][idx]].add(grps[jg][idx])
                if all(len(s) == 1 for s in d.values()):
                    sols[ig]["proper"] = False
                    break

        return sols

    def _solve_lap_with_label(self, lgp, label):
        cost, mask_mat_pair, next_info_pair = self._get_cost_matrix(
            lgp, label
        )
        sols = self._solve_lap(cost, mask_mat_pair)

        return sols, next_info_pair

    def _get_sorted_labels(self, lgp):

        l2i = lgp[0].label2idxs
        sorted_labels = sorted(l2i.keys(), key=lambda k: -len(l2i[k]))

        return sorted_labels

    def _update_labels(self, lgp, label, sol, next_info_pair):

        new_lgp = [lgp[0].copy(), lgp[1].copy()]

        g2l = {}
        for g, l in zip(sol["groups_pair"][1], next_info_pair[1].keys()):
            g2l[g] = l

        for i in range(2):
            merged_next_idxs = defaultdict(list)
            for g, l in zip(sol["groups_pair"][i], next_info_pair[i].keys()):
                merged_next_idxs[g2l[g]].extend(next_info_pair[i][l]["idxs"])

            del new_lgp[i].label2idxs[label]
            for l, idxs in merged_next_idxs.items():
                new_lgp[i].label2idxs[l] = sorted(idxs)
                for idx in idxs:
                    new_lgp[i].labels[idx] = l

        self._del_irred_labels(new_lgp,g2l.values())

        return new_lgp


    def _del_irred_labels(self,lgp,parents):

        neighs = defaultdict(set)
        for l in parents:
            for i in range(2):
                for idx in lgp[i].label2idxs[l]:
                    neighs[l].update(lgp[i].labels[idx_next]
                        for idx_next in lgp[i].graph[idx].keys())

        influenced = set(parents)
        for k0,s0 in neighs.items():
            for k1,s1 in neighs.items():
                if k0>k1:
                    influenced.update(s0&s1)

        for l in influenced:
            if l in lgp[0].irred_labels:
                del lgp[0].irred_labels[l]


def merge_lgp(lgp):
    nshift = len(lgp[0].labels)
    graph0 = {i: {j: w for j, w in d.items()} for i, d in lgp[0].graph.items()}
    shifted_graph1 = {
        i + nshift: {j + nshift: w for j, w in d.items()}
        for i, d in lgp[1].graph.items()
    }

    return LabeledGraph({**graph0, **shifted_graph1}, lgp[0].labels + lgp[1].labels)


def connect_mlg(lgp, ini_mlg):

    nshift = len(lgp[0].labels)
    mlg = ini_mlg.copy()

    for l, idxs0 in lgp[0].label2idxs.items():
        idxs1 = lgp[1].label2idxs[l]
        for i0 in idxs0:
            for i1 in idxs1:
                mlg.graph[i0][i1 + nshift] = 1
                mlg.graph[i1 + nshift][i0] = 1

    return mlg


def hash_mlg(lgp, ini_mlg):
    mlg = connect_mlg(lgp, ini_mlg)
    return hash(tuple(sorted(mlg.get_WL_labels())))


def hash_lgp(lgp):
    mlg = LabeledGraph({}, lgp[0].ini_labels + lgp[1].ini_labels)
    nshift = len(lgp[0].labels)
    for i,v in lgp[0].graph.items():
        mlg.graph[i] = {j:w for j,w in v.items()}
    for i,v in lgp[1].graph.items():
        mlg.graph[i+nshift] = {j+nshift:w for j,w in v.items()}
    for l, idxs0 in lgp[0].label2idxs.items():
        idxs1 = lgp[1].label2idxs[l]
        for i0 in idxs0:
            for i1 in idxs1:
                mlg.graph[i0][i1 + nshift] = 1
                mlg.graph[i1 + nshift][i0] = 1
    return hash(tuple(sorted(mlg.get_WL_labels())))


def simplify_labels(lgp):
    old2new = {}
    cnt = 1
    for l in lgp[0].labels:
        if l not in old2new:
            old2new[l] = cnt
            cnt += 1

    for lg in lgp:
        for i, l in enumerate(lg.labels):
            lg.labels[i] = old2new[l]
        lg.build_label2idxs()


def hash_current_labels(lgp):
    old2new = {}
    cnt = 1
    for l in lgp[0].labels:
        if l not in old2new:
            old2new[l] = cnt
            cnt += 1

    labels = []
    for lg in lgp:
        for i, l in enumerate(lg.labels):
            labels.append(old2new[l])

    return hash(tuple(labels))

