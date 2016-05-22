 /******************************************************************************
 * Copyright (C) 2015-2017 Darren Strash <strash@kit.edu>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 2 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#include "parallel_branch_and_reduce_algorithm.h"
// #include "mis_config.h"
#include "fast_set.h"
#include "parallel_modified.h"
#include "data_structure/graph_access.h"
#include "data_structure/priority_queues/maxNodeHeap.h"

#include <vector>
#include <set>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstdio>
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <algorithm>  // sort()
#include <deque>
#include <chrono>

#include "kaHIP_interface.h"
////#define debug(x) {fprintf(2, x);}

using namespace std;

int  parallel_branch_and_reduce_algorithm::REDUCTION   = 3;
int  parallel_branch_and_reduce_algorithm::LOWER_BOUND = 4;
int  parallel_branch_and_reduce_algorithm::BRANCHING   = 2;
bool parallel_branch_and_reduce_algorithm::outputLP    = false;
long parallel_branch_and_reduce_algorithm::nBranchings = 0;
int  parallel_branch_and_reduce_algorithm::debug       = 0;

parallel_branch_and_reduce_algorithm::parallel_branch_and_reduce_algorithm(vector<vector<int>> &_adj, int const _N, MISConfig mis_config)
: adj() 
, n(_adj.size())
, mis_config(mis_config)
{
////    srand(4327897);
    SHRINK = 0.5;
    depth = 0;
    maxDepth = 10;
    rootDepth = -1; // invalid value

    n = _adj.size();
    adj.swap(_adj);

    N = _N;
    y.resize(N, 0);
    for (int i = 0; i < n; i++) y[i] = 1;
    for (int i = n; i < N; i++) y[i] = 2;
    x.resize(N, 0);
    for (int i = 0; i < n; i++) x[i] = -1;
    for (int i = n; i < N; i++) x[i] = 2;
    rn = n;
    in.resize(n, -1);
    out.resize(n, -1);

    level.resize(mis_config.number_of_partitions);
    for(int i = 0; i < level.size(); ++i) {
        level[i].resize(n * 2, 0);
    }

    modTmp.resize(mis_config.number_of_partitions);
    for(int i = 0; i < modTmp.size(); ++i) {
        modTmp[i].resize(n, 0);
    }

    modifiedN.resize(mis_config.number_of_partitions, 0);
    modifieds.resize(mis_config.number_of_partitions);
    for(int i = 0; i < modifieds.size(); ++i) {
        modifieds[i].resize(N, shared_ptr<parallel_modified>());
    }

    partitions.resize(_N);
    partition_nodes.resize(mis_config.number_of_partitions);

    for(int i = 0; i < mis_config.number_of_partitions; ++i) {
        used.push_back(fast_set(n * 2));
    }
}

int parallel_branch_and_reduce_algorithm::deg(int v) {
    assert(x[v] < 0);
    int deg = 0;
    for (int u : adj[v]) if (x[u] < 0) deg++;
    return deg;
}

void parallel_branch_and_reduce_algorithm::set(int v, int a)
{
    assert(x[v] < 0);
    x[v] = a;
    if (a == 0) {
        for (int u : adj[v]) if (x[u] < 0) {
            x[u] = 1;
        }
    }
}

// methods that modify the graph

void parallel_branch_and_reduce_algorithm::compute_fold(vector<int> const &S, vector<int> const &NS, int partition) {
    assert(NS.size() == S.size() + 1);
    vector<int> removed(S.size() * 2);
    for (unsigned int i = 0; i < S.size(); i++) removed[i] = S[i];
    for (unsigned int i = 0; i < S.size(); i++) removed[S.size() + i] = NS[1 + i];
    int s = NS[0];
    fast_set &used_partition = used[partition];
    used_partition.clear();
    for (int v : S) used_partition.add(v);
    vector<int> &tmp = modTmp[partition];
    int p = 0;
    for (int v : NS) {
        assert(!used_partition.get(v));
        for (int u : adj[v]) if (x[u] < 0 && used_partition.add(u)) {
            tmp[p++] = u;
        }
    }
    vector<vector<int>> newAdj(p + 1);
    {
        vector<int> copyOfTmp(tmp.begin(), tmp.begin() + p);
        newAdj[0].swap(copyOfTmp);
    }
    std::sort(newAdj[0].begin(), newAdj[0].end());
    vector<int> vs(p + 1);
    vs[0] = s;
    used_partition.clear();
    for (int v : S) used_partition.add(v);
    for (int v : NS) used_partition.add(v);
    for (unsigned int i = 0; i < newAdj[0].size(); i++) {
        int v = newAdj[0][i];
        p = 0;
        bool add = false;
        for (int u : adj[v]) if (x[u] < 0 && !used_partition.get(u)) {
            if (!add && s < u) {
                tmp[p++] = s;
                add = true;
            }
            tmp[p++] = u;
        }
        if (!add) tmp[p++] = s;
        vs[1 + i] = v;

        {
            vector<int> copyOfTmp(tmp.begin(), tmp.begin() + p);
            newAdj[i+1].swap(copyOfTmp);
        }
    }
    modifieds[partition][modifiedN[partition]++] = make_shared<parallel_fold>(parallel_fold(S.size(), removed, vs, newAdj, this));
////    cout << __LINE__ << ", " << this << ", " << depth << ": Setting modifieds[" << modifiedN-1 << "]=" << modifieds[modifiedN-1] << endl << flush;
}

void parallel_branch_and_reduce_algorithm::reverse() {
    for (int partition = 0; partition < mis_config.number_of_partitions; ++partition) {
        for (int i = modifiedN[partition] - 1; i >= 0; i--) {
                modifieds[partition][i]->reverse(y);
            }
    }
}

bool parallel_branch_and_reduce_algorithm::fold2Reduction() {
    std::vector<char> changed_per_partition(mis_config.number_of_partitions, 0);
    #pragma omp parallel for
    for(int partition = 0; partition < mis_config.number_of_partitions; ++partition) {
        vector<int> &tmp = level[partition];
        for (int v : partition_nodes[partition]) if (x[v] < 0) {
            int p = 0;
            for (int u : adj[v]) {
                if (partitions[u] != partition) {
                    goto loop;
                }
                if (x[u] < 0) {
                tmp[p++] = u;
                if (p > 2) goto loop;
                }
            }
            if (p < 2) continue;
            for(int i = 0; i < 2; ++i) {
                for(int u : adj[tmp[i]]) {
                    if(partitions[u] != partition) {
                        goto loop;
                    }
                }
            }
            for (int u : adj[tmp[0]]) if (u == tmp[1]) {
                set(v, 0);
                goto loop;
            }
            {
            vector<int> copyOfTmp(tmp.begin(), tmp.begin() + 2);
            compute_fold(vector<int>{v}, copyOfTmp, partition);
            changed_per_partition[partition] = 1;
            }
loop:       ;
        }
    }
    
    for(char changed_partition : changed_per_partition) {
        if(changed_partition != 0) {
            return true;
        }
    }
    return false;
}

bool parallel_branch_and_reduce_algorithm::isolatedCliqueReduction() {
    std::vector<char> changed_per_partition(mis_config.number_of_partitions, 0);

    #pragma omp parallel for
    for(int partition = 0; partition < mis_config.number_of_partitions; ++partition) {
        for(int v : partition_nodes[partition]) {
            if(x[v] < 0) {
                bool changed_node = isolatedCliqueReduction(v, partition);
                if(changed_node) {
                    changed_per_partition[partition] = 1;
                }
            }
        }
    }


    for(char changed : changed_per_partition) {
        if(changed != 0) {
            return true;
        }
    }
    return false;
}

bool parallel_branch_and_reduce_algorithm::isolatedCliqueReduction(NodeID vertex, int partition) {
    auto degreeVertex = deg(vertex);
    for (int neighbor : adj[vertex]) {
        if (partitions[neighbor] != partition) {
            return false;
        }
        if(x[neighbor] < 0 && deg(neighbor) < degreeVertex) {
            return false;
        }
    }
    fast_set &used_partition = used[partition];
    for (int neighbor : adj[vertex]) {
        if(x[neighbor] < 0) {
            used_partition.clear();

            for (int nNeighbor : adj[neighbor]) {
                if(x[nNeighbor] < 0) {
                    used_partition.add(nNeighbor);
                }
            }
            used_partition.add(neighbor);

            for (int neighbor2 : adj[vertex]) {
                if(x[neighbor2] < 0 && !used_partition.get(neighbor2)) {
                    return false;
                }
            }
        }
    }
    set(vertex, 0);
    return true;
}

std::string parallel_branch_and_reduce_algorithm::debugString() const {
    stringstream ins;
#ifdef PUT_TIME
    time_t rawtime;
    struct tm *timeinfo;

    time (&rawtime);
    timeinfo = localtime (&rawtime);
    ins << std::put_time(timeinfo, "%T") << "  ";
#else
    std::locale::global(std::locale("ja_JP.utf8"));
    std::time_t t = std::time(NULL);
    char mbstr[100];
    if (std::strftime(mbstr, sizeof(mbstr), "%T", std::localtime(&t))) {
        std::cout << mbstr << '\n';
    }
#endif
    for (int i = 0; i < depth && i < maxDepth; ++i) {
        ins << " ";
    }
////    printf ("Current local time and date: %s", asctime(timeinfo));
    return ins.str();

}

void parallel_branch_and_reduce_algorithm::PrintState() const
{
    cout << "State(" << this << "):" << endl << flush;
    cout << "adj=" << endl << flush;
    for (unsigned int j = 0; j < adj.size(); ++j) {
        cout << j << " : ";
        for (int const k : adj[j]) {
            cout << k << " ";
        }
        cout << endl;
    }
    cout << "N  =" << N << endl << flush;
    cout << "in =";
    for (int const i : in) {
        cout << i << " ";
    }
    cout << endl << flush;
    cout << "out=";
    for (int const i : out) {
        cout << i << " ";
    }
    cout << endl << flush;
}


void parallel_branch_and_reduce_algorithm::reduce_graph()
{
    std::vector<int> xadj;
    std::vector<int> adjncy;
    for(auto node_adj_list : adj) {
        xadj.push_back(adjncy.size());
        for(auto neighbor : node_adj_list) {
            adjncy.push_back(neighbor);
        }
    }
    xadj.push_back(adjncy.size());
    int edgecut = 0;
    int number_of_partitions = (int)mis_config.number_of_partitions;
    double imbalance = mis_config.imbalance;
    kaffpa(&n, NULL, xadj.data(), NULL, adjncy.data(), &number_of_partitions, &imbalance, false, rand(), mis_config.kahip_mode, &edgecut, partitions.data());
    std::cout << "Edgecut: " << edgecut << std::endl;
    for(NodeID node = 0; node < N; ++node) {
        partition_nodes[partitions[node]].push_back(node);
    }
    clock_t begin = clock();
    for (;;) {
        // if (REDUCTION >= 0) deg1Reduction();
////        if (n > 100 && n * SHRINK >= rn && !outputLP && decompose()) return true;
        // if (REDUCTION >= 0 && REDUCTION < 2 && dominateReduction()) continue;
        // if (REDUCTION >= 2 && unconfinedReduction()) continue;
        // if (REDUCTION >= 1 && lpReduction()) continue;
        /*if (REDUCTION >= 3) {
            int r = packingReduction();
            if (r < 0) return;
            if (r > 0) continue;
        }*/
        if (REDUCTION >= 1 && fold2Reduction()) continue;
        if (REDUCTION >= 1 && isolatedCliqueReduction()) continue;
        // if (REDUCTION >= 2 && twinReduction()) continue;
        // if (REDUCTION >= 2 && funnelReduction()) continue;
        // if (REDUCTION >= 2 && deskReduction()) continue;
        break;
    }
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cout << "Parallel took " << elapsed_secs << " seconds" << endl;
    size_t low_degree_count(0);
    for (int v = 0; v < n; v++) if (x[v] < 0) {
        if (deg(v) <= 1) {
            low_degree_count++;
        }
    }
    cout << "There are " << low_degree_count << " degree 0 and 1 vertices left!" << endl << flush;
}

size_t parallel_branch_and_reduce_algorithm::get_current_is_size() const {

    vector<int> x2(x);

    for(int partition = 0; partition < mis_config.number_of_partitions; ++partition) {
        for (int i = modifiedN[partition] - 1; i >= 0; i--) {
            modifieds[partition][i]->reverse(x2);
        }
    }

    size_t current_is_size(0);
    for (unsigned int i = 0; i < adj.size(); ++i) {
        if (x2[i] == 0) {
            current_is_size++;
        }
    }

    return current_is_size;
}

size_t parallel_branch_and_reduce_algorithm::get_current_is_size_with_folds() const {

    size_t folded_vertex_count(0);
    size_t current_is_size(0);
    for (int const i : x) {
        if (i == 0) current_is_size++;
        if (i == 2) folded_vertex_count++;
    }

    return current_is_size + folded_vertex_count/2;
}

bool parallel_branch_and_reduce_algorithm::folded_vertices_exist() const {

    vector<int> x2(x);

    for(int partition = 0; partition < mis_config.number_of_partitions; ++partition) {
        for (int i = modifiedN[partition] - 1; i >= 0; i--) {
            modifieds[partition][i]->reverse(x2);
        }
    }

    for (int const i : x2) {
        if (i == 2) return true;
    }

    return false;
}

vector<int> parallel_branch_and_reduce_algorithm::compute_maximal_is() {

    int vertexToForceInIndependentSet(0);
    while (vertexToForceInIndependentSet != -1) {
        reduce_graph();

        vertexToForceInIndependentSet = -1;
        for (unsigned int i = 0; i < x.size(); ++i) {
            if (x[i] == -1) { // status not determined
                vertexToForceInIndependentSet = i;
                break;
            }
        }

        // add vertex to independent set
        if (vertexToForceInIndependentSet != -1) {
            set(vertexToForceInIndependentSet, 0);
        }
    }

    vector<int> x2(x);

    for(int partition = 0; partition < mis_config.number_of_partitions; ++partition) {
        for (int i = modifiedN[partition] - 1; i >= 0; i--) {
            modifieds[partition][i]->reverse(x2);
        }
    }

    size_t current_is_size(0);
    for (unsigned int i = 0; i < adj.size(); ++i) {
        if (x2[i] == 0) {
            current_is_size++;
        }
    }

    return x2;
}

size_t parallel_branch_and_reduce_algorithm::compute_alternative_maximal_is_size() {

    int vertexToForceInIndependentSet(0);
    while (vertexToForceInIndependentSet != -1) {
        reduce_graph();

        vertexToForceInIndependentSet = -1;
        for (int i = 0; i < static_cast<int>(x.size()); ++i) {
            if (x[i] == -1) { // status not determined
                vertexToForceInIndependentSet = i;
                break;
            }
        }

        // add vertex to independent set
        if (vertexToForceInIndependentSet != -1) {
            set(vertexToForceInIndependentSet, 0);
        }
    }

    size_t numberOfFoldedVertices(0);
    size_t sizeOfIS(0);
    for (int const i : x) {
        if (i == 0) sizeOfIS++;
        if (i == 2) numberOfFoldedVertices++;
    }

    return sizeOfIS + numberOfFoldedVertices/2;
}

size_t parallel_branch_and_reduce_algorithm::number_of_nodes_remaining() const {

    size_t node_count(0);
    for (int i : x) if (i == -1) node_count++;

    return node_count;
}

void parallel_branch_and_reduce_algorithm::force_into_independent_set(vector<NodeID> const &nodes) {

    for (NodeID const node : nodes) {
        assert(x[node] == -1); // should not have been assigned yet.
        if (x[node] != -1) {
            cout << "ERROR: invalid vertex selected for independent set!" << endl << flush;
        }
////        cout << "Switching node " << node << " from value " << x[node] << " to 0" << endl;
        set(node, 0);
    }
}

