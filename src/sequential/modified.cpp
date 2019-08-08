 /******************************************************************************
 * Copyright (C) 2015-2017 Darren Strash <strash@kit.edu>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *****************************************************************************/

// local includes
#include "modified.h"
#include "branch_and_reduce_algorithm.h"

// system includes
#include <vector>

modified::modified(int const _add, std::vector<int> &_removed, std::vector<int> &_vs, std::vector<std::vector<int>> &newAdj, branch_and_reduce_algorithm *_pAlg)
: add(_add)
, pAlg(_pAlg) 
{
    removed.swap(_removed);
    vs.swap(_vs);
    oldAdj.resize(vs.size());
    pAlg->crt += add;
    for (int i = 0; i < static_cast<int>(removed.size()); i++) pAlg->vRestore[--(pAlg->rn)] = -1;
    for (int v : removed) {
        assert(pAlg->x[v] < 0);
        pAlg->x[v] = 2;
    }
    for (int i = 0; i < static_cast<int>(vs.size()); i++) {
        oldAdj[i].swap(pAlg->adj[vs[i]]);
        pAlg->adj[vs[i]].swap(newAdj[i]);
    }
}

modified::modified(std::vector<int> &_removed, std::vector<int> &_vs, branch_and_reduce_algorithm *_pAlg)
: add(0)
, pAlg(_pAlg)
{
    removed.swap(_removed);
    vs.swap(_vs);
}

void modified::restore() {
    pAlg->crt -= add;
    pAlg->rn += removed.size();
    for (int v : removed) pAlg->x[v] = -1;
    for (int i = 0; i < static_cast<int>(vs.size()); i++) {
        pAlg->adj[vs[i]] = oldAdj[i];
        int inV = pAlg->in[vs[i]], outV = pAlg->out[vs[i]];
        for (int u : pAlg->adj[vs[i]]) {
            if (u == inV) inV = -1;
            if (u == outV) outV = -1;
        }
        if (inV >= 0) {
            pAlg->out[pAlg->in[vs[i]]] = -1;
            pAlg->in[vs[i]] = -1;
        }
        if (outV >= 0) {
            pAlg->in[pAlg->out[vs[i]]] = -1;
            pAlg->out[vs[i]] = -1;
        }
    }
}


