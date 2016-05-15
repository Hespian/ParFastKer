 /******************************************************************************
 * parallel_modified.h
 *
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

#ifndef PARALLEL_MODIFIED_H
#define PARALLEL_MODIFIED_H

#include <vector>
#include <cassert>

class parallel_branch_and_reduce_algorithm;

class parallel_modified {

public:
    int add;
    std::vector<int> removed;
    std::vector<int> vs;
    std::vector<std::vector<int>> oldAdj;
    parallel_branch_and_reduce_algorithm *pAlg;

public:
    parallel_modified(int const add, std::vector<int> &removed, std::vector<int> &vs, std::vector<std::vector<int>> &newAdj, parallel_branch_and_reduce_algorithm *_pAlg);

    parallel_modified(std::vector<int> &removed, std::vector<int> &vs, parallel_branch_and_reduce_algorithm *_pAlg);

    virtual ~parallel_modified() {};

    void restore();

    virtual void reverse(std::vector<int> &x) = 0;
};

class parallel_fold : public parallel_modified
{

public:
    parallel_fold(int const add, std::vector<int> &removed, std::vector<int> &vs, std::vector<std::vector<int>> &newAdj, parallel_branch_and_reduce_algorithm *_pAlg)
    : parallel_modified(add, removed, vs, newAdj, _pAlg)
    { }

    parallel_fold(std::vector<int> &removed, std::vector<int> &vs, parallel_branch_and_reduce_algorithm *_pAlg)
    : parallel_modified(removed, vs, _pAlg)
    { }

    virtual ~parallel_fold() {}

    virtual void reverse(std::vector<int> &x) {
        int k = removed.size() / 2;
        if (x[vs[0]] == 0) {
            for (int i = 0; i < k; i++) x[removed[i]] = 1;
            for (int i = 0; i < k; i++) x[removed[k + i]] = 0;
        } else if (x[vs[0]] == 1) {
            for (int i = 0; i < k; i++) x[removed[i]] = 0;
            for (int i = 0; i < k; i++) x[removed[k + i]] = 1;
        }
    }
};

class parallel_alternative : public parallel_modified {

public:
    int k;

    parallel_alternative(int const add, std::vector<int> &removed, std::vector<int> &vs, std::vector<std::vector<int>> &newAdj, parallel_branch_and_reduce_algorithm *_pAlg, int k)
    : parallel_modified(add, removed, vs, newAdj, _pAlg)
    {
        this->k = k;
    }

    parallel_alternative(std::vector<int> &removed, std::vector<int> &vs, parallel_branch_and_reduce_algorithm *_pAlg, int k)
    : parallel_modified(removed, vs, _pAlg)
    {
        this->k = k;
    }

    virtual ~parallel_alternative() {}

    void reverse(std::vector<int> &x) {
        bool A0 = false, A1 = true;
        bool B0 = false, B1 = true;
        for (int i = 0; i < k; i++) {
            if (x[vs[i]] == 0) A0 = true;
            if (x[vs[i]] != 1) A1 = false;
        }
        for (int i = k; i < static_cast<int>(vs.size()); i++) {
            if (x[vs[i]] == 0) B0 = true;
            if (x[vs[i]] != 1) B1 = false;
        }
        if (A1 || B0) {
            for (int i = 0; i < static_cast<int>(removed.size()/ 2); i++) x[removed[i]] = 0;
            for (int i = static_cast<int>(removed.size() / 2); i < static_cast<int>(removed.size()); i++) x[removed[i]] = 1;
        } else if (B1 || A0) {
            for (int i = 0; i < static_cast<int>(removed.size() / 2); i++) x[removed[i]] = 1;
            for (int i = static_cast<int>(removed.size() / 2); i < static_cast<int>(removed.size()); i++) x[removed[i]] = 0;
        }
    }
};

#endif // parallel_MODIFIED_H
