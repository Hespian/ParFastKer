 /******************************************************************************
 * modified.h
 *
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

#ifndef MODIFIED_H
#define MODIFIED_H

#include <vector>
#include <cassert>

class branch_and_reduce_algorithm;

class modified {

public:
    int add;
    std::vector<int> removed;
    std::vector<int> vs;
    std::vector<std::vector<int>> oldAdj;
    branch_and_reduce_algorithm *pAlg;

public:
    modified(int const add, std::vector<int> &removed, std::vector<int> &vs, std::vector<std::vector<int>> &newAdj, branch_and_reduce_algorithm *_pAlg);

    modified(std::vector<int> &removed, std::vector<int> &vs, branch_and_reduce_algorithm *_pAlg);

    virtual ~modified() {};

    void restore();

    virtual void reverse(std::vector<int> &x) = 0;
};

class fold : public modified
{

public:
    fold(int const add, std::vector<int> &removed, std::vector<int> &vs, std::vector<std::vector<int>> &newAdj, branch_and_reduce_algorithm *_pAlg)
    : modified(add, removed, vs, newAdj, _pAlg)
    { }

    fold(std::vector<int> &removed, std::vector<int> &vs, branch_and_reduce_algorithm *_pAlg)
    : modified(removed, vs, _pAlg)
    { }

    virtual ~fold() {}

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

class alternative : public modified {

public:
    int k;

    alternative(int const add, std::vector<int> &removed, std::vector<int> &vs, std::vector<std::vector<int>> &newAdj, branch_and_reduce_algorithm *_pAlg, int k)
    : modified(add, removed, vs, newAdj, _pAlg)
    {
        this->k = k;
    }

    alternative(std::vector<int> &removed, std::vector<int> &vs, branch_and_reduce_algorithm *_pAlg, int k)
    : modified(removed, vs, _pAlg)
    {
        this->k = k;
    }

    virtual ~alternative() {}

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

#endif // MODIFIED_H
