/******************************************************************************
 * macros_common.h               
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 ******************************************************************************
 * Copyright (C) 2013-2015 Christian Schulz <christian.schulz@kit.edu>
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

#ifndef MACROS_COMMON_H 
#define MACROS_COMMON_H

// If DEBUG has not been defined yet then define it as false.
#ifndef DEBUG
# define DEBUG 0
#endif

#if DEBUG
# define VLOG(x) x
# define VVLOG(x) x
#else
# define VLOG(x)
# define VVLOG(x)
#endif


// A macro to define empty copy constructors and assignment operators.
// Use this after private: to hide them to the outside.
#define DISALLOW_COPY_AND_ASSIGN(TypeName) \
  TypeName(const TypeName&);               \
  void operator=(const TypeName&)


// A macro to disallow the assignment operators.
#define DISALLOW_ASSIGN(TypeName) \
  void operator=(const TypeName&)


// Helper macro for STR.
#define ASSERT_H_XSTR(x) (#x)

// This macro allows to convert an expression to a string.
#define STR(x) ASSERT_H_XSTR(x)


#endif // ifndef MACROS_COMMON_H

