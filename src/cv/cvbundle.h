/*
 * Copyright (c) 2025
 * See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 * 
 * @Author: Dr. Guanghong Zuo
 * @Date: 2025-12-06 Saturday 16:24:17
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2025-12-10 Wednesday 14:01:19
 */

#ifndef CVBUNDLE_H
#define CVBUNDLE_H

#include <vector>

#include "cv.h"

typedef vector<CVvec> CVBundle;

string readCVBundle(const string&, CVBundle&);
void writeCVBundle(const string&, const CVBundle&);

#endif // CVBUNDLE_H