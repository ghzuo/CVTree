/*
 * Copyright (c) 2024
 * Wenzhou Institute, University of Chinese Academy of Sciences.
 * See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 * 
 * @Author: Dr. Guanghong Zuo
 * @Date: 2024-05-11 14:47:45
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2024-05-11 19:07:13
 */

#ifndef SAMPLING_H
#define SAMPLING_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>

#include "sampleMeth.h"
#include "genome.h"

void sampleGenome(SampleMeth*, const string&, const vector<string>&);

#endif // SAMPLING_H

