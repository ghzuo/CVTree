/*
 * Copyright (c) 2025
 * See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2025-12-14 Sunday 13:59:59
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2025-12-14 Sunday 22:58:35
 */

#ifndef ALIGN_H
#define ALIGN_H

#include <algorithm>
#include <functional>
#include <stack>
#include <vector>
using namespace std;

// the sequence align function
template <typename T, typename R>
void seqAlign(const vector<T> &seq1, const vector<T> &seq2, R &res,
              function<void(const T &, const T &, R &)> act) {
  typename vector<T>::const_iterator it1 = seq1.begin();
  typename vector<T>::const_iterator it2 = seq2.begin();
  for (;;) {
    if (*it1 < *it2) {
      if (++it1 == seq1.end())
        break;
    } else if (*it2 < *it1) {
      if (++it2 == seq2.end())
        break;
    } else {
      act(*it1, *it2, 0);
      if (++it1 == seq1.end() || ++it2 == seq2.end())
        break;
    }
  }
};

// the binary align function

// 处理序列对的辅助函数模板
template <typename T, typename R>
void processPair(typename vector<T>::const_iterator shortBeg,
                 typename vector<T>::const_iterator shortEnd,
                 typename vector<T>::const_iterator longBeg,
                 typename vector<T>::const_iterator longEnd,
                 stack<pair<pair<typename vector<T>::const_iterator,
                                 typename vector<T>::const_iterator>,
                            pair<typename vector<T>::const_iterator,
                                 typename vector<T>::const_iterator>>> &stack,
                 function<void(const T &, const T &, R &)> act, R &res) {

  // 从较短序列的中间取出一个元素, 并在较长序列中寻找不小于该元素的位置
  auto mid = shortBeg + (shortEnd - shortBeg) / 2;
  auto lowerBound = lower_bound(longBeg, longEnd, *mid);

  // 将序列分为两段，分别对齐；如果没有找到则将较长序列剩余部分与较短前部对齐
  // 如果找到相等的元素，执行act操作, 并将相等的元素从较长序列中移除
  // 检查不要将空序列压入栈
  if (lowerBound == longEnd && shortEnd != mid) {
    stack.push({{shortBeg, mid}, {longBeg, longEnd}});
  } else if (*mid < *lowerBound) {
    if (mid != shortBeg && longBeg != lowerBound)
      stack.push({{shortBeg, mid}, {longBeg, lowerBound}});
    if (mid + 1 != shortEnd)
      stack.push({{mid + 1, shortEnd}, {lowerBound, longEnd}});
  } else {
    act(*mid, *lowerBound, res);
    if (mid != shortBeg && longBeg != lowerBound)
      stack.push({{shortBeg, mid}, {longBeg, lowerBound}});
    if (mid + 1 != shortEnd && lowerBound + 1 != longEnd)
      stack.push({{mid + 1, shortEnd}, {lowerBound + 1, longEnd}});
  }
};

template <typename T, typename R>
void binaryAlign(const vector<T> &seq1, const vector<T> &seq2, R &res,
                 function<void(const T &, const T &, R &)> act) {

  // check the input not empty
  if (seq1.empty() || seq2.empty())
    return;
  
  // Initial the stack
  using Iter = typename vector<T>::const_iterator;
  stack<pair<pair<Iter, Iter>, pair<Iter, Iter>>> stack;
  stack.push({{seq1.begin(), seq1.end()}, {seq2.begin(), seq2.end()}});

  while (!stack.empty()) {
    // 从栈中取出一对序列
    auto seqs = stack.top();
    Iter beg1 = seqs.first.first;
    Iter end1 = seqs.first.second;
    Iter beg2 = seqs.second.first;
    Iter end2 = seqs.second.second;
    stack.pop();

    if (end1 - beg1 < end2 - beg2) {
      processPair<T, R>(beg1, end1, beg2, end2, stack, act, res);
    } else {
      processPair<T, R>(beg2, end2, beg1, end1, stack, act, res);
    }
  }
};

#endif