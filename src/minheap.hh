#ifndef MINHEAP_HH_SHYJ79GD
#define MINHEAP_HH_SHYJ79GD

#include <vector>
#include <queue>

using std::vector;
using std::priority_queue;

template<typename T>
class MinHeap {
public:

    MinHeap(size_t max_size) : _maxsz(max_size) {}

    inline bool push(T value)
    {
        if (_q.size() < _maxsz) {
            _q.push(value);
        } else if (value < _q.top()) {
            _q.pop();
            _q.push(value);
        } else {
            return false;
        }
        return true;
    }

    vector<T> finalize() {
      std::vector<T> result(_q.size());
      while (_q.size()) {
        result[_q.size() - 1] = _q.top();
        _q.pop();
      }
      return result;
    }

    inline T top() { return _q.top(); }
    inline T pop() { return _q.pop(); }
    inline size_t size() { return _q.size(); }
    inline bool full() { return _q.size() == _maxsz; }
    const priority_queue<T> & queue() { return _q; }

private:
    size_t _maxsz;
    priority_queue<T> _q;
};

#if 0
template<typename T>
class MinHeap {
public:

    MinHeap(size_t max_size) : _maxsz(max_size) {}

    inline bool push(T value)
    {
        if (_q.size() < _maxsz) {
            _q.push(value);
        } else if (value < _q.top()) {
            _q.pop();
            _q.push(value);
        } else {
            return false;
        }
        return true;
    }

    vector<T> finalize() {
      std::vector<T> result(_q.size());
      while (_q.size()) {
        result[_q.size() - 1] = _q.top();
        _q.pop();
      }
      return result;
    }

    inline T top() { return _q.top(); }
    inline T pop() { return _q.pop(); }
    inline size_t size() { return _q.size(); }
    inline bool full() { return _q.size() == _maxsz; }
    const priority_queue<T> & queue() { return _q; }

private:
    size_t _maxsz;
    priority_queue<T> _q;
};
#endif


#endif /* end of include guard: MINHEAP_HH_SHYJ79GD */
