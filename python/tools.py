#!/usr/bin/env python

def parse_enum(cpp_file, enum_name):
    """Poor man's C++ enum parsing."""
    import os
    src_text = open(os.path.join(os.getenv('CMSSW_BASE'), cpp_file)).read()

    n = src_text.index(enum_name)
    b = src_text.index('{', n) + 1
    e = src_text.index('}', n)
    src_text = strip_comments(src_text[b:e])
    
    enum = {}
    i = -1

    items = []
    for item in src_text.split(','):
        item = item.strip()
        if not item or item.startswith('//'):
            continue
        if '=' in item:
            item, value = item.split('=')
            i = eval(value.strip())
        else:
            i += 1
        enum[item.strip()] = i
    return enum
    
def rec_level_code(rec):
    """Returns a unique one-character code for the given rec level. Currently maps 0-9 onto '0'-'9' and 10- to 'A'-."""
    if rec < 10: return chr(rec + 48)
    else: return chr(rec + 55)

def rec_levels():
    """Return the dict obtained from parsing the rec level enum from RecLevelHelper.h, and a list of the levels in order."""
    level_dict = parse_enum('src/SUSYBSMAnalysis/Zprime2muAnalysis/src/RecLevelHelper.h', 'RecLevel')
    levels = [pair[0][1:] for pair in sorted(level_dict.items(), key=lambda pair: pair[1])][:-1] # don't include lLAST in this list
    return level_dict, levels
    
def replace_all(s, start, end, rep=''):
    n = s.find(start)
    while n >= 0:
        s = s[:n] + rep + s[(s.find(end, n)+1):]
        n = s.find(start)
    return s

def strip_comments(cpp_src):
    """Poor man's C++ comment parsing."""
    cpp_src = replace_all(cpp_src, '//', '\n', '\n')
    cpp_src = replace_all(cpp_src, '/*', '*/')
    return cpp_src

# Useful backports from python version > 2.4.
try:
    from collections import defaultdict
except:
    # http://code.activestate.com/recipes/523034/
    class defaultdict(dict):
        def __init__(self, default_factory=None, *a, **kw):
            if (default_factory is not None and
                not hasattr(default_factory, '__call__')):
                raise TypeError('first argument must be callable')
            dict.__init__(self, *a, **kw)
            self.default_factory = default_factory
        def __getitem__(self, key):
            try:
                return dict.__getitem__(self, key)
            except KeyError:
                return self.__missing__(key)
        def __missing__(self, key):
            if self.default_factory is None:
                raise KeyError(key)
            self[key] = value = self.default_factory()
            return value
        def __reduce__(self):
            if self.default_factory is None:
                args = tuple()
            else:
                args = self.default_factory,
            return type(self), args, None, None, self.items()
        def copy(self):
            return self.__copy__()
        def __copy__(self):
            return type(self)(self.default_factory, self)
        def __deepcopy__(self, memo):
            import copy
            return type(self)(self.default_factory,
                              copy.deepcopy(self.items()))
        def __repr__(self):
            return 'defaultdict(%s, %s)' % (self.default_factory,
                                            dict.__repr__(self))

try:
    from itertools import product
except ImportError:
    # http://docs.python.org/library/itertools.html#itertools.product
    def product(*args, **kwds):
        pools = map(tuple, args) * kwds.get('repeat', 1)
        result = [[]]
        for pool in pools:
            result = [x+[y] for x in result for y in pool]
        for prod in result:
            yield tuple(prod)
    
__all__ = ['defaultdict', 'parse_enum', 'product', 'rec_level_code', 'rec_levels', 'replace_all', 'strip_comments']

if __name__ == '__main__':
    print 'test replace_all:'
    print replace_all("""
    line 1; // a comment
    line 2;
    line 3  a sdg     g;
    """, '//', '\n', '\n')
    print 'parse_enum test', parse_enum('src/SUSYBSMAnalysis/Zprime2muAnalysis/src/RecLevelHelper.h', 'RecLevel')
    print 'parse_enum test', parse_enum('src/SUSYBSMAnalysis/Zprime2muAnalysis/src/CutHelper.h',      'CutResult')
   
