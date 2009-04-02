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

__all__ = ['parse_enum', 'rec_level_code', 'rec_levels', 'replace_all', 'strip_comments']

if __name__ == '__main__':
    print 'test replace_all:'
    print replace_all("""
    line 1; // a comment
    line 2;
    line 3  a sdg     g;
    """, '//', '\n', '\n')
    print 'parse_enum test', parse_enum('src/SUSYBSMAnalysis/Zprime2muAnalysis/src/RecLevelHelper.h', 'RecLevel')
    print 'parse_enum test', parse_enum('src/SUSYBSMAnalysis/Zprime2muAnalysis/src/CutHelper.h',      'CutResult')
   
