#!/usr/bin/env python

import os, cPickle, gzip

def big_warn(s):
    x = '#' * len(s)
    print x
    print x
    print x
    print s
    print x
    print x
    print x

def files_from_dbs(dataset, ana02=True):
    # could use DBSAPI but this is easier
    url = '--url https://cmsdbsprod.cern.ch:8443/cms_dbs_ph_analysis_02_writer/servlet/DBSServlet' if ana02 else ''
    cmd = 'dbs search %(url)s --query "find file where dataset=%(dataset)s"' % locals()
    cmdout = os.popen(cmd).readlines()
    ret = [y.strip('\n') for y in cmdout if '.root' in y]
    if not ret:
        raise RuntimeError('no files for %s (ana02: %s) found. dbs command output:\n' % (dataset, ana02) + ''.join(cmdout))
    return ret

def parse_enum(cpp_file, enum_name, as_list=False, drop_last=False, remove_substr=None):
    """Poor man's C++ enum parsing."""
    import os
    src_text = open(os.path.join(os.getenv('CMSSW_BASE'), cpp_file)).read()

    n = src_text.index(enum_name)
    b = src_text.index('{', n) + 1
    e = src_text.index('}', n)
    src_text = strip_comments(src_text[b:e])
    
    enum = {}
    i = -1

    last = None
    
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
            
        item = item.strip()
        if remove_substr is not None:
            item = item.replace(remove_substr, '')

        if last is None or i >= last[1]:
            last = (item, i)

        enum[item] = i

    if drop_last:
        del enum[last[0]]

    if as_list:
        enum = enum.items()
        enum.sort(key=lambda x: x[1])
        enum = [x[0] for x in enum]

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

def to_pickle(obj, fn, proto=-1, comp=False):
    if comp or '.gzpickle' in fn:
        f = gzip.GzipFile(fn, 'wb')
    else:
        f = open(fn, 'wb')
    cPickle.dump(obj, f, proto)

__all__ = ['big_warn', 'files_from_dbs', 'parse_enum', 'rec_level_code', 'rec_levels', 'replace_all', 'strip_comments']

if __name__ == '__main__':
    print 'test replace_all:'
    print replace_all("""
    line 1; // a comment
    line 2;
    line 3  a sdg     g;
    """, '//', '\n', '\n')
    print 'parse_enum test', parse_enum('src/SUSYBSMAnalysis/Zprime2muAnalysis/src/RecLevelHelper.h', 'RecLevel')
    print 'parse_enum test', parse_enum('src/SUSYBSMAnalysis/Zprime2muAnalysis/src/CutHelper.h',      'CutResult')
   
