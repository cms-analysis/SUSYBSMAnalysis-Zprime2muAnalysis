#!/usr/bin/env python

import os, subprocess

def hadd_ex(new_name, files):
    """Use ROOT's hadd tool to merge files into a new file with path
    new_name. This is a simple wrapper that suppresses the stdout from
    hadd, only reporting a summary line of how many files were
    merged. We check that the number of files reported merged by hadd
    is the same as the number in the input list, or if there were any
    other problems reported by hadd. If so, we print an error to
    stdout (with reverse video to make it stand out) and return
    False. On success, return True.
    """
    
    l = len(files)
    print 'hadding %i files to %s' % (l, new_name)
    args = ['hadd', new_name] + files

    p = subprocess.Popen(args=args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    stdout, stderr = p.communicate()
    open(new_name + '.haddlog', 'wt').write(stdout)
    assert stderr is None

    if p.returncode != 0:
        print '\033[36;7m PROBLEM hadding \033[m', new_name
        #print p.stdout.read()
        return False

    max_file_num = 0
    for line in stdout.split('\n'):
        if 'Source file' in line:
            max_file_num = max(max_file_num, int(line.split(':')[0].split(' ')[-1]))
    print max_file_num, 'files merged to', new_name
    if max_file_num != l:
        print '\033[36;7m PROBLEM hadding', new_name
        return False

    return True

def hadd(new_name, files, chunk_size=900):
    """Use hadd_ex above to merge files into new_name, but in chunks
    to get around hadd's limit of 999 input files. Returns whether
    operation was a success.
    """
    
    if len(files) <= chunk_size:
        return hadd_ex(new_name, files)
    
    files = files[:]
    new_files = []
    while files:
        these = files[:chunk_size]
        files = files[chunk_size:]

        this_fn = new_name + '_%i' % len(new_files)
        new_files.append(this_fn)

        if not hadd_ex(this_fn, these):
            print '\033[36;7m PROBLEM hadding \033[m', new_name, 'in chunks of', chunk_size, 'on', this_fn
            return False

    assert len(new_files) < chunk_size

    ok = hadd_ex(new_name, new_files)
    if not ok:
        print '\033[36;7m PROBLEM hadding', new_name, 'in chunks of', chunk_size, 'assembling final file'
        return False

    for fn in new_files:
        os.remove(fn)

    return True

__all__ = [
    'hadd',
    'hadd_ex',
    ]

if __name__ == '__main__':
    import sys
    hadd(sys.argv[1], sys.argv[2:])
