import glob, os, re
from collections import defaultdict

def crab_status(working_dir, debug=False):
    d = defaultdict(list)
    
    cmd = 'crab -c %s -status' % working_dir
    if debug: print cmd
    s = os.popen4(cmd)[1].read()
    
    if 'Total Jobs' not in s:
        raise RuntimeError('unable to get status for working_dir=' + working_dir)
    
    for x in s.split('\n'):
        x = [y.strip() for y in x.split(' ') if y.strip()]
        if len(x) < 2: continue

        try:
            x[0] = int(x[0])
        except ValueError:
            pass
        else:
            id, status = x[:2]
            if len(x) > 3:
                codes = x[3:]
                
            if status == 'Retrieved':
                key = '%s_%s' % (status, '_'.join(codes))
            else:
                key = status
            d[key].append(id)

    if debug:
        print s
        for k in sorted(d.keys()):
            print '%s: %s' % (k.ljust(25), crabify_list(d[k]))
            
    return d

def crabify_list(l):
    return ','.join(str(x) for x in sorted(l))

def files_from_crab_dir(crab_dir):
    fjrs = glob.glob(os.path.join(crab_dir, 'res', 'crab_fjr*xml'))
    fjrs.sort(key = lambda x: int(x.split('_')[-1].split('.xml')[0]))
    files = []

    # Fragile xml parsing!
    wrapper_re = re.compile(r'<FrameworkError ExitStatus="(.+)" Type="WrapperExitCode"/>')
    exe_re = re.compile(r'<FrameworkError ExitStatus="(.+)" Type="ExeExitCode"/>')
    filename_re = re.compile(r'[ \t](/store/user.*root)')
    for fjr in fjrs:
        s = open(fjr).read()
        wrapper_mo = wrapper_re.search(s)
        exe_mo = exe_re.search(s)
        filename_mo = filename_re.search(s)
        if wrapper_mo is None or exe_mo is None or filename_mo is None:
            raise RuntimeError('cannot parse %s for exit codes and output filename (wrapper_mo %s exe_mo %s filename_mo %s)' % (fjr, wrapper_mo, exe_mo, filename_mo))
        if wrapper_mo.group(1) != '0' or exe_mo.group(1) != '0':
            raise RuntimeError('exit codes for %s not 0 (wrapper %s, exe %s)' % (fjr, wrapper_mo.group(1), exe_mo.group(1)))
        files.append(filename_mo.group(1))

    return files

def dataset_from_publish_log(publish_log_fn):
    ad = [x.strip().replace('=== dataset ', '') for x in open(publish_log_fn).readlines() if x.startswith('=== dataset') and not x.startswith('=== dataset description')]
    assert(len(ad) == 1)
    return ad[0]

def last_crab_dir(crab_dir_base='crab', filter=''):
    return sorted([x for x in glob.glob(os.path.join(crab_dir_base, '*')) if filter in x and os.path.isdir(x) and os.path.isdir(os.path.join(x, 'res'))], key = lambda x: os.stat(x).st_ctime)[-1]

__all__ = [
    'crab_status',
    'crabify_list',
    'dataset_from_publish_log',
    'files_from_crab_dir',
    'last_crab_dir'
    ]
