import glob, os, re
from collections import defaultdict

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

__all__ = [
    'dataset_from_publish_log',
    'files_from_crab_dir',
    ]
